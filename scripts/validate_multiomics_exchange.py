#!/usr/bin/env python3
"""Strict, offline validation of a MultiOmics downstream exchange v1 package."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def validate(root):
    root = Path(root).resolve(strict=True)
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    errors = []
    if manifest.get("schema_version") != "1.0.0":
        errors.append("schema_version must be 1.0.0")
    files = manifest.get("files", {})
    matrices = manifest.get("matrices", [])
    if not matrices:
        errors.append("At least one matrix is required")

    def read_verified(spec):
        if not isinstance(spec, dict):
            raise ValueError("File manifest entry is not an object")
        rel = spec.get("path", "")
        if not rel or Path(rel).is_absolute() or ".." in Path(rel).parts:
            raise ValueError("Unsafe relative path: " + repr(rel))
        path = (root / rel).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError("Missing or escaped file: " + str(rel))
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        if digest != spec.get("sha256"):
            raise ValueError("SHA-256 mismatch: " + str(rel))
        return path

    def read_csv(spec):
        path = read_verified(spec)
        with path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            header = reader.fieldnames or []
            if len(header) != len(set(header)):
                raise ValueError("Duplicate CSV headers: " + str(path))
            rows = list(reader)
            if any(None in row or any(value is None for value in row.values()) for row in rows):
                raise ValueError("Ragged CSV rows: " + str(path))
            return header, rows

    try:
        sh, samples = read_csv(files["samples"])
        fh, features = read_csv(files["features"])
        for field in ("assay_id", "sample_id", "subject_id", "omic"):
            if field not in sh:
                errors.append("samples.csv missing " + field)
        for field in ("omic", "feature_id", "id_namespace"):
            if field not in fh:
                errors.append("features.csv missing " + field)
        if errors:
            return errors
        specimen_owners = {}
        specimen_context = {}
        for sample in samples:
            sid, subject = sample["sample_id"], sample["subject_id"]
            if sid in specimen_owners and specimen_owners[sid] != subject:
                errors.append("One sample_id assigned to multiple biological subjects")
            specimen_owners[sid] = subject
            previous = specimen_context.get(sid, {})
            for field in ('tissue', 'collection_time', 'collection_time_unit', 'timepoint', 'timepoint_unit', 'time_origin'):
                current_value = sample.get(field, '') or ''
                previous_value = previous.get(field, '')
                if current_value and previous_value and current_value != previous_value:
                    errors.append('Conflicting ' + field + ' for shared specimen ' + sid)
            specimen_context[sid] = {
                field: sample.get(field, '') or previous.get(field, '')
                for field in ('tissue', 'collection_time', 'collection_time_unit', 'timepoint', 'timepoint_unit', 'time_origin')
            }
        # Optional individual signatures are never inferred from group q-values.
        # Validate measured assay identity and declared score/model provenance.
        # SHA syntax is not proof of an independently reproduced frozen fit.
        if "individual_signatures" in files:
            _, signatures = read_csv(files["individual_signatures"])
            mandatory = ("subject_id", "sample_id", "assay_id", "omic",
                         "signature_id", "score", "score_scale", "fit_scope",
                         "reference_model_sha256", "feature_definition_sha256",
                         "source_exchange_analysis_id", "score_provenance")
            if not signatures:
                errors.append("Empty optional individual signatures file")
            by_assay = {(s["omic"], s["assay_id"]): s for s in samples}
            duplicates = set()
            for score in signatures:
                if any(not score.get(field, "").strip() for field in mandatory):
                    errors.append("Incomplete individual signature row")
                    continue
                assay = by_assay.get((score["omic"], score["assay_id"]))
                if not assay or score["subject_id"] != assay["subject_id"] or score["sample_id"] != assay["sample_id"]:
                    errors.append("Signature assay/specimen/subject mismatch")
                try:
                    if not math.isfinite(float(score["score"])):
                        errors.append("Nonfinite individual signature")
                except (ValueError, TypeError):
                    errors.append("Nonnumeric individual signature")
                if score["score_scale"].lower() in ("raw_counts", "group_logfc", "log2foldchange", "p_value", "q_value", "unknown"):
                    errors.append("Group or raw counts cannot be declared individual signature scale")
                if score["fit_scope"] not in ("training_frozen_test_projection", "out_of_fold_training"):
                    errors.append("Invalid signature fitting/projection scope")
                if score["score_provenance"] not in ("native_r_individual_score", "external_frozen_projection"):
                    errors.append("Missing eligible signature provenance")
                if score["source_exchange_analysis_id"] != manifest.get("analysis_id"):
                    errors.append("Signature source analysis ID mismatch")
                for key in ("reference_model_sha256", "feature_definition_sha256"):
                    digest = score[key]
                    if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                        errors.append("Signature has invalid model/feature SHA-256")
                unique = (score["omic"], score["assay_id"], score["signature_id"])
                if unique in duplicates:
                    errors.append("Duplicate assay/signature")
                duplicates.add(unique)
        assay_keys = [(r["omic"], r["assay_id"]) for r in samples]
        feature_keys = [(r["omic"], r["feature_id"]) for r in features]
        if len(set(assay_keys)) != len(assay_keys):
            errors.append("Duplicate (omic, assay_id) in samples")
        if len(set(feature_keys)) != len(feature_keys):
            errors.append("Duplicate (omic, feature_id) in features")
        if any(not all(r.get(x, "").strip() for x in ("assay_id", "sample_id", "subject_id", "omic")) for r in samples):
            errors.append("Blank mandatory sample identifier")
        if any(not all(r.get(x, "").strip() for x in ("omic", "feature_id", "id_namespace")) for r in features):
            errors.append("Blank mandatory feature identifier")
        matrix_omics = [m.get('omic') for m in matrices]
        if set(r['omic'] for r in samples) != set(matrix_omics):
            errors.append('Declared sample omics do not match available matrices')
        seen_omics = set()
        for mat in matrices:
            omic = mat.get("omic")
            if omic in seen_omics:
                errors.append("Duplicate matrix for omic " + str(omic))
            seen_omics.add(omic)
            if mat.get("orientation") != "feature_rows_sample_columns":
                errors.append("Invalid matrix orientation")
            if mat.get("missing_value") != "NA":
                errors.append("Matrix missing value must be declared NA")
            h, rows = read_csv(mat["file"])
            if not h or h[0] != "feature_id" or len(h) < 2:
                errors.append("Invalid matrix header for " + str(omic))
                continue
            ids = h[1:]
            expected_samples = {id for o, id in assay_keys if o == omic}
            if len(set(ids)) != len(ids):
                errors.append("Duplicate matrix assay_id column for " + str(omic))
            if set(ids) != expected_samples or len(ids) != len(expected_samples):
                errors.append("Matrix/sample assay_id mismatch for " + str(omic))
            measured = set()
            known = {id for o, id in feature_keys if o == omic}
            for r in rows:
                fid = r["feature_id"]
                if not fid or fid in measured:
                    errors.append("Missing/duplicate feature in matrix " + str(omic))
                measured.add(fid)
                for aid in ids:
                    value = r.get(aid, "")
                    if value == "NA":
                        continue
                    if value in ("", None):
                        errors.append("Undeclared empty numeric value for " + str(omic))
                        continue
                    try:
                        num = float(value)
                        if not math.isfinite(num):
                            errors.append("Non-finite numeric value for " + str(omic))
                        if omic == 'transcriptomics' and mat.get('value_type') == 'raw_counts' and (num < 0 or not num.is_integer() or num > 9007199254740991):
                            errors.append('RNA raw_counts must be non-negative safe integers')
                    except ValueError:
                        errors.append("Non-numeric matrix cell for " + str(omic))
            if measured - known:
                errors.append("Matrix contains feature IDs absent from features.csv for " + str(omic))
            if not rows:
                errors.append("Empty matrix for " + str(omic))
    except (OSError, KeyError, TypeError, ValueError, csv.Error) as exc:
        errors.append(str(exc))
    return errors


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 scripts/validate_multiomics_exchange.py EXCHANGE_DIRECTORY")
    try:
        failures = validate(sys.argv[1])
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        failures = [str(exc)]
    if failures:
        for failure in failures:
            print("FAIL: " + failure, file=sys.stderr)
        sys.exit(1)
    print("PASS: MultiOmics exchange structure, identifiers, values and SHA-256")
