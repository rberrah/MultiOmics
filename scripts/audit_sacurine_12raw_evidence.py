#!/usr/bin/env python3
"""Independent, package-free audit of genuine native-xcms LC-MS study evidence.

Checks original-source manifest, native result counts and *every feature's*
published R-QA flags numerically. Does not see / reconstruct original human
spectra and never certifies instrument calibration or a candidate biomarker.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


def _read_json(path: Path) -> dict:
    x = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(x, dict):
        raise ValueError("Expected structured JSON report: " + path.name)
    return x


def _number(value: str) -> float | None:
    if value == "":
        return None
    n = float(value)
    if not math.isfinite(n) or n < 0:
        raise ValueError("Invalid measured quantitative QA number")
    return n


def _flag(value: str) -> bool:
    if value not in ("TRUE", "FALSE"):
        raise ValueError("Boolean scientific QA flags must be explicit, not missing")
    return value == "TRUE"


def _count(value: str, maximum: int) -> int:
    if not value.isdigit() or int(value) > maximum:
        raise ValueError("Invalid number of authentic finite sample measurements")
    return int(value)


def verify(source_report: Path, native_report: Path,
           qa_report: Path, feature_csv: Path) -> dict:
    source = _read_json(source_report)
    native = _read_json(native_report)
    qa = _read_json(qa_report)
    if (source.get("record") != "https://zenodo.org/records/10130758"
        or source.get("raw_source_verified") is not True
        or source.get("metadata_verified") is not True
        or source.get("raw_count") != 12
        or source.get("source_injection_orders_verified") is not True):
        raise ValueError("Source data were not independently authenticated")
    originals = source.get("files", [])
    if len(originals) != 12 or len({f.get("name") for f in originals}) != 12:
        raise ValueError("Twelve distinct genuine mzML inputs are required")
    expected_roles = {"biological": 6, "pooled_qc": 3, "blank": 3}
    if dict(Counter(f.get("role") for f in originals)) != expected_roles:
        raise ValueError("Original sources do not contain 6/3/3 biological/QC/blank files")
    for f in originals:
        if not f["name"].endswith(".mzML") or f.get("bytes", 0) <= 0 or (
            len(f.get("published_md5", "")) != 32
            or len(f.get("sha256", "")) != 64
        ):
            raise ValueError("Original raw file identity/hash evidence missing")
    if (native.get("status") != "ok" or native.get("samples") != 12
        or native.get("features", 0) < 1
        or native.get("chromatographic_peaks", 0) < 1
        or native.get("blank_samples_excluded_from_feature_definition") != 3):
        raise ValueError("Native R xcms fitting or authentic blank grouping incomplete")
    missing_fraction = native.get("matrix_missing_fraction")
    if not isinstance(missing_fraction, (int, float)) or not 0 <= missing_fraction <= 1:
        raise ValueError("Native xcms measurement missingness invalid")
    if (qa.get("real_biological") != 6
        or qa.get("real_pooled_qc") != 3
        or qa.get("real_blanks") != 3
        or qa.get("raw_spectra_md5_verified") is not True
        or qa.get("source_published_md5_verified") is not True
        or qa.get("numeric_export_replay_equal") is not True
        or qa.get("native_grouped_features") != native["features"]
        or qa.get("native_chromatographic_peaks") != native["chromatographic_peaks"]):
        raise ValueError("Native R scientific evidence disagrees with authenticated source")
    orders = qa.get("original_acquisition_positions")
    if (not isinstance(orders, list) or len(orders) != 12
        or len(set(orders)) != 12 or min(orders) != 27 or max(orders) != 185
        or sorted(orders) == list(range(1, 13))
        or qa.get("all_original_injection_orders_verified") is not True):
        raise ValueError("Original real injection order replaced or lost")
    if source.get("source_head") != qa.get("source_sha"):
        raise ValueError("Native R and source hash verification used different source revisions")
    if any(qa.get(key) is not False for key in (
        "scientific_certified", "clinical_certification",
        "biological_group_test_performed", "feature_identification_validated",
        "omics_pk_subject_linkage_validated"
    )) or qa.get("qc_injection_count_too_low_for_drift_screen") is not True:
        raise ValueError("Evidence overclaims inferential, biological, drift or clinical validation")

    ids: set[str] = set()
    stats = Counter()
    with feature_csv.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {
            "feature_id", "finite_biological_values",
            "finite_pooled_qc_values", "finite_blank_values",
            "biological_missing_fraction", "pooled_qc_rsd",
            "pooled_qc_precision_status", "blank_over_biological_ratio",
            "blank_background_status", "pooled_qc_drift_status",
            "flag_qc_variation", "flag_blank_signal",
            "flag_blank_not_evaluable", "flag_qc_not_evaluable",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Missing complete native feature-level QC audit columns")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Duplicate native QC CSV field names")
        for row in reader:
            name = row["feature_id"]
            if not name or name in ids or None in row:
                raise ValueError("Duplicate, missing or malformed feature provenance")
            ids.add(name)
            bio = _count(row["finite_biological_values"], 6)
            qc = _count(row["finite_pooled_qc_values"], 3)
            blank = _count(row["finite_blank_values"], 3)
            miss = _number(row["biological_missing_fraction"])
            if miss is None or abs(miss - (6 - bio) / 6) > 1e-9:
                raise ValueError("Feature biological missingness disagrees with real six acquisitions")
            rsd = _number(row["pooled_qc_rsd"])
            ratio = _number(row["blank_over_biological_ratio"])
            qc_evaluable = rsd is not None
            blank_evaluable = ratio is not None
            if qc_evaluable and (qc != 3 or row["pooled_qc_precision_status"] != "rsd_evaluable"):
                raise ValueError("QC RSD observed without exactly three measured pooled-QC values")
            if qc < 3 and (qc_evaluable or not _flag(row["flag_qc_not_evaluable"])):
                raise ValueError("Insufficient measured pooled QCs reported as satisfactory")
            if blank_evaluable and (blank < 2 or bio < 2 or
                                    row["blank_background_status"] != "ratio_evaluable"):
                raise ValueError("Blank ratio present without enough original measured intensities")
            if _flag(row["flag_qc_variation"]) != (qc_evaluable and rsd > .30):
                raise ValueError("Per-feature QC variability flag does not match prespecified 30% screen")
            flag_blank = (blank_evaluable and ratio >= .20) or (
                row["blank_background_status"] == "positive_blank_with_zero_biological_median"
            )
            if _flag(row["flag_blank_signal"]) != flag_blank:
                raise ValueError("Per-feature blank background flag does not match 20% screen")
            if _flag(row["flag_blank_not_evaluable"]) != (
                not blank_evaluable and
                row["blank_background_status"] != "positive_blank_with_zero_biological_median"
            ):
                raise ValueError("Absent blank data wrongly interpreted as a zero-background pass")
            if row["pooled_qc_drift_status"] != "too_few_pooled_qc_injections":
                raise ValueError("Three true pooled QCs wrongly declared drift-calibrated")
            stats["features"] += 1
            stats["qc_evaluable"] += int(qc_evaluable)
            stats["qc_rsd_over_30pct"] += int(qc_evaluable and rsd > .30)
            stats["blank_evaluable"] += int(blank_evaluable)
            stats["blank_over_20pct"] += int(blank_evaluable and ratio >= .20)
            stats["blank_not_evaluable"] += int(_flag(row["flag_blank_not_evaluable"]))
            stats["both_evaluable"] += int(qc_evaluable and blank_evaluable)
            stats["both_below_descriptive_screens"] += int(
                qc_evaluable and blank_evaluable and rsd <= .30 and ratio < .20)
            stats["biological_missing_over_20pct"] += int(miss > .20)
    if stats["features"] != native["features"]:
        raise ValueError("Native grouped features are missing from exported QA data")
    crosscheck = {
        "qc_rsd_features_evaluable": "qc_evaluable",
        "qc_rsd_gt30pct": "qc_rsd_over_30pct",
        "blank_ratio_features_evaluable": "blank_evaluable",
        "blank_ratio_ge20pct": "blank_over_20pct",
    }
    if any(qa.get(json_key) != stats[value_key]
           for json_key, value_key in crosscheck.items()):
        raise ValueError("Native R QA summary counters disagree with independent per-feature replay")
    denominator = stats["features"]
    result = {
        "schema": "multiomics.sacurine-v3-12raw-independent-evidence-replay.v1",
        "source": source["record"], "source_sha": qa["source_sha"],
        "original_spectra_authenticated": True,
        "real_biological": 6, "real_pooled_qc": 3, "real_blanks": 3,
        "native_chromatographic_peaks": native["chromatographic_peaks"],
        "native_grouped_features": denominator,
        "native_matrix_missing_fraction": missing_fraction,
        "observed_feature_counts": dict(stats),
        "qc_evaluable_fraction_of_all_features": stats["qc_evaluable"] / denominator,
        "qc_high_rsd_fraction_of_evaluable_features": (
            stats["qc_rsd_over_30pct"] / stats["qc_evaluable"]
            if stats["qc_evaluable"] else None
        ),
        "blank_ratio_evaluable_fraction_of_all_features": (
            stats["blank_evaluable"] / denominator
        ),
        "blank_high_fraction_of_evaluable_features": (
            stats["blank_over_20pct"] / stats["blank_evaluable"]
            if stats["blank_evaluable"] else None
        ),
        "review_status": "scientific_QC_review_required",
        "feature_exclusions_automatically_applied": False,
        "raw_MS_spectral_replay_performed_by_this_Python_script": False,
        "biological_validation": False,
        "clinical_validation": False,
        "omics_to_PK_validation": False,
        "limitations": (
            "Software and instrument preprocessing integrity only; 3 sparse genuine QC "
            "injections do not validate drift correction; non-evaluable blank values are "
            "unknown, not zero contamination; QC/blank screens are descriptive, not "
            "calibrated feature-selection, molecular identification or statistical FDR."
        ),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-report", required=True, type=Path)
    parser.add_argument("--native-report", required=True, type=Path)
    parser.add_argument("--qa-report", required=True, type=Path)
    parser.add_argument("--feature-csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    result = verify(args.source_report, args.native_report,
                    args.qa_report, args.feature_csv)
    if args.out.exists() and (args.out.is_file() or any(args.out.iterdir())):
        raise SystemExit("Do not overwrite previous scientific qualification artifacts")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "sacurine_independent_qc_evidence_review.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    counts = result["observed_feature_counts"]
    lines = [
        "# Independent quantitative review of actual 12-injection Sacurine source",
        "",
        f"Source: {result['source']} — code SHA {result['source_sha']}",
        "",
        f"- Authenticated 6 urine / 3 pooled-QC / 3 actual blanks, "
        f"{result['native_chromatographic_peaks']} peaks / {result['native_grouped_features']} features",
        f"- QC RSD available: {counts['qc_evaluable']}; RSD >30%: {counts['qc_rsd_over_30pct']}",
        f"- Blank/urine ratio available: {counts['blank_evaluable']}; >=20%: {counts['blank_over_20pct']}",
        f"- Blank not evaluable: {counts['blank_not_evaluable']}; **not** blank-negative",
        f"- Both screens evaluable: {counts['both_evaluable']}; "
        f"below both descriptive screens: {counts['both_below_descriptive_screens']}",
        "",
        "**Research decision: SCIENTIFIC QC REVIEW REQUIRED.**",
        "",
        "The feature screens were prespecified but are not an optimised molecular",
        "identification filter. No feature was removed by this audit, and no downstream",
        "hypothesis test or PK link was run. A green compute workflow is NOT proof",
        "that the features satisfy an experimental QC acceptance plan.",
        "",
    ]
    (args.out / "SACURINE_INDEPENDENT_QC_REVIEW.md").write_text(
        "\n".join(lines), encoding="utf-8")
    print("PASS: independently reconciled all native xcms QC flags; scientific review REQUIRED")


if __name__ == "__main__":
    main()
