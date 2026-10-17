#!/usr/bin/env python3
"""Genuine GSE22871 native Cyp3a11 probe intensities and within-array ranks.

This is NOT Nutrimouse model projection, cross-platform normalization, a
validated gene signature or PK evidence. Selection uses source identifiers
alone, never the genotype/dose or expression response. One row per mouse.
"""
import csv
import gzip
import hashlib
import json
import math
import re
import sys
from pathlib import Path

from qualify_geo_gse22871 import PINNED_GEO_SOURCE_SHA256, SAMPLE_IDS, require
from qualify_gpl1261_annotation import PINNED_GPL1261_SHA256

N_PROBES = 45101
ANCHOR = "cyp3a11"


def source_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def symbols(raw):
    return {v.casefold() for v in re.split(r"[|;,/\s]+", raw or "") if v}


def select_native_candidates(rows):
    """Read ALL external IDs first; prohibit phenotype-driven probe picking."""
    candidates = [r for r in rows
                  if ANCHOR in symbols(r["geo_symbol_raw"]) or
                  ANCHOR in symbols(r["bioc_symbols_existing_audit"])]
    require(0 < len(candidates) <= 10, "Missing or unexpectedly many predeclared Cyp3a11 candidates")
    require(len({r["probe_id"] for r in candidates}) == len(candidates),
            "Duplicated external candidate probe identities")
    return sorted(candidates, key=lambda r: r["probe_id"])


def expression_rows(source):
    """Stream original processed GEO table without rewriting/renormalizing data."""
    with gzip.open(source, "rt", encoding="utf-8", newline="") as f:
        table = False
        header = None
        ended = False
        for line in f:
            if not table:
                if line.strip() == "!series_matrix_table_begin":
                    table = True
                continue
            if header is None:
                header = next(csv.reader([line], delimiter="\t"))
                require(header[0] == "ID_REF" and header[1:] == list(SAMPLE_IDS),
                        "Genuine matrix GSM header or order changed")
                continue
            if line.strip() == "!series_matrix_table_end":
                ended = True
                break
            row = next(csv.reader([line], delimiter="\t"))
            require(len(row) == 31, "Truncated original probe expression row")
            try:
                values = [float(x) for x in row[1:]]
            except ValueError as e:
                raise ValueError("Original source expression non-numeric") from e
            require(all(math.isfinite(v) for v in values), "Nonfinite probe intensities")
            yield row[0], values
        require(ended, "Original GEO matrix terminator missing")


def native_evidence(source, candidates):
    target = {x["probe_id"] for x in candidates}
    measured = {}
    seen = set()
    total = 0
    for name, values in expression_rows(source):
        require(name not in seen, "Duplicated original GEO probe ID")
        seen.add(name)
        total += 1
        if name in target:
            measured[name] = values
    require(total == N_PROBES and set(measured) == target,
            "Original probe total or predeclared candidates changed")
    # A raw percentile is a DESCRIPTIVE within-array rank among all 45,101
    # hybridization probes. It is not a normalized across-platform gene score.
    fewer = {name: [0] * 30 for name in target}
    ties = {name: [0] * 30 for name in target}
    recounted = 0
    for _, values in expression_rows(source):
        recounted += 1
        for name in target:
            for i, value in enumerate(values):
                if value < measured[name][i]:
                    fewer[name][i] += 1
                elif value == measured[name][i]:
                    ties[name][i] += 1
    require(recounted == total, "Independent second pass did not read the full source")
    ranks = {}
    for name in target:
        require(all(n >= 1 for n in ties[name]), "Missing self-match in within-array rank")
        ranks[name] = [(fewer[name][i] + (ties[name][i] - 1)/2) / (total - 1)
                       for i in range(30)]
        require(all(0 <= x <= 1 for x in ranks[name]), "Invalid rank percentile")
    return measured, ranks


def run(source_dir, bioc_dir):
    source_dir, bioc_dir = Path(source_dir), Path(bioc_dir)
    source = source_dir / "GSE22871_series_matrix.txt.gz"
    require(source_sha(source) == PINNED_GEO_SOURCE_SHA256, "Original GEO input SHA mismatch")
    intake = json.loads((source_dir / "external_geo_intake_audit.json").read_text())
    agree = json.loads((bioc_dir / "original_NCBI_vs_Bioc_Entrez_ID" /
                        "GPL1261_entrez_id_concordance_audit.json").read_text())
    require(intake["n_biological_subjects"] == 30 and
            intake["n_processed_affymetrix_probe_ids"] == N_PROBES,
            "Unverified external biological cohort")
    require(agree["source_original_geo_series_sha256"] == PINNED_GEO_SOURCE_SHA256 and
            agree["source_original_gpl1261_sha256"] == PINNED_GPL1261_SHA256 and
            agree["source_probe_measurements"] == N_PROBES,
            "External annotation evidence is not pinned to original source")
    path = bioc_dir / "original_NCBI_vs_Bioc_Entrez_ID" / "GPL1261_all_probe_entrez_concordance.csv"
    with path.open(encoding="utf-8", newline="") as f:
        reference = list(csv.DictReader(f))
    require(len(reference) == N_PROBES, "Incomplete original independent stable-ID crosswalk")
    candidates = select_native_candidates(reference)
    measured, rank = native_evidence(source, candidates)
    out = bioc_dir / "native_external_probe_observations"
    out.mkdir(parents=True, exist_ok=True)
    subjects = intake["samples"]
    require([s["assay_id"] for s in subjects] == list(SAMPLE_IDS),
            "Native numerical assays not paired to original real sample metadata")
    results = []
    for candidate in candidates:
        probe = candidate["probe_id"]
        for i, sub in enumerate(subjects):
            results.append({
                "subject_id": sub["subject_id"], "assay_id": sub["assay_id"],
                "probe_id": probe, "candidate_gene": "Cyp3a11",
                "native_processed_log2_RMA_intensity": format(measured[probe][i], ".17g"),
                "within_array_probe_intensity_percentile": format(rank[probe][i], ".17g"),
                "genotype": sub["genotype"], "pfos_dose_mg_kg_day": sub["dose_mg_per_kg_day"],
                "entrez_concordance": candidate["entrez_concordance"],
                "symbol_concordance": candidate["symbol_concordance"],
                "cross_platform_model_score": "NOT_COMPUTED",
            })
    fname = out / "native_Cyp3a11_probe_measurements_not_model_scores.csv"
    with fname.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)
    high_confidence = len(candidates) == 1 and all(
        c["entrez_concordance"] == "exact_entrez_set" and
        c["symbol_concordance"] == "exact_symbol_set_agreement" and
        c["multi_entrez_ambiguous"] == "False" for c in candidates
    )
    report = {
        "study": "GSE22871",
        "original_geo_sha256": PINNED_GEO_SOURCE_SHA256,
        "original_annotation_sha256": PINNED_GPL1261_SHA256,
        "numeric_table_sha256": source_sha(fname),
        "n_independent_mice": len(subjects),
        "n_original_measured_probes": N_PROBES,
        "n_anchor_candidates": len(candidates),
        "n_native_probe_measurement_rows": len(results),
        "candidate_probe_ids": [c["probe_id"] for c in candidates],
        "two_reference_unique_anchor_mapping": high_confidence,
        "data_scale": "source_processed_RMA_log2_affymetrix",
        "rank_definition": "midrank among all 45101 processed probes within each individual array",
        "selection_used_dose_or_genotype": False,
        "train_or_test_cross_platform_scaling_fitted": False,
        "nutrimouse_frozen_model_applied": False,
        "external_biological_signature_scores_exported": 0,
        "real_PK_concentration_outcomes": False,
        "status": "EXTERNAL_NATIVE_MEASUREMENT_ONLY_PROJECTION_BLOCKED",
        "limits": [
            "Cross-platform macroarray versus Affymetrix expression levels are not interchangeable",
            "Probe signal and probe percentile are descriptive array measurements, not CYP activity",
            "The genotype/PFOS factors are retained for future prespecified exploratory checks only",
            "No inference of independent biological or clinical PK validation",
        ],
    }
    (out / "external_native_measurement_audit.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("PASS: genuine native probe intensities on 30 unique GEO subjects;",
          len(candidates), "source-driven Cyp3a11 candidates; two-pass rank replay;",
          "cross-platform model score BLOCKED")


def selftest():
    c = [
        {"probe_id": "a", "geo_symbol_raw": "Cyp3a11", "bioc_symbols_existing_audit": "",
         "entrez_concordance": "exact_entrez_set", "symbol_concordance": "exact_symbol_set_agreement",
         "multi_entrez_ambiguous": "False"},
        {"probe_id": "b", "geo_symbol_raw": "Cyp3a1", "bioc_symbols_existing_audit": "Other",
         "entrez_concordance": "exact_entrez_set", "symbol_concordance": "exact_symbol_set_agreement",
         "multi_entrez_ambiguous": "False"},
    ]
    assert [x["probe_id"] for x in select_native_candidates(c)] == ["a"]
    assert "cyp3a11" not in symbols("Cyp3a110")
    try:
        select_native_candidates([c[0], c[0]])
        raise AssertionError("Duplicate candidate silently accepted")
    except ValueError:
        pass
    try:
        select_native_candidates([c[1]])
        raise AssertionError("Missing prespecified anchor silently accepted")
    except ValueError:
        pass
    x = [4., 7., 9., 10., 15.]
    rank = lambda a, v: (sum(z < v for z in a) +
                          (sum(z == v for z in a) - 1)/2)/(len(a)-1)
    assert rank(x, 9) == rank([v*3+100 for v in x], 127) == .5
    assert rank([1,2,2,3,4], 2) == .375
    print("PASS: candidate identity, absent/duplicate negative controls, and monotone-rank self-tests")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        selftest()
    elif len(sys.argv) == 3:
        run(sys.argv[1], sys.argv[2])
    else:
        sys.exit("Usage: qualify_geo_native_probe_evidence.py --self-test | GEO_DIR BIOC_AUDIT_DIR")
