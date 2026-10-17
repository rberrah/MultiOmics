#!/usr/bin/env python3
"""Exploratory predeclared native Cyp3a11 high/zero PFOS contrasts in 30 real GEO mice.

Uses one uniquely annotated Affymetrix probe, source log2-RMA intensities,
a fixed 2x3 design, exact 10-choose-5 randomisation tests and two BH q values.
No transfer of Nutrimouse model, no PK endpoints or biological signature.
"""
import csv
import hashlib
import itertools
import json
import math
import statistics
import sys
from pathlib import Path

from qualify_geo_gse22871 import require, PINNED_GEO_SOURCE_SHA256
from qualify_gpl1261_annotation import PINNED_GPL1261_SHA256

PRIMARY = (("wild-type", 10, 0), ("PPARalpha-null", 10, 0))
EXPECTED_ROWS = 30
N_PER_GROUP = 5


def exact_p(high, zero):
    """Exact two-sided permutation using all 252 relabellings of 10 mice."""
    combined = list(high) + list(zero)
    require(len(high) == len(zero) == 5 and len(set(combined)) >= 2,
            "Five real animals required per exposure group and variable probe expression")
    observed = abs(statistics.mean(high) - statistics.mean(zero))
    diffs = []
    for chosen in itertools.combinations(range(10), 5):
        chosen = set(chosen)
        a = [v for i, v in enumerate(combined) if i in chosen]
        b = [v for i, v in enumerate(combined) if i not in chosen]
        diffs.append(abs(statistics.mean(a) - statistics.mean(b)))
    require(len(diffs) == 252, "Incomplete exact permutation universe")
    eps = 1e-12 * max(1, observed)
    return sum(x >= observed - eps for x in diffs) / len(diffs)


def bh(pvals):
    m = len(pvals)
    out = [0.0] * m
    last = 1.0
    for i in sorted(range(m), key=lambda k: pvals[k], reverse=True):
        rank = 1 + sum(pvals[k] < pvals[i] for k in range(m))
        # Equal p-values use largest rank, as in standard BH.
        rank = max(j + 1 for j in range(m) if pvals[j] <= pvals[i])
        last = min(last, pvals[i] * m / rank)
        out[i] = min(1.0, last)
    return out


def compute(rows):
    require(len(rows) == EXPECTED_ROWS, "Not exactly 30 original source animals")
    require(len({r["assay_id"] for r in rows}) == EXPECTED_ROWS,
            "Duplicate biological mouse/assay in external study")
    require(len({r["probe_id"] for r in rows}) == 1,
            "No unambiguous single measured Cyp3a11 probe for contrast")
    require(all(r["cross_platform_model_score"] == "NOT_COMPUTED" for r in rows),
            "A falsely represented external biological score was inserted")
    cells = {}
    for genotype in ("wild-type", "PPARalpha-null"):
        for dose in (0, 3, 10):
            items = [r for r in rows if r["genotype"] == genotype and
                     int(r["pfos_dose_mg_kg_day"]) == dose]
            require(len(items) == N_PER_GROUP and
                    len({r["subject_id"] for r in items}) == N_PER_GROUP,
                    f"Incorrect true replicate count for {genotype} dose {dose}")
            values = [float(r["native_processed_log2_RMA_intensity"]) for r in items]
            require(all(math.isfinite(v) for v in values),
                    "Nonfinite original processed intensity")
            cells[(genotype, dose)] = values
    results = []
    for genotype, high, zero in PRIMARY:
        higher, baseline = cells[(genotype, high)], cells[(genotype, zero)]
        results.append({
            "contrast": f"{genotype}:PFOS_10_vs_0_mg_kg_day",
            "genotype": genotype,
            "n_high": len(higher),
            "n_zero": len(baseline),
            "high_mean_original_log2_RMA": statistics.mean(higher),
            "zero_mean_original_log2_RMA": statistics.mean(baseline),
            "high_minus_zero_native_log2_RMA": statistics.mean(higher) - statistics.mean(baseline),
            "exact_two_sided_permutation_p": exact_p(higher, baseline),
            "permutations": 252,
            "clinical_PK_endpoint": False,
        })
    qvals = bh([r["exact_two_sided_permutation_p"] for r in results])
    for r, q in zip(results, qvals):
        r["BH_q_two_predefined_contrasts"] = q
    return results


def selftest():
    x = [1,2,3,4,5]
    assert exact_p(x, x) == 1
    y = [9,10,11,12,13]
    assert abs(exact_p(x,y) - 2/252) < 1e-14
    assert bh([.01,.03]) == [.02,.03]
    try:
        exact_p([1,2,3,4], y)
    except ValueError:
        pass
    else:
        raise AssertionError("Malformed true biological replicate counts silently accepted")
    print("PASS: 252 exact-label assignments, BH, null case and missing-mouse negative control")


def run(root):
    root = Path(root)
    reference = root / "independent-bioc-concordance" / "native_external_probe_observations"
    source = root / "GSE22871_series_matrix.txt.gz"
    report = json.loads((reference / "external_native_measurement_audit.json").read_text())
    require(hashlib.sha256(source.read_bytes()).hexdigest() == PINNED_GEO_SOURCE_SHA256,
            "Genuine source gzip changed")
    require(report["original_geo_sha256"] == PINNED_GEO_SOURCE_SHA256 and
            report["original_annotation_sha256"] == PINNED_GPL1261_SHA256,
            "Native observations not linked to frozen dual-reference audit")
    require(report["two_reference_unique_anchor_mapping"] and
            report["n_anchor_candidates"] == 1 and
            report["n_independent_mice"] == 30 and
            report["external_biological_signature_scores_exported"] == 0,
            "Only a unique, verified unambiguous real source probe is admissible")
    path = reference / "native_Cyp3a11_probe_measurements_not_model_scores.csv"
    require(hashlib.sha256(path.read_bytes()).hexdigest() == report["numeric_table_sha256"],
            "Frozen upstream native probe evidence modified")
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    results = compute(rows)
    outdir = root / "external_native_PFOS_exploratory_contrasts"
    outdir.mkdir(exist_ok=True)
    csv_path = outdir / "Cyp3a11_native_primary_contrasts.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    report = {
        "source_accession": "GSE22871",
        "original_source_sha256": PINNED_GEO_SOURCE_SHA256,
        "original_GPL1261_annotation_sha256": PINNED_GPL1261_SHA256,
        "native_measurement_input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "n_mice": 30,
        "assay": "source_original_processed_RMA_log2_Affymetrix_probe",
        "design": "predeclared_wildtype_and_PPARalpha_null_10_vs_0mgkg",
        "n_primary_contrasts": 2,
        "test": "exact_2_sided_10_choose_5_permutation_within_genotype",
        "multiple_testing": "BH_on_two_primary_contrasts",
        "results": results,
        "status": "EXPLORATORY_NATIVE_PROBE_EXPOSURE_ASSOCIATION_NOT_CROSS_PLATFORM_VALIDATION",
        "requires_cross_platform_bridge_for_Nutrimouse_projection": True,
        "frozen_Nutrimouse_model_scores_generated": 0,
        "biologically_validated_PK_scores": 0,
        "PK_concentrations_available": False,
        "claim_boundaries": [
            "Exposure-related raw probe intensity is not externally validated CYP enzyme activity",
            "The original paper reports the same cohort; this is an analysis of the public observations, not independent replication of the paper",
            "These two endpoints were defined from the known design and remain exploratory",
            "No test result, whether positive or negative, enables Nutrimouse model transfer",
        ],
    }
    (outdir / "native_contrast_report.json").write_text(json.dumps(report, indent=2)+"\n")
    for rec in results:
        print("GSE22871 genuine native expression", rec["contrast"],
              "difference:", round(rec["high_minus_zero_native_log2_RMA"],6),
              "exact_p:", round(rec["exact_two_sided_permutation_p"],6),
              "BH_q:", round(rec["BH_q_two_predefined_contrasts"],6))
    print("PASS: tested two honest individual-mouse exploratory external expression contrasts; no PK score")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        selftest()
    elif len(sys.argv) == 2:
        run(sys.argv[1])
    else:
        sys.exit("Usage: qualify_geo_native_pfos_contrasts.py --self-test | GSE22871_ORIGINAL_DIR")
