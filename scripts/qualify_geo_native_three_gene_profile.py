#!/usr/bin/env python3
"""Independent GSE22871 *target-native* three-gene measurements.

A descriptive midrank index, NOT the four-gene Nutrimouse frozen model,
NOT calibrated macroarray-to-Affymetrix projection, and NOT a PK covariate.
Gene/probe selection uses only the independently checked identity crosswalk.
"""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

from qualify_geo_gse22871 import PINNED_GEO_SOURCE_SHA256, SAMPLE_IDS, require
from qualify_geo_native_probe_evidence import expression_rows, N_PROBES

SOURCE_GENES = ("CYP3A11", "FAS", "CYP4A14")
EXCLUDED_GENE = "CYP4A10"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def locked_probes(rows):
    """Never recover a missing gene via partial matches or dose association."""
    require(len(rows) == 4, "Frozen source panel must have exactly four genes")
    lookup = {r["source_gene_symbol"].upper(): r for r in rows}
    require(len(lookup) == 4 and set(lookup) == set(SOURCE_GENES) | {EXCLUDED_GENE},
            "Frozen source genes changed")
    require(rows[0]["source_gene_symbol"].upper() == "CYP3A11",
            "Original training-selected anchor changed")
    require(lookup["FAS"]["resolved_target_gene_symbol"] == "Fasn" and
            lookup["FAS"]["source_entrez_ids"] == "14104",
            "FAS/Fasn original RefSeq identity not proved")
    require(lookup[EXCLUDED_GENE]["strict_probe_count"] == "0" and
            not lookup[EXCLUDED_GENE]["matching_geo_probe_ids"],
            "Unexpected CYP4A10 mapping: requalify original frozen model before proceeding")
    output = {}
    for gene in SOURCE_GENES:
        r = lookup[gene]
        probe = r["matching_geo_probe_ids"]
        require(r["strict_probe_count"] == "1" and
                r["technical_crosswalk_status"] == "gene_to_target_probe_strict_unique" and
                probe and "|" not in probe,
                "No unique independently annotated target probe for " + gene)
        output[gene] = probe
    require(len(set(output.values())) == 3, "Multiple genes mapped to one target probe")
    return output


def profile(original, targets):
    """Two full passes: original 45,101 rows, per-mouse midranks, no fitted model."""
    seen = set()
    observed = {}
    n = 0
    for probe, values in expression_rows(original):
        require(probe not in seen, "Duplicate original probe")
        seen.add(probe)
        n += 1
        if probe in targets:
            observed[probe] = values
    require(n == N_PROBES and set(observed) == set(targets),
            "Original complete GEO source or three frozen probes missing")
    less = {p: [0] * 30 for p in targets}
    equal = {p: [0] * 30 for p in targets}
    count = 0
    for _, vals in expression_rows(original):
        count += 1
        for probe in targets:
            for i, value in enumerate(vals):
                if value < observed[probe][i]:
                    less[probe][i] += 1
                elif value == observed[probe][i]:
                    equal[probe][i] += 1
    require(count == N_PROBES, "Second independent full-matrix pass incomplete")
    ranks = {}
    for probe in targets:
        require(all(n >= 1 for n in equal[probe]), "Target probe absent in rank calculation")
        ranks[probe] = [(less[probe][i] + (equal[probe][i] - 1) / 2) / (N_PROBES - 1)
                        for i in range(30)]
        require(all(0 <= x <= 1 and math.isfinite(x) for x in ranks[probe]),
                "Invalid within-array probe percentile")
    return observed, ranks


def write_csv(path, records):
    require(records, "Missing individual assay records")
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def run(root):
    root = Path(root)
    original = root / "GSE22871_series_matrix.txt.gz"
    meta = root / "external_geo_intake_audit.json"
    bridge = (root / "independent-bioc-concordance" /
              "real_source_to_target_gene_bridge")
    crosswalk = bridge / "locked_panel_to_external_probe_crosswalk.csv"
    readiness = json.loads((bridge / "frozen_panel_external_transport_readiness.json").read_text())
    intake = json.loads(meta.read_text())
    require(sha256(original) == PINNED_GEO_SOURCE_SHA256 and
            readiness["external_original_gse22871_sha256"] == PINNED_GEO_SOURCE_SHA256 and
            readiness["candidate_table_sha256"] == sha256(crosswalk),
            "Source or immutable probe-identity crosswalk SHA mismatch")
    require(readiness["strictly_mapped_feature_count"] == 3 and
            readiness["panel_n_frozen_genes"] == 4 and
            not readiness["cross_platform_reference_bridge_fitted"] and
            not readiness["direct_source_zscore_to_target_is_justified"] and
            readiness["n_external_signature_scores_exported"] == 0 and
            readiness["n_PK_concentration_outcomes"] == 0,
            "Cross-platform claim or frozen panel mapping status changed")
    subjects = intake["samples"]
    require(intake["n_biological_subjects"] == 30 and
            [s["assay_id"] for s in subjects] == list(SAMPLE_IDS) and
            len({s["subject_id"] for s in subjects}) == 30,
            "Original mouse assay identities or expected sample count changed")
    with crosswalk.open(newline="", encoding="utf-8") as f:
        targets = locked_probes(list(csv.DictReader(f)))
    measured, ranks = profile(original, set(targets.values()))
    probes = []
    index = []
    for i, sub in enumerate(subjects):
        values = []
        for gene in SOURCE_GENES:
            probe = targets[gene]
            rank = ranks[probe][i]
            values.append(rank)
            probes.append({
                "subject_id": sub["subject_id"], "assay_id": sub["assay_id"],
                "original_frozen_feature_label": gene,
                "target_platform_gene": "Fasn" if gene == "FAS" else gene,
                "target_probe_id": probe,
                "original_processed_log2_RMA": format(measured[probe][i], ".17g"),
                "within_array_midrank_percentile": format(rank, ".17g"),
            })
        index.append({
            "subject_id": sub["subject_id"], "assay_id": sub["assay_id"],
            "target_native_three_probe_mean_midrank": format(sum(values) / 3, ".17g"),
            "scale": "within_array_probe_midrank_0_to_1",
            "classification": "descriptive_not_transferred_model_score",
        })
    out = root / "external_native_three_gene_descriptive_profile"
    out.mkdir(parents=True, exist_ok=True)
    probe_csv = out / "individual_native_three_gene_probe_measurements.csv"
    index_csv = out / "individual_native_three_gene_descriptive_index.csv"
    write_csv(probe_csv, probes)
    write_csv(index_csv, index)
    audit = {
        "study": "GSE22871",
        "original_gse22871_sha256": PINNED_GEO_SOURCE_SHA256,
        "crosswalk_sha256": sha256(crosswalk),
        "source_frozen_model_sha256": readiness["source_train_only_model_sha256"],
        "original_technical_panel_genes": 4,
        "unambiguously_mapped_native_probe_genes": 3,
        "excluded_source_feature": EXCLUDED_GENE,
        "target_probe_ids_by_original_gene": targets,
        "n_independent_subjects": 30,
        "n_native_probe_rows": len(probes),
        "n_descriptive_per_assay_index_rows": len(index),
        "raw_native_probe_csv_sha256": sha256(probe_csv),
        "descriptive_index_csv_sha256": sha256(index_csv),
        "selection_uses_external_dose_genotype_outcomes": False,
        "native_index_preprocessing_fits_across_subjects": False,
        "nutrimouse_frozen_four_gene_model_applied": False,
        "frozen_model_cross_platform_projections": 0,
        "pk_covariates_exported": 0,
        "status": "TARGET_NATIVE_DESCRIPTIVE_INDEX_ONLY_NOT_EXTERNALLY_VALIDATED_SIGNATURE",
        "limits": [
            "Exactly three mapped genes; CYP4A10 excluded, not imputed",
            "Per-array midrank over 45101 probes is platform-relative expression, not enzyme activity",
            "This arithmetic three-gene index is a DIFFERENT estimator from the original frozen four-gene Nutrimouse z score",
            "No source-to-target bridge, paired technology calibration, external biological replication or PK concentrations",
            "Do not write these exploratory rows as individual_signatures.csv or label fit_scope as a frozen model projection",
        ],
    }
    (out / "native_three_gene_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("PASS: 30 genuine mice, 3 exact mapped probes, 90 native measurements; "
          "descriptive index independently replayable; 0 transported Nutrimouse/PK scores")


def selftest():
    mk = lambda gene, probe, count: {
        "source_gene_symbol": gene,
        "resolved_target_gene_symbol": "Fasn" if gene == "FAS" else gene,
        "source_entrez_ids": "14104" if gene == "FAS" else "1",
        "strict_probe_count": str(count),
        "matching_geo_probe_ids": probe,
        "technical_crosswalk_status": ("gene_to_target_probe_strict_unique"
                                       if count == 1 else
                                       "no_unambiguous_two_reference_target_probe"),
    }
    base = [mk("CYP3A11", "p1", 1), mk("FAS", "p2", 1),
            mk("CYP4A14", "p3", 1), mk("CYP4A10", "", 0)]
    assert locked_probes(base) == {"CYP3A11": "p1", "FAS": "p2", "CYP4A14": "p3"}
    for changes in (
        {1: {"source_entrez_ids": "14102"}},
        {2: {"matching_geo_probe_ids": "p1"}},
        {3: {"matching_geo_probe_ids": "p4", "strict_probe_count": "1"}},
        {0: {"strict_probe_count": "0"}},
    ):
        damaged = [dict(x) for x in base]
        for i, patch in changes.items():
            damaged[i].update(patch)
        try:
            locked_probes(damaged)
        except ValueError:
            pass
        else:
            raise AssertionError("An ambiguous/falsified gene-to-probe mapping was accepted")
    print("PASS: FAS alias, duplicate gene probe, CYP4A10 rescue and absent anchor rejected")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        selftest()
    elif len(sys.argv) == 2:
        run(sys.argv[1])
    else:
        sys.exit("Usage: qualify_geo_native_three_gene_profile.py --self-test | GEO_OUTPUT_DIR")
