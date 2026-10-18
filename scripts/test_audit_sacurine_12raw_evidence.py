#!/usr/bin/env python3
"""Deterministic negative controls for independent Sacurine instrument QA replay."""
from __future__ import annotations
import csv
import json
import tempfile
from pathlib import Path
from audit_sacurine_12raw_evidence import verify

names = (
    [(f"HU_neg_{i:03d}.mzML", "biological") for i in (48, 90, 123, 157, 173, 192)]
    + [(f"QC1_{i:03d}.mzML", "pooled_qc") for i in (2, 8, 14)]
    + [(f"Blanc{i:02d}.mzML", "blank") for i in (5, 10, 16)]
)
source = {
    "record": "https://zenodo.org/records/10130758",
    "raw_source_verified": True,
    "metadata_verified": True,
    "raw_count": 12,
    "source_injection_orders_verified": True,
    "source_head": "fixture_source_only",
    "files": [
        {"name": n, "role": role, "bytes": 100, "published_md5": "a"*32,
         "sha256": "b"*64}
        for n, role in names
    ]
}
native = {"status": "ok", "samples": 12, "features": 2,
          "chromatographic_peaks": 20,
          "blank_samples_excluded_from_feature_definition": 3,
          "matrix_missing_fraction": 0.1}
qa = {
    "real_biological": 6, "real_pooled_qc": 3, "real_blanks": 3,
    "raw_spectra_md5_verified": True, "source_published_md5_verified": True,
    "numeric_export_replay_equal": True,
    "native_grouped_features": 2, "native_chromatographic_peaks": 20,
    "all_original_injection_orders_verified": True,
    "original_acquisition_positions": [185,105,27,165,148,137,100,75,39,173,94,29],
    "source_sha": "fixture_source_only",
    "scientific_certified": False, "clinical_certification": False,
    "biological_group_test_performed": False,
    "feature_identification_validated": False,
    "omics_pk_subject_linkage_validated": False,
    "qc_injection_count_too_low_for_drift_screen": True,
    "qc_rsd_features_evaluable": 2, "qc_rsd_gt30pct": 1,
    "blank_ratio_features_evaluable": 1, "blank_ratio_ge20pct": 0,
}
rows = [
    {"feature_id":"FT001","finite_biological_values":"6",
     "finite_pooled_qc_values":"3","finite_blank_values":"3",
     "biological_missing_fraction":"0","pooled_qc_rsd":"0.2",
     "pooled_qc_precision_status":"rsd_evaluable",
     "blank_over_biological_ratio":"0.1",
     "blank_background_status":"ratio_evaluable",
     "pooled_qc_drift_status":"too_few_pooled_qc_injections",
     "flag_qc_variation":"FALSE","flag_blank_signal":"FALSE",
     "flag_blank_not_evaluable":"FALSE","flag_qc_not_evaluable":"FALSE"},
    {"feature_id":"FT002","finite_biological_values":"6",
     "finite_pooled_qc_values":"3","finite_blank_values":"0",
     "biological_missing_fraction":"0","pooled_qc_rsd":"0.4",
     "pooled_qc_precision_status":"rsd_evaluable",
     "blank_over_biological_ratio":"",
     "blank_background_status":"too_few_finite_blank_values",
     "pooled_qc_drift_status":"too_few_pooled_qc_injections",
     "flag_qc_variation":"TRUE","flag_blank_signal":"FALSE",
     "flag_blank_not_evaluable":"TRUE","flag_qc_not_evaluable":"FALSE"},
]


def must_fail(func):
    try:
        func()
    except (ValueError, KeyError, TypeError):
        return
    raise AssertionError("Corrupted scientific QA evidence was accepted")


with tempfile.TemporaryDirectory() as dirname:
    folder = Path(dirname)
    files = [folder / name for name in ("source.json","native.json","qa.json","features.csv")]
    def save():
        for path,data in zip(files[:3],(source,native,qa)):
            path.write_text(json.dumps(data)+"\n",encoding="utf-8")
        with files[3].open("w",newline="",encoding="utf-8") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    def check():
        save()
        return verify(*files)

    good = check()
    assert good["observed_feature_counts"]["features"] == 2
    assert good["observed_feature_counts"]["qc_rsd_over_30pct"] == 1
    assert good["observed_feature_counts"]["blank_not_evaluable"] == 1
    assert good["observed_feature_counts"]["both_evaluable"] == 1
    assert good["review_status"] == "scientific_QC_review_required"
    assert not good["feature_exclusions_automatically_applied"]

    rows[1]["flag_blank_not_evaluable"]="FALSE"
    must_fail(check)
    rows[1]["flag_blank_not_evaluable"]="TRUE"
    rows[1]["flag_qc_variation"]="FALSE"
    must_fail(check)
    rows[1]["flag_qc_variation"]="TRUE"
    rows[0]["biological_missing_fraction"]="0.2"
    must_fail(check)
    rows[0]["biological_missing_fraction"]="0"
    qa["qc_rsd_gt30pct"]=0
    must_fail(check)
    qa["qc_rsd_gt30pct"]=1
    qa["scientific_certified"]=True
    must_fail(check)
    qa["scientific_certified"]=False
    qa["original_acquisition_positions"]=list(range(1,13))
    must_fail(check)
    qa["original_acquisition_positions"]=[185,105,27,165,148,137,100,75,39,173,94,29]
    source["files"][0]["role"]="blank"
    must_fail(check)
    source["files"][0]["role"]="biological"
    assert check()["observed_feature_counts"]["features"] == 2

print("PASS: independent Python numerical LC-MS QA counters, unobserved blanks and fabricated validation negative controls")
