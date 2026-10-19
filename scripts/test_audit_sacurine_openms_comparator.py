#!/usr/bin/env python3
"""Package-free negative tests for the Sacurine OpenMS comparator helpers."""
from __future__ import annotations
import importlib.util,tempfile,csv
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("audit",HERE/"audit_sacurine_openms_comparator.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

ref=[("R1",100.0,60.0),("R2",150.0,120.0),("R3",200.0,180.0)]
alt=[("A1",100.001,61.0),("A2",150.001,121.0),("A3",202.0,181.0)]
rows,used=m.reciprocal_matches(ref,alt)
assert [x["match_status"] for x in rows]==["reciprocal_unique","reciprocal_unique","no_coordinate_candidate"]
assert len(used)==2
amb=alt+[("A4",100.0011,60.8)]
rows,_=m.reciprocal_matches(ref,amb)
assert rows[0]["match_status"]=="ambiguous_multiple_candidates"

with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp)/"openms.tsv"
    p.write_text(
      "#MAP\tid\tfilename\tlabel\tsize\n"
      "MAP\t0\t/work/S1.featureXML\t\t2\n"
      "MAP\t1\t/work/S2.featureXML\t\t2\n"
      "#CONSENSUS\trt_cf\tmz_cf\tintensity_cf\tcharge_cf\twidth_cf\tquality_cf\trt_0\tmz_0\tintensity_0\tcharge_0\twidth_0\trt_1\tmz_1\tintensity_1\tcharge_1\twidth_1\n"
      "CONSENSUS\t60\t100\t11\t1\t4\t1\t60\t100\t10\t1\t4\t61\t100.001\t12\t1\t4\n",
      encoding="utf-8")
    maps,cons=m.parse_openms_text(p)
    assert set(maps)=={"0","1"} and len(cons)==1
    assert m.sample_from_map_filename(maps["0"]["filename"])=="S1"

refqc=[{"feature_id":"R1","pooled_qc_rsd":"0.1","blank_over_biological_ratio":"0.25"},
       {"feature_id":"R2","pooled_qc_rsd":"0.5","blank_over_biological_ratio":""}]
altqc=[{"feature_id":"A1","pooled_qc_rsd":0.2,"blank_over_biological_ratio":0.30},
       {"feature_id":"A2","pooled_qc_rsd":0.2,"blank_over_biological_ratio":0.1}]
pairs=[{"reference_feature":"R1","alternative_feature":"A1","match_status":"reciprocal_unique"},
       {"reference_feature":"R2","alternative_feature":"A2","match_status":"reciprocal_unique"}]
a=m.paired_agreement(pairs,refqc,altqc)
assert a["qc_rsd_gt_30pct"]["n_evaluable_both"]==2
assert a["qc_rsd_gt_30pct"]["agreement_fraction"]==0.5
assert a["blank_over_biological_ge_20pct"]["n_evaluable_both"]==1
assert a["blank_over_biological_ge_20pct"]["agreement_fraction"]==1.0
print("PASS: OpenMS comparator parser, reciprocal uniqueness and paired QC negative controls")
