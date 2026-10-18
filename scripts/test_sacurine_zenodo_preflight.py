#!/usr/bin/env python3
"""Offline negative controls; NO fabricated dataset is labelled external evidence."""
import csv
import hashlib
import importlib.util
import io
import pathlib
import tempfile

SCRIPT = pathlib.Path("scripts/qualify_sacurine_toy_zenodo.py")
spec = importlib.util.spec_from_file_location("sacurine_zenodo_preflight", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

ids = [
    ("QC1_014", "QC", "pool", 185),
    ("QC1_008", "QC", "pool", 105),
    ("QC1_002", "QC", "pool", 27),
    ("HU_neg_192", "sample", "sample", 165),
    ("HU_neg_173", "sample", "sample", 148),
    ("HU_neg_157", "sample", "sample", 137),
    ("HU_neg_123", "sample", "sample", 100),
    ("HU_neg_090", "sample", "sample", 75),
    ("HU_neg_048", "sample", "sample", 39),
    ("Blanc16", "blk", "blank", 173),
    ("Blanc10", "blk", "blank", 94),
    ("Blanc05", "blk", "blank", 29),
]
headers = ["sample_name", "class", "sampleType", "injectionOrder", "batch"]
rows = [dict(zip(headers, [name, group, role, str(order), "ne1"]))
        for name, group, role, order in ids]


def encode(items):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=headers, delimiter="\t")
    writer.writeheader()
    writer.writerows(items)
    return output.getvalue().encode()


def expect_error(data, strict=False):
    try:
        mod.validate_source(data, require_source_md5=strict)
    except ValueError:
        return
    raise AssertionError("Unsafe source or tampered metadata unexpectedly accepted.")


manifest = mod.validate_source(encode(rows), require_source_md5=False)
assert manifest["roles"] == {"biological": 6, "pooled_qc": 3, "blank": 3}
assert manifest["source_injection_orders_verified_from_metadata"]
assert manifest["source_injection_order_range"] == [27, 185]
assert manifest["qc_drift_claim"].startswith("NOT_VALIDATED")
assert not manifest["metadata_md5_verified"]
assert not manifest["clinical_certification"]
assert not manifest["independent_external_biological_validation"]
assert len(manifest["all_file_md5_published"]) == 12
assert all(len(x) == 32 for x in manifest["all_file_md5_published"].values())
expect_error(encode(rows), strict=True)

broken = [dict(r) for r in rows]
broken[1]["injectionOrder"] = broken[0]["injectionOrder"]
expect_error(encode(broken))
broken = [dict(r) for r in rows]
broken[1]["class"] = "sample"
expect_error(encode(broken))
broken = [dict(r) for r in rows]
broken[5]["sampleType"] = "pool"
expect_error(encode(broken))
broken = [dict(r) for r in rows]
broken[3]["batch"] = "fabricated_batch"
expect_error(encode(broken))
broken = [dict(r) for r in rows]
broken[4]["sample_name"] = "UNKNOWN"
expect_error(encode(broken))
broken = [dict(r) for r in rows]
broken.pop()
expect_error(encode(broken))
broken = [dict(r) for r in rows]
broken[5]["injectionOrder"] = ""
expect_error(encode(broken))
broken = [dict(r) for r in rows]
broken[5]["injectionOrder"] = "5.5"
expect_error(encode(broken))
expect_error(encode(rows).replace(b"ne1", b"ne2", 1))
expect_error(encode(rows) + b"\x00")

with tempfile.TemporaryDirectory() as temp:
    directory = pathlib.Path(temp)
    try:
        mod.verify_raw_dir(directory)
        raise AssertionError("Incomplete raw instrument archive passed.")
    except ValueError:
        pass
    for name in mod.RAW_MD5:
        (directory / name).write_bytes(b"NOT genuine raw LC-MS data")
    try:
        mod.verify_raw_dir(directory)
        raise AssertionError("Falsified raw chromatograms passed MD5 verification.")
    except ValueError:
        pass

print("PASS: source-controlled Sacurine v3 roles, real-order gate, strong MD5 mismatch and missing-QC negative controls")
