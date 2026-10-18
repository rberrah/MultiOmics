#!/usr/bin/env python3
"""Source authenticity preflight for Zenodo Sacurine toy v3 (DOI 10.5281/zenodo.10130758).

This script checks actual source-provided metadata/MD5; it does NOT perform
xcms peak extraction, quantify carryover/drift, or certify a clinical cohort.
Do not synthesize injection chronology from file lists.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import pathlib
import re
import sys

ZENODO_RECORD = "https://zenodo.org/records/10130758"
DOI = "10.5281/zenodo.10130758"
METADATA_FILE = "sampleMetadata_12samp_completed.tsv"
METADATA_MD5 = "f88044818627dee797c75d5a310f5c8e"
# Published on the immutable Zenodo v3 record, not guessed from file names.
RAW_MD5 = {
    "Blanc05.mzML": "2051f46e74d449990fc24d269ec7be9f",
    "Blanc10.mzML": "bbe6c60ae6e2e3a471f729dcd291a792",
    "Blanc16.mzML": "58f4285e99a000f44005805ba7b727ea",
    "HU_neg_048.mzML": "563edf9f1005ccb4f594011643b64cf3",
    "HU_neg_090.mzML": "0a59b1a64a99e17577e015fcf9693463",
    "HU_neg_123.mzML": "fbbaba73ef7160ecbea72f7061130518",
    "HU_neg_157.mzML": "8da4a03622ad373076c13e1ed2df9dcb",
    "HU_neg_173.mzML": "15714f24e2ff191916b0c57ffbede6c7",
    "HU_neg_192.mzML": "4f4419b678bcf9473c85317931b2cb07",
    "QC1_002.mzML": "86f69f207c7233136c5964fe65b02f3a",
    "QC1_008.mzML": "1f854194a405818091247b90dbf8523e",
    "QC1_014.mzML": "927b4f3e8045d32d62e6ed7d4a6978ee",
}
ROLES = {"QC": ("pool", "pooled_qc"), "blk": ("blank", "blank"),
         "sample": ("sample", "biological")}
REQUIRED = ["sample_name", "class", "sampleType", "injectionOrder", "batch"]
STRICT_META_SHA = True


def md5_bytes(raw: bytes) -> str:
    return hashlib.md5(raw).hexdigest()


def md5_file(path: pathlib.Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_source(raw: bytes, require_source_md5: bool = True) -> dict:
    if not raw or len(raw) > 20000 or b"\x00" in raw:
        raise ValueError("Original metadata must be a small, nonempty text TSV.")
    actual_md5 = md5_bytes(raw)
    if require_source_md5 and actual_md5 != METADATA_MD5:
        raise ValueError("Zenodo source metadata checksum differs from published v3 MD5.")
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("Metadata is not UTF-8.") from exc
    # DictReader would silently collapse repeated header names in a dict.
    header = next(csv.reader(io.StringIO(text, newline=""), delimiter="\t"), [])
    if len(header) != len(set(header)) or any(not x for x in header):
        raise ValueError("Duplicated or empty source metadata header names.")
    rows = list(csv.DictReader(io.StringIO(text, newline=""), delimiter="\t"))
    if len(rows) != 12 or not rows:
        raise ValueError("Unexpected original study size: need 6 biological, 3 pooled QC, 3 blanks.")
    fields = list(rows[0])
    if any(name not in fields for name in REQUIRED) or len(fields) != len(set(fields)):
        raise ValueError("Original required metadata columns are missing or duplicated.")
    if any(row.get(None) is not None for row in rows):
        raise ValueError("Malformed source TSV with surplus values.")
    sample_names = []
    orders = []
    classes = []
    for row in rows:
        name = (row["sample_name"] or "").strip()
        typ = (row["class"] or "").strip()
        kind = (row["sampleType"] or "").strip()
        if name not in {path.removesuffix(".mzML") for path in RAW_MD5}:
            raise ValueError("Unrecognized source sample identity: " + name)
        if typ not in ROLES or kind != ROLES[typ][0]:
            raise ValueError("Published sample class/type mismatch for " + name)
        if (row["batch"] or "").strip() != "ne1":
            raise ValueError("Unexpected source batch; never invent or collapse real batches.")
        entry = (row["injectionOrder"] or "").strip()
        if not re.fullmatch(r"[1-9][0-9]*", entry):
            raise ValueError("Missing or non-integer source injection order for " + name)
        sample_names.append(name)
        orders.append(int(entry))
        classes.append(ROLES[typ][1])
    if len(set(sample_names)) != 12 or len(set(orders)) != 12:
        raise ValueError("Duplicate sample identity or published injection order.")
    if {name + ".mzML" for name in sample_names} != set(RAW_MD5):
        raise ValueError("Metadata sample IDs differ from the complete 12-file Zenodo v3 manifest.")
    counts = {name: classes.count(name) for name in set(classes)}
    if counts != {"biological": 6, "pooled_qc": 3, "blank": 3}:
        raise ValueError("Biological / QC / blank classes are not as published.")
    return {
        "source": ZENODO_RECORD,
        "doi": DOI,
        "record_version": "v3",
        "original_metadata_name": METADATA_FILE,
        "metadata_md5_published": METADATA_MD5,
        "metadata_md5_observed": actual_md5,
        "metadata_md5_verified": (actual_md5 == METADATA_MD5),
        "all_file_md5_published": RAW_MD5,
        "sample_count": len(rows),
        "roles": {"biological": 6, "pooled_qc": 3, "blank": 3},
        "source_batch": "ne1",
        "source_injection_orders_verified_from_metadata": True,
        "source_injection_order_range": [min(orders), max(orders)],
        "injection_chronology_policy": (
            "Original Zenodo values preserved. These are sparse subset injections "
            "from a larger run, not 1..12 contiguous injections."
        ),
        "qc_drift_claim": "NOT_VALIDATED_three_sparse_pool_QC_injections",
        "raw_peak_processing": "NOT_EXECUTED_source_preflight_only",
        "independent_external_biological_validation": False,
        "clinical_certification": False,
    }


def verify_raw_dir(directory: pathlib.Path) -> dict:
    if not directory.is_dir():
        raise ValueError("Raw directory does not exist.")
    available = [path.name for path in directory.iterdir() if path.is_file()]
    if set(available) != set(RAW_MD5) or len(available) != 12:
        raise ValueError("Original raw input directory must contain exactly 12 published mzML files.")
    for name, expected in RAW_MD5.items():
        if md5_file(directory / name) != expected:
            raise ValueError(f"Raw mzML differs from Zenodo published md5: {name}")
    return {"raw_file_count": 12, "raw_file_md5_verified": True}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--raw-dir", type=pathlib.Path)
    args = parser.parse_args()
    if args.out.exists() and (args.out.is_file() or any(args.out.iterdir())):
        raise SystemExit("Output location is not empty: cannot overwrite prior scientific evidence.")
    result = validate_source(args.metadata.read_bytes())
    if args.raw_dir is not None:
        result.update(verify_raw_dir(args.raw_dir))
    else:
        result.update({"raw_file_count": 0, "raw_file_md5_verified": False})
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "zenodo_sacurine_toy_v3_source_preflight.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("PASS: real Zenodo metadata role/injection-order integrity; raw mzML processing NOT claimed.")


if __name__ == "__main__":
    main()
