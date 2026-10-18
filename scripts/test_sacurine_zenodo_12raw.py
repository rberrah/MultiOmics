#!/usr/bin/env python3
"""Offline adversarial tests, not original human LC-MS experiments."""
import importlib.util
import pathlib
import sys
import tempfile
import types

ROOT = pathlib.Path("scripts").resolve()
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location(
    "sacurine_twelve_original", ROOT / "verify_sacurine_zenodo_12raw.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
from qualify_sacurine_toy_zenodo import RAW_MD5, METADATA_FILE, METADATA_MD5


def fail(f):
    try:
        f()
    except (ValueError, OSError):
        return
    raise AssertionError("Tampered or incomplete source was incorrectly accepted.")


size_map = {name: 20_000_000 for name in RAW_MD5}
size_map[METADATA_FILE] = 645
rec = {"id": m.RECORD_ID, "files": [
    {"key": name, "checksum": "md5:" + checksum, "size": size_map[name]}
    for name, checksum in {**RAW_MD5, METADATA_FILE: METADATA_MD5}.items()
]}
sizes = m.inspect_record(rec)
assert len(sizes) == 13 and sizes[METADATA_FILE] == 645
for changed in [
    {**rec, "id": "wrong-record"},
    {**rec, "files": rec["files"][:-1]},
    {**rec, "files": rec["files"] + [rec["files"][0]]},
    {**rec, "files": [{**rec["files"][0], "checksum": "md5:" + "a" * 32}]
                   + rec["files"][1:]},
    {**rec, "files": [{**rec["files"][0], "size": 200_000_000}]
                   + rec["files"][1:]},
]:
    fail(lambda source=changed: m.inspect_record(source))

with tempfile.TemporaryDirectory() as dirname:
    root = pathlib.Path(dirname)
    missing = root / "example.mzML"
    fail(lambda: m.download_verified(
        "unknown.mzML", missing, 100, "a" * 32, attempts=1))
    missing.write_bytes(b"already exists")
    fail(lambda: m.download_verified(
        next(iter(RAW_MD5)), missing, 100, "a" * 32, attempts=1))

assert sum(name.startswith("Blanc") for name in RAW_MD5) == 3
assert sum(name.startswith("QC1_") for name in RAW_MD5) == 3
assert sum(name.startswith("HU_neg_") for name in RAW_MD5) == 6
assert m.MAX_TOTAL < 400_000_000
assert m.MAX_FILE < 50_000_000
print("PASS: authentic 12-original Spectra Zenodo v3 manifest, publisher MD5 and fail-closed source budget")
