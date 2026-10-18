#!/usr/bin/env python3
"""Acquire the authentic 12 Sacurine v3 mzML injections from fixed Zenodo record.

Strictly source verification and on-disk acquisition. No sample-level original
chromatograms or clinical attributes are exported as GitHub Actions artifacts.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

from qualify_sacurine_toy_zenodo import (
    RAW_MD5, METADATA_FILE, METADATA_MD5, validate_source, md5_file,
)

RECORD_ID = 10130758
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"
BASE = f"https://zenodo.org/records/{RECORD_ID}/files"
MAX_JSON = 2_000_000
MAX_FILE = 45_000_000
MAX_TOTAL = 330_000_000
TIMEOUT = 100


def inspect_record(record: dict) -> dict[str, int]:
    if str(record.get("id")) != str(RECORD_ID):
        raise ValueError("Zenodo returned the wrong immutable v3 record.")
    entries = {}
    for row in record.get("files", []):
        name = row.get("key") or row.get("filename")
        if not isinstance(name, str) or name in entries:
            raise ValueError("Missing/duplicate official Zenodo record filename.")
        entries[name] = row
    expected = dict(RAW_MD5)
    expected[METADATA_FILE] = METADATA_MD5
    sizes = {}
    for name, md5 in expected.items():
        row = entries.get(name)
        if not isinstance(row, dict) or row.get("checksum") != "md5:" + md5:
            raise ValueError("Original publisher MD5 changed for " + name)
        size = row.get("size")
        if not isinstance(size, int) or size < 1 or size > MAX_FILE:
            raise ValueError("Original publisher size is absent/invalid for " + name)
        sizes[name] = size
    if sizes[METADATA_FILE] != 645 or sum(sizes.values()) > MAX_TOTAL:
        raise ValueError("Unexpected source size or research download budget.")
    return sizes


def fetch_bytes(url: str, max_bytes: int) -> bytes:
    req = Request(url, headers={"User-Agent": "MultiOmics-research-source-replay/1.0"})
    with urlopen(req, timeout=TIMEOUT) as response:
        if not response.url.startswith("https://"):
            raise ValueError("Unsafe redirect in public science source.")
        data = response.read(max_bytes + 1)
    if not data or len(data) > max_bytes:
        raise ValueError("External scientific file exceeds the prespecified bound.")
    return data


def download_verified(name: str, dest: Path, expected_size: int,
                      expected_md5: str, attempts: int = 3) -> dict:
    if name not in RAW_MD5:
        raise ValueError("Source filename is not part of fixed study.")
    if dest.exists():
        raise ValueError("Refusing to reuse unverified or stale raw LC-MS input.")
    url = BASE + "/" + name + "?download=1"
    partial = dest.with_name(dest.name + ".partial")
    for attempt in range(attempts):
        partial.unlink(missing_ok=True)
        h_md5 = hashlib.md5()
        h_sha = hashlib.sha256()
        size = 0
        try:
            req = Request(url, headers={"User-Agent": "MultiOmics-research-source-replay/1.0"})
            with urlopen(req, timeout=TIMEOUT) as response, partial.open("wb") as sink:
                if not response.url.startswith("https://"):
                    raise ValueError("Unsafe raw-source redirect.")
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    size += len(chunk)
                    if size > MAX_FILE or size > expected_size:
                        raise ValueError("Unexpected raw-file size / changed publisher data.")
                    h_md5.update(chunk)
                    h_sha.update(chunk)
                    sink.write(chunk)
            if size != expected_size or h_md5.hexdigest() != expected_md5:
                raise ValueError("Real raw data differ from fixed Zenodo size or MD5: " + name)
            # Atomic promotion only after complete source checksum verification.
            partial.replace(dest)
            return {"name": name, "bytes": size,
                    "published_md5": expected_md5, "sha256": h_sha.hexdigest(),
                    "role": ("blank" if name.startswith("Blanc") else
                             "pooled_qc" if name.startswith("QC1_") else "biological")}
        except (HTTPError, URLError, TimeoutError):
            partial.unlink(missing_ok=True)
            if attempt == attempts - 1:
                raise
            time.sleep(attempt + 1)
        finally:
            partial.unlink(missing_ok=True)
    raise RuntimeError("Never continue without verified published raw source.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--local-raw-directory", type=Path,
                        help="Offline institutionally acquired originals, still verified via fixed publisher MD5.")
    parser.add_argument("--local-metadata", type=Path)
    args = parser.parse_args()
    if (args.local_raw_directory is None) != (args.local_metadata is None):
        parser.error("Local original spectra and local original metadata must be supplied together.")
    output = args.output
    if output.exists() and any(output.iterdir()):
        raise SystemExit("Refusing to overwrite existing source evidence.")
    output.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": "multiomics.sacurine-zenodo-v3-12raw-provenance.v1",
        "record": f"https://zenodo.org/records/{RECORD_ID}",
        "source_head": os.environ.get("QUALIFIED_SOURCE_SHA",
                        os.environ.get("GITHUB_SHA", "local_unverified")),
        "study": "6 biological urine, 3 pooled QC, 3 true analytic blanks; sparse original injection order",
        "raw_count": 0,
        "metadata_verified": False,
        "raw_source_verified": False,
        "xcms_processing_completed": False,
        "biological_validation": False,
        "clinical_certification": False,
        "files": [],
    }
    try:
        if args.local_raw_directory is not None:
            # Size is not sufficient as authentication; published hashes still required.
            remote_sizes = None
            raw_metadata = args.local_metadata.read_bytes()
        else:
            response = fetch_bytes(API_URL, MAX_JSON)
            remote_sizes = inspect_record(json.loads(response))
            raw_metadata = fetch_bytes(BASE + "/" + METADATA_FILE + "?download=1", 20000)
        source = validate_source(raw_metadata)
        if hashlib.md5(raw_metadata).hexdigest() != METADATA_MD5:
            raise ValueError("Published 12-sample source TSV checksum changed.")
        (output / METADATA_FILE).write_bytes(raw_metadata)
        manifest["metadata_verified"] = True
        manifest["source_injection_orders_verified"] = source[
            "source_injection_orders_verified_from_metadata"]
        for name, md5 in RAW_MD5.items():
            if args.local_raw_directory is not None:
                incoming = args.local_raw_directory / name
                if not incoming.is_file() or incoming.stat().st_size > MAX_FILE:
                    raise ValueError("Missing or oversized original local source: " + name)
                verified_md5 = md5_file(incoming)
                if verified_md5 != md5:
                    raise ValueError("Local raw file mismatches fixed Zenodo MD5: " + name)
                path = output / name
                # Avoid copying 257MB of raw clinical spectra between managed folders.
                path.symlink_to(incoming.resolve())
                row = dict(name=name, bytes=incoming.stat().st_size,
                           published_md5=md5, sha256=_sha256(incoming),
                           role=_role(name))
            else:
                row = download_verified(name, output / name,
                                        remote_sizes[name], md5)
            manifest["files"].append(row)
            manifest["raw_count"] += 1
            print("Verified original research mzML:", name, flush=True)
        if manifest["raw_count"] != 12:
            raise ValueError("Incomplete source acquisition.")
        manifest["raw_source_verified"] = True
        manifest["status"] = "12_original_mzml_publisher_md5_verified"
    finally:
        (output / "sacurine_v3_12_original_source_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")
    print("PASS: 12 authentic Sacurine mzML sources, 6:3:3 roles, published source MD5 and chronology; no xcms claims", flush=True)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _role(name: str) -> str:
    return ("blank" if name.startswith("Blanc") else
            "pooled_qc" if name.startswith("QC1_") else "biological")


if __name__ == "__main__":
    main()
