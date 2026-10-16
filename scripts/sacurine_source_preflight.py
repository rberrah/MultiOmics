#!/usr/bin/env python3
"""Sacurine MTBLS404 external raw-data preflight.

Do not trust remote ZIP names or infer chromatographic acquisition order.
No publication-grade acceptance until the archive SHA is separately pinned.
The archive is never committed or uploaded as a GitHub Actions artifact.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from zipfile import ZipFile, BadZipFile

ORIGIN = "https://download.workflow4metabolomics.org/datasets/"
ARCHIVE = ORIGIN + "sacurine-neg-subset_mzXML.zip"
METADATA = ORIGIN + "sacurine-neg_sampleMetadata.tsv"
MAX_ARCHIVE_BYTES = 300_000_000
MAX_EXTRACTED_BYTES = 2_000_000_000
EXPECTED = {"biological": 20, "pooled_qc": 10, "blank": 6}
SUPPORT_FILE_EXTENSIONS = {".mzxml", ".mzml"}
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def download(url: str, target: Path, limit: int) -> str:
    """Strict upper bound, HTTP(S) only, no partial files retained as success."""
    if not url.startswith("https://"):
        raise ValueError("Only verified research HTTPS sources are allowed")
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".incomplete")
    for attempt in range(3):
        digest = hashlib.sha256()
        written = 0
        try:
            request = Request(url, headers={"User-Agent": "MultiOmics-Research-Qualification/1.0"})
            with urlopen(request, timeout=150) as response, partial.open("wb") as output:
                if not response.url.startswith("https://"):
                    raise ValueError("Unencrypted research data redirect")
                for buf in iter(lambda: response.read(1024 * 1024), b""):
                    written += len(buf)
                    if written > limit:
                        raise ValueError("Remote scientific file exceeds fixed size budget")
                    digest.update(buf)
                    output.write(buf)
            if written == 0:
                raise ValueError("Empty scientific source file")
            partial.replace(target)
            return digest.hexdigest()
        except (HTTPError, URLError, TimeoutError) as exc:
            partial.unlink(missing_ok=True)
            if attempt == 2:
                raise RuntimeError(f"Remote research data unavailable after 3 attempts: {url}") from exc
            time.sleep(2 ** attempt)
    raise AssertionError("Unreachable")


def identify_role(basename: str) -> str:
    # These are provisional, explicitly validated again with source sample sheet.
    stem = Path(basename).stem
    if re.match(r"(?i)^blanc([0-9_]|$)", stem):
        return "blank"
    if re.match(r"(?i)^qc([0-9_]|$)", stem):
        return "pooled_qc"
    if re.match(r"(?i)^hu[_-]", stem):
        return "biological"
    raise ValueError(f"Cannot safely identify sample purpose from filename: {basename}")


def inspect_archive(archive: Path) -> tuple[list[dict], list[str]]:
    rows: list[dict] = []
    other: list[str] = []
    total = 0
    seen: set[str] = set()
    with ZipFile(archive) as zipfile:
        for entry in zipfile.infolist():
            if entry.is_dir():
                continue
            name = entry.filename
            path = PurePosixPath(name)
            if path.is_absolute() or ".." in path.parts or "\\" in name:
                raise ValueError(f"Unsafe ZIP path: {name}")
            if any(part.startswith(".") for part in path.parts):
                raise ValueError(f"Unexpected hidden ZIP entry: {name}")
            # Reject symlink and arbitrary stored executable objects.
            mode = (entry.external_attr >> 16) & 0o170000
            if mode not in (0, 0o100000):
                raise ValueError(f"Unexpected ZIP entry type: {name}")
            total += entry.file_size
            if total > MAX_EXTRACTED_BYTES:
                raise ValueError("Uncompressed archive exceeds research-data budget")
            if path.suffix.lower() not in SUPPORT_FILE_EXTENSIONS:
                other.append(name)
                continue
            short = path.name
            if not NAME_RE.fullmatch(short):
                raise ValueError(f"Unsupported raw scientific filename: {short}")
            key = short.casefold()
            if key in seen:
                raise ValueError(f"Duplicate / case-colliding chromatogram: {short}")
            seen.add(key)
            rows.append({
                "filename": short,
                "zip_member": name,
                "role_provisional": identify_role(short),
                "compressed_bytes": entry.compress_size,
                "uncompressed_bytes": entry.file_size,
                "crc32": f"{entry.CRC:08x}"
            })
    counts = {role: sum(r["role_provisional"] == role for r in rows) for role in EXPECTED}
    if counts != EXPECTED or len(rows) != 36:
        raise ValueError(f"Unexpected raw sample balance {counts}; expected {EXPECTED}")
    return sorted(rows, key=lambda r: r["filename"].lower()), other


def inspect_metadata(meta_path: Path, rows: list[dict]) -> dict:
    with meta_path.open(newline="", encoding="utf-8-sig", errors="replace") as stream:
        head = stream.readline()
        if "\t" not in head:
            raise ValueError("Official metadata does not appear to be TSV")
        stream.seek(0)
        reader = csv.DictReader(stream, delimiter="\t")
        header = reader.fieldnames or []
        matrix = list(reader)
    sample_candidates = [
        key for key in header if re.search(r"sample|file|assay", key, re.IGNORECASE)
    ]
    role_candidates = [
        key for key in header if re.search(r"class|type|group", key, re.IGNORECASE)
    ]
    if not matrix or not header:
        raise ValueError("Official sample metadata is blank")
    stemlist = {Path(r["filename"]).stem.casefold() for r in rows}
    matched = {}
    for key in sample_candidates:
        val = {Path(str(r.get(key, "")).strip()).stem.casefold() for r in matrix}
        matched[key] = len(stemlist & val)
    return {
        "schema": "multiomics.sacurine-metadata-preflight.v1",
        "metadata_columns": header,
        "metadata_rows": len(matrix),
        "source_name_candidates": sample_candidates,
        "role_candidates": role_candidates,
        "exact_sample_name_matches": matched,
        "official_roles_confirmed": False,
        "comment": "Filename classes remain provisional until raw roles and complete source metadata are reconciled. Nothing here verifies QC injection sequence."
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("tmp/sacurine-preflight"))
    args = parser.parse_args()
    out = args.output
    if out.exists() and any(out.iterdir()):
        raise RuntimeError("Refusing to overwrite prior Sacurine run artifacts")
    out.mkdir(parents=True, exist_ok=True)
    private_archive = Path(os.environ.get("RUNNER_TEMP", str(out))) / "sacurine-neg-subset_mzXML.zip"
    private_metadata = out / "sacurine-neg_sampleMetadata.tsv"
    state = {
        "schema": "multiomics.sacurine-raw-external-preflight.v1",
        "study": "MetaboLights MTBLS404; Workflow4Metabolomics Sacurine",
        "reference": "Thévenot et al. 2015, DOI:10.1021/acs.jproteome.5b00354",
        "archive_url": ARCHIVE,
        "metadata_url": METADATA,
        "expected_sample_roles": EXPECTED,
        "raw_file_sha256_lock": "not yet established; first source qualification",
        "source_commit": os.environ.get("GITHUB_SHA", "local"),
        "raw_peak_picking_executed": False,
        "scientific_certification": False,
    }
    try:
        state["sample_metadata_sha256"] = download(METADATA, private_metadata, 8_000_000)
        state["archive_sha256"] = download(ARCHIVE, private_archive, MAX_ARCHIVE_BYTES)
        state["archive_size_bytes"] = private_archive.stat().st_size
        rows, others = inspect_archive(private_archive)
        state["sample_roles_provisional"] = EXPECTED
        state["unexpected_nondata_entries"] = others
        state["metadata_assessment"] = inspect_metadata(private_metadata, rows)
        with (out / "archive_chromatogram_index.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        state["status"] = "source_acquired_and_archive_structure_valid"
        print("PASS: official public Sacurine archive structure and expected 20 biological / 10 QC / 6 blanks", flush=True)
        print("ARCHIVE SHA256:", state["archive_sha256"], "bytes:", state["archive_size_bytes"], flush=True)
        print("SOURCE METADATA COLUMNS:", state["metadata_assessment"]["metadata_columns"], flush=True)
        print("SAMPLE JOIN MATCHES:", state["metadata_assessment"]["exact_sample_name_matches"], flush=True)
    except Exception as exc:
        state["status"] = "failed_preflight"
        state["failure_type"] = type(exc).__name__
        state["failure_summary"] = str(exc)
        raise
    finally:
        (out / "sacurine_source_provenance.json").write_text(
            json.dumps(state, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        # The actual data remain private to the CI machine and never become
        # a GH Actions artifact, even on failure.
        private_archive.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
