#!/usr/bin/env python3
"""Independent source-availability probe for MTBLS404 ISA metadata only.

No raw study data downloaded, and no individual-level records archived.
Use only as source-discovery evidence; not chromatographic validation.
"""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

SOURCES = (
    "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS404/s_sacurine.txt",
    "https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS404/s_sacurine.txt",
)


def run():
    output = Path("tmp/mtbls404-metadata-probe")
    output.mkdir(parents=True, exist_ok=True)
    record = {
        "schema": "multiomics.mtbls404-primary-metadata-source-check.v1",
        "source_sha": os.environ.get("GITHUB_SHA", "local"),
        "candidate_urls": list(SOURCES),
        "raw_data_retrieved": False,
        "xcms_executed": False,
        "status": "no_verified_source"
    }
    try:
        evidence = None
        for url in SOURCES:
            for attempt in range(2):
                try:
                    request = Request(url, headers={"User-Agent": "MultiOmics-External-Validation/1.0"})
                    with urlopen(request, timeout=32) as response:
                        result = response.read(5_000_001)
                        if len(result) > 5_000_000 or not result:
                            raise ValueError("Invalid or oversized source metadata")
                        body = result.decode("utf-8-sig", errors="strict")
                        reader = csv.DictReader(io.StringIO(body), delimiter="\t")
                        cols = reader.fieldnames or []
                        if len(cols) < 3:
                            raise ValueError("ISA-Tab study source does not have enough columns")
                        rows = list(reader)
                        if len(rows) < 20:
                            raise ValueError("Unexpectedly few study records")
                        names = [s for s in cols if ("sample" in s.lower() or
                             "subset" in s.lower() or "comment" in s.lower() or
                             "type" in s.lower())]
                        evidence = {
                            "retrieved_from": response.url,
                            "metadata_bytes": len(result),
                            "metadata_sha256": hashlib.sha256(result).hexdigest(),
                            "metadata_columns": cols,
                            "record_count": len(rows),
                            "candidate_identification_columns": names,
                            "contains_subset_indicator": any("subset" in s.lower() for s in cols),
                            "sample_records_archived": False,
                        }
                        break
                except (TimeoutError, HTTPError, URLError, UnicodeError, ValueError) as error:
                    record.setdefault("attempts", []).append({
                        "url": url, "error_type": type(error).__name__,
                        "http_status": getattr(error, "code", None)
                    })
                    if attempt == 0:
                        time.sleep(1)
            if evidence:
                break
        if not evidence:
            raise RuntimeError("Primary EBI study metadata is not accessible from this runner")
        record["observed_source"] = evidence
        record["status"] = "primary_metadata_available_but_no_raw_sacurine_validation"
        print("PASS: official MetaboLights study metadata accessible",
              evidence["record_count"], "source records",
              evidence["metadata_sha256"])
    finally:
        (output / "metabolights_source_probe.json").write_text(
            json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
    if record["status"] == "no_verified_source":
        raise RuntimeError("EBI primary metadata source unavailable; see archived evidence")


if __name__ == "__main__":
    run()
