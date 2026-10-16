#!/usr/bin/env python3
"""Acquisition audited of public Sacurine toy raw files, Zenodo record 1346742.

Independent real-source test for blank injections. These are 2 biological
urine specimens + 2 water blanks, NOT pooled QC, NOT a population study.
Original mzXML files are not committed or exported from CI as artifacts.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

RECORD="1346742"
BASE="https://zenodo.org/records/1346742/files/"
SOURCES={
 "Blanc12.mzXML":("blank","bd61be730ad155c81c9dab2725e09f86"),
 "Blanc17.mzXML":("blank","985839a47c0a3643639439f8e9a52754"),
 "HU_neg_099.mzXML":("biological","69c5e3c80588d66a7faa2bb4f2649903"),
 "HU_neg_178.mzXML":("biological","967aa2fae55ff9f9b9a72295aacf669e"),
}
SOURCE_MD5="0097d890313953a52a86d00360a1cf81"
MAX_BYTES=18_000_000

def inspect_local_file(path, expected_md5):
    if not path.is_file() or not 1_000_000 < path.stat().st_size < MAX_BYTES:
        raise ValueError("Unexpected or missing original mzXML file: "+path.name)
    md5=hashlib.md5()
    sha=hashlib.sha256()
    with path.open("rb") as inp:
        head=inp.read(2048)
        if b"mzXML" not in head and b"mzxml" not in head:
            raise ValueError("Source does not appear to be the expected mzXML")
        inp.seek(0)
        for chunk in iter(lambda:inp.read(1024*1024),b""):
            md5.update(chunk)
            sha.update(chunk)
    if md5.hexdigest()!=expected_md5:
        raise ValueError("Research-provider MD5 changed: "+path.name)
    return {"bytes":path.stat().st_size,"published_md5":expected_md5,
            "observed_sha256":sha.hexdigest()}

def acquire_one(url, out, expected_md5):
    for attempt in range(3):
        pending=out.with_name(out.name+".partial")
        pending.unlink(missing_ok=True)
        try:
            request=Request(url,headers={"User-Agent":"MultiOmics-Scientific-Replay/1.0"})
            with urlopen(request,timeout=65) as source, pending.open("wb") as sink:
                if not source.url.startswith("https://"):
                    raise ValueError("Unexpected insecure redirect for public LC-MS archive")
                written=0
                for buf in iter(lambda:source.read(512*1024),b""):
                    written+=len(buf)
                    if written>=MAX_BYTES:
                        raise ValueError("Unexpected oversize Sacurine source")
                    sink.write(buf)
            pending.replace(out)
            return inspect_local_file(out,expected_md5)
        except (HTTPError,URLError,TimeoutError) as e:
            pending.unlink(missing_ok=True)
            if attempt==2:
                raise RuntimeError("Public Zenodo Sacurine file inaccessible: "+out.name) from e
            time.sleep(attempt+1)
    raise RuntimeError("unreachable")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("tmp/sacurine-zenodo-toy"))
    parser.add_argument("--local-mzxml-directory",type=Path,default=None)
    args=parser.parse_args()
    folder=args.output
    folder.mkdir(parents=True,exist_ok=True)
    report={
      "schema":"multiomics.external-sacurine-zenodo-raw.v1",
      "source":"https://zenodo.org/records/1346742",
      "doi":"10.5281/zenodo.1346742",
      "published_md5_authority":"Zenodo public 1346742 record file checksums",
      "source_commit":os.environ.get("QUALIFIED_SOURCE_SHA",os.environ.get("GITHUB_SHA","local")),
      "samples":[{"filename":k,"sample_type":v[0]} for k,v in SOURCES.items()],
      "research_design":"2 real biological urine + 2 true blanks; zero pooled QC",
      "native_xcms_executed":False,
      "statistically_validated":False,
      "status":"unverified"
    }
    try:
        for name,(kind,expected_md5) in SOURCES.items():
            source=folder/name if args.local_mzxml_directory is None else args.local_mzxml_directory/name
            evidence=(acquire_one(BASE+name+"?download=1",source,expected_md5)
                      if args.local_mzxml_directory is None else inspect_local_file(source,expected_md5))
            row=next(x for x in report["samples"] if x["filename"]==name)
            row.update(evidence)
            print("VERIFIED",name,evidence["observed_sha256"],flush=True)
        report["status"]="four_real_raw_files_integrity_verified"
    finally:
        (folder/"sacurine_zenodo_source_manifest.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("PASS: four real public Sacurine mzXML files match original provider checksums",flush=True)

if __name__=="__main__":
    main()
