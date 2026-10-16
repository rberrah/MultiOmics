#!/usr/bin/env python3
"""Verified Zenodo/Galaxy Sacurine real pooled-QC training mzML acquisitions.

Galaxy GTN references Zenodo record 3244991 for exactly 6 original urine
acquisitions and 3 genuine QC injections. Require record API-published MD5,
record version and per-file SHA256; no silent unverified raw source.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

BASE="https://zenodo.org"
RECORD="3244991"
RAW_NAMES=[
 "HU_neg_048.mzML","HU_neg_090.mzML","HU_neg_123.mzML",
 "HU_neg_157.mzML","HU_neg_173.mzML","HU_neg_192.mzML",
 "QC1_002.mzML","QC1_008.mzML","QC1_014.mzML",
]
MAX_JSON=2_000_000
MAX_FILE=60_000_000

def remote_record():
    url=f"{BASE}/api/records/{RECORD}"
    request=Request(url,headers={"User-Agent":"MultiOmics-External-LCMS-Study/1.0"})
    with urlopen(request,timeout=75) as response:
        data=response.read(MAX_JSON+1)
        if not 50<len(data)<=MAX_JSON:
            raise ValueError("Invalid Zenodo API record response size")
    parsed=json.loads(data)
    if str(parsed.get("id"))!=RECORD:
        raise ValueError("Zenodo did not return exact immutable scientific record")
    by_key={}
    for entry in parsed.get("files",[]):
        name=entry.get("key") or entry.get("filename")
        checksum=entry.get("checksum","")
        if not name or not isinstance(checksum,str):
            continue
        if name in by_key:raise ValueError("Duplicate Zenodo record file name")
        by_key[name]=entry
    missing=sorted(set(RAW_NAMES)-set(by_key))
    if missing:
        raise ValueError("Galaxy training public file absent from Zenodo record: "+",".join(missing))
    return by_key

def verify_copy(url,destination,md5):
    pending=destination.with_name(destination.name+".incomplete")
    pending.unlink(missing_ok=True)
    try:
        req=Request(url,headers={"User-Agent":"MultiOmics-External-LCMS-Study/1.0"})
        with urlopen(req,timeout=95) as remote,pending.open("wb") as output:
            hmd5=hashlib.md5()
            hsha=hashlib.sha256()
            n=0
            for buff in iter(lambda:remote.read(1024*1024),b""):
                n+=len(buff)
                if n>MAX_FILE:raise ValueError("File exceeds fixed raw-study safety cap")
                hmd5.update(buff)
                hsha.update(buff)
                output.write(buff)
        if n<1_000_000 or hmd5.hexdigest()!=md5:
            raise ValueError("Real Zenodo raw sample MD5 mismatch or truncation")
        pending.replace(destination)
        return dict(bytes=n,md5=md5,sha256=hsha.hexdigest())
    finally:
        pending.unlink(missing_ok=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("tmp/sacurine-zenodo-qc-source"))
    parser.add_argument("--metadata-only",action="store_true",help="Only query official published raw file sizes/checksums; do not download any spectra")
    args=parser.parse_args()
    out=args.output
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Refusing to reuse or overwrite real source folder")
    report=dict(schema="multiomics.real-zenodo-pooled-qc-source.v1",
                record_id=RECORD,
                record_url=BASE+"/records/"+RECORD,
                source="Galaxy Training Network LC-MS tutorial",
                expected_real_biological=6,expected_real_pooled_qc=3,expected_blanks=0,
                source_sha=os.getenv("QUALIFIED_SOURCE_SHA",os.getenv("GITHUB_SHA","local")),
                raw_processing_completed=False,
                data_provenance_verified=False,
                file_evidence=[])
    try:
        entries=remote_record()
        if args.metadata_only:
            originals=[]
            for name in RAW_NAMES:
                e=entries[name]
                originals.append(dict(name=name,reported_bytes=e.get("size"),published_checksum=e.get("checksum")))
            report["published_file_inventory"]=originals
            report["reported_total_bytes"]=sum(int(e["reported_bytes"]) for e in originals)
            report["status"]="zenodo_api_original_file_metadata_verified_only"
            print(json.dumps({"total_source_bytes":report["reported_total_bytes"],
                              "file_sizes":originals},indent=2),flush=True)
            return
        for name in RAW_NAMES:
            entry=entries[name]
            check=entry.get("checksum","")
            if not check.startswith("md5:") or len(check)!=36:
                raise ValueError("No original Zenodo MD5 for "+name)
            published_md5=check[4:]
            # Zenodo checksum and exact immutable record ID are source authority.
            url=f"{BASE}/records/{RECORD}/files/{name}?download=1"
            info=None
            for attempt in range(3):
                try:
                    info=verify_copy(url,out/name,published_md5)
                    break
                except (HTTPError,URLError,TimeoutError) as err:
                    if attempt==2:raise RuntimeError("Public Zenodo source inaccessible: "+name) from err
                    time.sleep(2+attempt)
            report["file_evidence"].append(dict(name=name,
                      role="pooled_qc" if name.startswith("QC1_") else "biological",
                      **info))
            print("verified original pooled-QC study file:",name,flush=True)
        report["data_provenance_verified"]=True
        report["status"]="9_true_raw_files_zenodo_md5_checked"
    finally:
        (out/"zenodo_galaxy_true_qc_source_manifest.json").write_text(
           json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("PASS: six original urine + three original pooled QC spectra verified")

if __name__=="__main__":
    main()
