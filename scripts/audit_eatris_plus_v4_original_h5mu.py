#!/usr/bin/env python3
"""Source-original privacy-aware EATRIS-Plus v4 h5mu investigation.

This research file is deliberately a SOURCE METADATA EXPLORATION,
not an RNA/protein prediction benchmark or a clinical PK study.
Downloading a 141.7MB published h5mu is opt-in via GitHub Actions
and no participant-row or assay matrix is redistributed.
"""
from __future__ import annotations
import hashlib
import json
import os
import tempfile
from collections import Counter
from urllib.request import Request,urlopen

RECORD="https://zenodo.org/api/records/19136013"
EXPECTED_NAME="EATRIS_PLUS_mudata.h5mu"
EXPECTED_MD5="bbf0342797adcbcd6a3c5b242ca64573"
EXPECTED_RECORD_ID=19136013
MAX_BYTES=160_000_000

class EATRISSourceError(ValueError): pass

def _source_metadata():
    req=Request(RECORD,headers={"User-Agent":"MultiOmics-original-EATRIS-v4-qualification/1.0",
                                 "Accept":"application/json"})
    with urlopen(req,timeout=30) as resp:
        meta=json.load(resp)
    if meta.get("id")!=EXPECTED_RECORD_ID:
        raise EATRISSourceError("Zenodo original v4 record ID drift")
    files=meta.get("files") or []
    matches=[f for f in files if f.get("key")==EXPECTED_NAME]
    if len(matches)!=1:
        raise EATRISSourceError("Expected original MuData h5mu file not in v4")
    f=matches[0];url=f.get("links",{}).get("self")
    if (f.get("checksum")!="md5:"+EXPECTED_MD5 or not isinstance(url,str)
            or not url.startswith("https://")):
        raise EATRISSourceError("Zenodo file content hash or https download link mismatch")
    size=f.get("size")
    if not isinstance(size,int) or size<100_000_000 or size>MAX_BYTES:
        raise EATRISSourceError("Official h5mu original size inconsistent")
    return url,size

def download_original(path):
    url,size=_source_metadata()
    m=hashlib.md5();total=0
    req=Request(url,headers={"User-Agent":"MultiOmics-source-v4-MuData-replay/1.0"})
    with urlopen(req,timeout=80) as response,open(path,"wb") as dst:
        while True:
            chunk=response.read(2_000_000)
            if not chunk:break
            total+=len(chunk)
            if total>MAX_BYTES:raise EATRISSourceError("Original source exceeds budget")
            m.update(chunk);dst.write(chunk)
    if total!=size or m.hexdigest()!=EXPECTED_MD5:
        raise EATRISSourceError("Original MuData bytes failed published size/hash check")
    return total

def audit_h5mu(path):
    import h5py
    # EATRIS source structure includes measurement rows: emit schema
    # only, NEVER variable, gene, phenotype or subject identifiers.
    with h5py.File(path,"r") as f:
        def describe(g,depth=0):
            d={}
            if depth>=3:return d
            for key,item in g.items():
                if not isinstance(item,(h5py.Group,h5py.Dataset)):continue
                # Expose public structure of h5mu HDF5 format but no
                # sample-specific labels, donor IDs or source measurement.
                if isinstance(item,h5py.Dataset):
                    d[key]={"kind":"dataset","shape":list(item.shape),
                            "dtype_family":item.dtype.kind}
                else:
                    d[key]={"kind":"group","sub":describe(item,depth+1)}
            return d
        if "mod" not in f or not isinstance(f["mod"],h5py.Group):
            raise EATRISSourceError("Official data has no MuData mod group")
        modalities={}
        for key,g in f["mod"].items():
            # These are original public source modality names, not
            # individual data or identifiable metadata.
            modalities[key]={
                "HDF5_group_schema":describe(g,1),
                "encoding_type":str(g.attrs.get("encoding-type","UNKNOWN")),
            }
        return {
            "status":"ORIGINAL_EATRIS_PLUS_V4_H5MU_AUTHENTICATED_SCHEMA_ONLY",
            "zenodo_record_id":EXPECTED_RECORD_ID,
            "published_version":"v4",
            "source_file_name":EXPECTED_NAME,
            "source_checksum_MD5":EXPECTED_MD5,
            "n_original_modality_groups":len(modalities),
            "original_modality_group_schema":modalities,
            "overall_h5mu_group_top_level":sorted(f.keys()),
            "original_individual_donor_labels_or_values_emitted":False,
            "clinical_PK_or_germline_PGx_measured_in_open_MuData":False,
        }

def original_source():
    with tempfile.TemporaryDirectory() as tmp:
        p=os.path.join(tmp,"original.h5mu")
        bytes_count=download_original(p)
        report=audit_h5mu(p)
        report["source_bytes_verified"]=bytes_count
        return report

if __name__=="__main__":
    try:print(json.dumps(original_source(),indent=2,sort_keys=True))
    except (EATRISSourceError,OSError,TimeoutError,ValueError) as exc:
        print(json.dumps({"status":"ORIGINAL_EATRIS_V4_SOURCE_UNAVAILABLE_ABSTAIN",
                          "error_type":type(exc).__name__,
                          "source_original_qualified":False}))
        raise SystemExit(2)
