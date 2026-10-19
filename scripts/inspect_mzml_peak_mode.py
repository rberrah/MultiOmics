#!/usr/bin/env python3
"""Inspect mzML MS1 peak representation without interpreting downstream yield."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import xml.etree.ElementTree as ET

CENTROID="MS:1000127"
PROFILE="MS:1000128"
MS_LEVEL="MS:1000511"

def local(tag:str)->str:
    return tag.rsplit("}",1)[-1]

def inspect(path:Path)->dict:
    groups={}
    counts={"centroid":0,"profile":0,"unknown":0}
    ms1=0
    for event,elem in ET.iterparse(path,events=("end",)):
        tag=local(elem.tag)
        if tag=="referenceableParamGroup":
            ident=elem.attrib.get("id","")
            accessions={x.attrib.get("accession","") for x in elem.iter()
                        if local(x.tag)=="cvParam"}
            if ident: groups[ident]=accessions
        elif tag=="spectrum":
            accessions=set()
            ms_level=None
            for x in elem.iter():
                xt=local(x.tag)
                if xt=="cvParam":
                    acc=x.attrib.get("accession","")
                    accessions.add(acc)
                    if acc==MS_LEVEL: ms_level=x.attrib.get("value")
                elif xt=="referenceableParamGroupRef":
                    accessions.update(groups.get(x.attrib.get("ref",""),set()))
            if ms_level=="1":
                ms1+=1
                centroid=CENTROID in accessions
                profile=PROFILE in accessions
                if centroid==profile:
                    counts["unknown"]+=1
                elif centroid:
                    counts["centroid"]+=1
                else:
                    counts["profile"]+=1
            elem.clear()
    if ms1<1:
        raise ValueError(f"No MS1 spectra found in {path.name}")
    represented=[key for key in ("centroid","profile") if counts[key]>0]
    if counts["unknown"]>0 or len(represented)!=1:
        raise ValueError(
            f"Ambiguous/mixed MS1 peak representation in {path.name}: {counts}"
        )
    return {"file":path.name,"ms1_spectra":ms1,"mode":represented[0],
            "counts":counts}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--raw-dir",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    files=sorted(args.raw_dir.glob("*.mzML"))
    if len(files)!=12:
        raise SystemExit(f"Expected exactly 12 mzML files, found {len(files)}")
    rows=[inspect(p) for p in files]
    modes={x["mode"] for x in rows}
    if len(modes)!=1:
        raise SystemExit("The 12 authenticated runs do not share one peak representation: "+repr(modes))
    report={
        "schema":"multiomics.sacurine-openms-peak-mode.v1",
        "common_mode":next(iter(modes)),
        "files":rows,
        "decision":"PeakPickerHiRes required only when common_mode == profile",
        "scientific_certification":False
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("PASS:",report["common_mode"],"MS1 representation across 12 authenticated mzML")

if __name__=="__main__":
    main()
