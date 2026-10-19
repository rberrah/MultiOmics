#!/usr/bin/env python3
"""Compare frozen xcms feature centroids with independent OpenMS consensus centroids."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import xml.etree.ElementTree as ET

PPM=20.0
RT_SECONDS=30.0

def read_xcms(path: Path):
    with path.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    required={"feature_id","mzmed","rtmed"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("xcms centroid CSV lacks feature_id/mzmed/rtmed")
    out=[]
    ids=set()
    for r in rows:
        fid=r["feature_id"].strip()
        mz=float(r["mzmed"]); rt=float(r["rtmed"])
        if not fid or fid in ids or not math.isfinite(mz) or mz<=0 or not math.isfinite(rt) or rt<0:
            raise ValueError("Invalid xcms centroid row")
        ids.add(fid); out.append((fid,mz,rt))
    return out

def local(tag: str)->str:
    return tag.rsplit("}",1)[-1]

def read_openms_consensus(path: Path):
    root=ET.parse(path).getroot()
    out=[]; ids=set()
    for node in root.iter():
        if local(node.tag)!="consensusElement":
            continue
        fid=str(node.attrib.get("id","")).strip()
        centroid=None
        for child in node:
            if local(child.tag)=="centroid":
                centroid=child; break
        if centroid is None:
            raise ValueError("OpenMS consensusElement without centroid")
        mz=float(centroid.attrib["mz"]); rt=float(centroid.attrib["rt"])
        if not fid or fid in ids or not math.isfinite(mz) or mz<=0 or not math.isfinite(rt) or rt<0:
            raise ValueError("Invalid OpenMS consensus centroid")
        ids.add(fid); out.append((fid,mz,rt))
    if not out: raise ValueError("No OpenMS consensus centroids found")
    return out

def candidates(a,b):
    ans=[]
    for _,mz,rt in a:
        hits=[]
        for j,(_,mz2,rt2) in enumerate(b):
            ppm=abs(mz2-mz)/mz*1e6
            if ppm<=PPM and abs(rt2-rt)<=RT_SECONDS:
                hits.append(j)
        ans.append(hits)
    return ans

def compare(xcms,openms):
    fwd=candidates(xcms,openms); rev=candidates(openms,xcms)
    rows=[]
    for i,(fid,mz,rt) in enumerate(xcms):
        hits=fwd[i]
        status="no_coordinate_candidate"
        oid=""; ppm=""; drt=""
        if len(hits)>1:
            status="ambiguous_multiple_candidates"
        elif len(hits)==1:
            j=hits[0]
            if len(rev[j])==1 and rev[j][0]==i:
                status="reciprocal_unique"
                oid,omz,ort=openms[j]
                ppm=(omz-mz)/mz*1e6
                drt=ort-rt
            else:
                status="not_reciprocally_unique"
        rows.append(dict(
            xcms_feature=fid,openms_feature=oid,status=status,
            xcms_mz=mz,xcms_rt_seconds=rt,
            mass_error_ppm=ppm,rt_error_seconds=drt,
            forward_candidate_count=len(hits)))
    matched=sum(r["status"]=="reciprocal_unique" for r in rows)
    status_counts={k:sum(r["status"]==k for r in rows) for k in [
        "reciprocal_unique","no_coordinate_candidate",
        "ambiguous_multiple_candidates","not_reciprocally_unique"]}
    return rows,dict(
        schema="multiomics.sacurine-openms-independent-coordinate-comparator.v1",
        matching=dict(mz_tolerance_ppm=PPM,rt_tolerance_seconds=RT_SECONDS,
                      reciprocal_unique_required=True,nearest_neighbor_rescue=False),
        xcms_features=len(xcms),openms_consensus_features=len(openms),
        reciprocal_unique_pairs=matched,
        fraction_xcms_reciprocal=matched/len(xcms),
        fraction_openms_reciprocal=matched/len(openms),
        xcms_status_counts=status_counts,
        chemical_identity_validated=False,
        quantitative_concentration_validated=False,
        interlaboratory_validation=False,
        pipeline_superiority_claimed=False,
        interpretation="Independent preprocessing coordinate concordance only; unmatched and ambiguous features are retained."
    )

def selftest():
    x=[("X1",100,60),("X2",200,120),("X3",300,180)]
    o=[("O1",100.001,61),("O2",200.001,121),("O3",350,180)]
    rows,s=compare(x,o)
    assert s["reciprocal_unique_pairs"]==2
    assert rows[2]["status"]=="no_coordinate_candidate"
    amb=o+[("O4",100.0011,60.5)]
    rows,s=compare(x,amb)
    assert rows[0]["status"]=="ambiguous_multiple_candidates"
    xml='<consensusXML><consensusElementList count="2"><consensusElement id="1"><centroid rt="60" mz="100.1" it="2"/></consensusElement><consensusElement id="2"><centroid rt="80" mz="200.2" it="3"/></consensusElement></consensusElementList></consensusXML>'
    p=Path("/tmp/openms-comparator-selftest.consensusXML");p.write_text(xml)
    assert len(read_openms_consensus(p))==2
    p.unlink()
    print("PASS: reciprocal unique xcms/OpenMS coordinate comparator and XML parser")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--xcms",type=Path)
    ap.add_argument("--openms",type=Path)
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    if args.self_test:
        selftest(); return
    if not all([args.xcms,args.openms,args.out]):
        ap.error("--xcms --openms --out required")
    args.out.mkdir(parents=True,exist_ok=True)
    x=read_xcms(args.xcms); o=read_openms_consensus(args.openms)
    rows,summary=compare(x,o)
    with (args.out/"xcms_to_openms_coordinate_matches.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print("PASS:",summary["reciprocal_unique_pairs"],"reciprocal unique coordinate pairs;",
          "no chemical identity or pipeline superiority claim")

if __name__=="__main__":
    main()
