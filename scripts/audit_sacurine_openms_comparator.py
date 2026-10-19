#!/usr/bin/env python3
"""Independent OpenMS-vs-xcms audit on authentic Sacurine v3 LC-MS features.

No overlap threshold defines success. A successful run establishes only that
the prespecified comparison executed with exact source identity and complete
numerical accounting.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from bisect import bisect_left,bisect_right
from pathlib import Path

def finite(value):
    try:
        x=float(value)
    except (TypeError,ValueError):
        return None
    return x if math.isfinite(x) else None

def read_source(path:Path):
    with path.open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    if len(rows)!=12:
        raise ValueError("Sacurine source metadata must have exactly 12 rows")
    role_map={"sample":"biological","QC":"pooled_qc","blk":"blank"}
    out={}
    orders=set()
    for row in rows:
        name=(row.get("sample_name") or "").strip()
        role=role_map.get((row.get("class") or "").strip())
        try: order=int(row.get("injectionOrder") or "")
        except ValueError as e: raise ValueError("Invalid source injection order") from e
        if not name or role is None or name in out or order in orders:
            raise ValueError("Invalid/duplicate authentic source metadata")
        out[name]={"role":role,"order":order}
        orders.add(order)
    if sorted(x["role"] for x in out.values()).count("biological")!=6:
        raise ValueError("Expected six biological samples")
    if sum(x["role"]=="pooled_qc" for x in out.values())!=3 or sum(x["role"]=="blank" for x in out.values())!=3:
        raise ValueError("Expected three pooled QC and three blanks")
    return out

def parse_openms_text(path:Path):
    map_header=None;cons_header=None;maps={};cons=[]
    with path.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.reader(f,delimiter="\t")
        for raw in reader:
            if not raw: continue
            kind=raw[0].lstrip("#").strip()
            if raw[0].startswith("#") and kind=="MAP":
                map_header=raw[1:];continue
            if raw[0].startswith("#") and kind=="CONSENSUS":
                cons_header=raw[1:];continue
            if kind=="MAP" and not raw[0].startswith("#"):
                if not map_header or len(raw)-1<len(map_header):
                    raise ValueError("Malformed OpenMS MAP row")
                row=dict(zip(map_header,raw[1:]))
                mid=str(row.get("id",""))
                if not mid or mid in maps: raise ValueError("Duplicate/missing OpenMS map id")
                maps[mid]=row
            elif kind=="CONSENSUS" and not raw[0].startswith("#"):
                if not cons_header or len(raw)-1<len(cons_header):
                    raise ValueError("Malformed OpenMS CONSENSUS row")
                cons.append(dict(zip(cons_header,raw[1:])))
    if not maps or not cons:
        raise ValueError("OpenMS TextExporter output lacks MAP or CONSENSUS records")
    return maps,cons

def sample_from_map_filename(filename:str)->str:
    name=Path(filename).name
    for suffix in (".featureXML",".featurexml"):
        if name.endswith(suffix): return name[:-len(suffix)]
    return Path(name).stem

def openms_table(text_path:Path,source:dict):
    maps,rows=parse_openms_text(text_path)
    map_to_sample={}
    for mid,row in maps.items():
        sample=sample_from_map_filename(row.get("filename",""))
        if sample not in source:
            raise ValueError("OpenMS map cannot be reconciled to authentic sample: "+sample)
        if sample in map_to_sample.values(): raise ValueError("Duplicate OpenMS sample map")
        map_to_sample[mid]=sample
    if set(map_to_sample.values())!=set(source):
        raise ValueError("OpenMS consensus does not contain exactly the 12 authentic runs")
    ordered=sorted(source,key=lambda s:source[s]["order"])
    features=[]
    for i,row in enumerate(rows,1):
        mz=finite(row.get("mz_cf"));rt=finite(row.get("rt_cf"))
        if mz is None or mz<=0 or rt is None or rt<0:
            raise ValueError("Invalid OpenMS consensus m/z or RT")
        values={}
        for mid,sample in map_to_sample.items():
            val=finite(row.get("intensity_"+mid))
            if val is not None and val<0: raise ValueError("Negative OpenMS intensity")
            values[sample]=val
        # Blank-only signals do not define the comparison feature universe.
        nonblank=[values[s] for s in ordered if source[s]["role"]!="blank"]
        if not any(x is not None and x>0 for x in nonblank):
            continue
        features.append({"feature_id":f"OPENMS_CF{i:06d}","mzmed":mz,"rtmed":rt,
                         "values":values})
    if not features: raise ValueError("No nonblank-defined OpenMS consensus features")
    return ordered,features

def qc_rows(features,ordered,source):
    qcs=[s for s in ordered if source[s]["role"]=="pooled_qc"]
    bios=[s for s in ordered if source[s]["role"]=="biological"]
    blanks=[s for s in ordered if source[s]["role"]=="blank"]
    out=[]
    for f in features:
        q=[f["values"][s] for s in qcs if f["values"][s] is not None]
        b=[f["values"][s] for s in blanks if f["values"][s] is not None]
        x=[f["values"][s] for s in bios if f["values"][s] is not None]
        rsd=None
        if len(q)==3 and statistics.mean(q)>0:
            rsd=statistics.stdev(q)/statistics.mean(q)
        blank=None
        if len(b)>=2 and len(x)>=2 and statistics.median(x)>0:
            blank=statistics.median(b)/statistics.median(x)
        out.append({"feature_id":f["feature_id"],"pooled_qc_rsd":rsd,
                    "blank_over_biological_ratio":blank})
    return out

def read_csv_dict(path:Path):
    with path.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))

def coordinates(rows):
    out=[]
    seen=set()
    for row in rows:
        fid=str(row.get("feature_id",""))
        mz=finite(row.get("mzmed"));rt=finite(row.get("rtmed"))
        if not fid or fid in seen or mz is None or mz<=0 or rt is None or rt<0:
            raise ValueError("Invalid xcms/OpenMS feature coordinate table")
        seen.add(fid);out.append((fid,mz,rt))
    if not out: raise ValueError("Empty coordinate table")
    return out

def candidate_lists(source,target,ppm=20.0,rt_seconds=30.0):
    sorted_target=sorted((mz,idx,rt) for idx,(_,mz,rt) in enumerate(target))
    masses=[x[0] for x in sorted_target]
    answer=[]
    for _,mz,rt in source:
        delta=mz*ppm/1e6
        lo=bisect_left(masses,mz-delta);hi=bisect_right(masses,mz+delta)
        hits=[]
        for tmz,tidx,trt in sorted_target[lo:hi]:
            if abs(tmz-mz)/mz*1e6<=ppm and abs(trt-rt)<=rt_seconds:hits.append(tidx)
        answer.append(hits)
    return answer

def reciprocal_matches(ref,alt,ppm=20.0,rt_seconds=30.0):
    fwd=candidate_lists(ref,alt,ppm,rt_seconds)
    rev=candidate_lists(alt,ref,ppm,rt_seconds)
    rows=[];used=set()
    for i,(fid,mz,rt) in enumerate(ref):
        hits=fwd[i];status="no_coordinate_candidate";aid=None;ppm_err=None;rt_err=None
        if len(hits)>1: status="ambiguous_multiple_candidates"
        elif len(hits)==1:
            j=hits[0]
            if len(rev[j])==1 and rev[j][0]==i:
                status="reciprocal_unique"
                aid=alt[j][0];used.add(j)
                ppm_err=(alt[j][1]-mz)/mz*1e6;rt_err=alt[j][2]-rt
            else: status="not_reciprocally_unique"
        rows.append({"reference_feature":fid,"alternative_feature":aid or "",
                     "candidate_count":len(hits),"match_status":status,
                     "matched_mass_error_ppm":ppm_err,
                     "matched_rt_error_seconds":rt_err})
    return rows,used

def write_csv(path,fieldnames,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames,extrasaction="ignore")
        w.writeheader()
        for row in rows:w.writerow({k:("" if row.get(k) is None else row.get(k)) for k in fieldnames})

def paired_agreement(matches,ref_qc,alt_qc):
    ref={r["feature_id"]:r for r in ref_qc};alt={r["feature_id"]:r for r in alt_qc}
    paired=[m for m in matches if m["match_status"]=="reciprocal_unique"]
    def metric(key,cutoff,inclusive=False):
        values=[]
        for m in paired:
            a=finite(ref.get(m["reference_feature"],{}).get(key))
            b=finite(alt.get(m["alternative_feature"],{}).get(key))
            if a is None or b is None:continue
            fa=a>=cutoff if inclusive else a>cutoff
            fb=b>=cutoff if inclusive else b>cutoff
            values.append(fa==fb)
        return {"n_evaluable_both":len(values),
                "agreement_fraction":sum(values)/len(values) if values else None}
    return {"qc_rsd_gt_30pct":metric("pooled_qc_rsd",.30,False),
            "blank_over_biological_ge_20pct":metric("blank_over_biological_ratio",.20,True)}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--metadata",required=True,type=Path)
    ap.add_argument("--openms-text",required=True,type=Path)
    ap.add_argument("--xcms-coordinates",required=True,type=Path)
    ap.add_argument("--xcms-qc",required=True,type=Path)
    ap.add_argument("--openms-image-digest",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    source=read_source(args.metadata)
    ordered,features=openms_table(args.openms_text,source)
    oqc=qc_rows(features,ordered,source)
    open_coords=[{"feature_id":f["feature_id"],"mzmed":f["mzmed"],"rtmed":f["rtmed"]} for f in features]
    ref=coordinates(read_csv_dict(args.xcms_coordinates));alt=coordinates(open_coords)
    matches,used=reciprocal_matches(ref,alt)
    xqc=read_csv_dict(args.xcms_qc)
    agreement=paired_agreement(matches,xqc,oqc)
    args.output.mkdir(parents=True,exist_ok=True)
    write_csv(args.output/"openms_feature_coordinates.csv",["feature_id","mzmed","rtmed"],open_coords)
    matrix=[]
    for f in features:
        row={"feature_id":f["feature_id"]}
        row.update({s:f["values"][s] for s in ordered})
        matrix.append(row)
    write_csv(args.output/"openms_feature_matrix.csv",["feature_id",*ordered],matrix)
    write_csv(args.output/"openms_qc_blank_screens.csv",
              ["feature_id","pooled_qc_rsd","blank_over_biological_ratio"],oqc)
    write_csv(args.output/"xcms_openms_reciprocal_matches.csv",
              ["reference_feature","alternative_feature","candidate_count","match_status",
               "matched_mass_error_ppm","matched_rt_error_seconds"],matches)
    matched=sum(m["match_status"]=="reciprocal_unique" for m in matches)
    status_counts={k:sum(m["match_status"]==k for m in matches) for k in
      ("reciprocal_unique","no_coordinate_candidate","ambiguous_multiple_candidates",
       "not_reciprocally_unique")}
    evaluable_rsd=sum(finite(r["pooled_qc_rsd"]) is not None for r in oqc)
    high_rsd=sum((finite(r["pooled_qc_rsd"]) or -1)>.30 for r in oqc)
    evaluable_blank=sum(finite(r["blank_over_biological_ratio"]) is not None for r in oqc)
    high_blank=sum((finite(r["blank_over_biological_ratio"]) or -1)>=.20 for r in oqc)
    report={
      "schema":"multiomics.sacurine-openms-independent-comparator.v1",
      "source":"https://zenodo.org/records/10130758",
      "openms_version":"3.6.0",
      "openms_image_digest":args.openms_image_digest.read_text(encoding="utf-8").strip(),
      "comparison_policy":{"ppm":20,"rt_seconds":30,"matching":"reciprocal_unique_only",
        "minimum_success_threshold":None,"posthoc_parameter_tuning":False},
      "authenticated_runs":12,"roles":{"biological":6,"pooled_qc":3,"blank":3},
      "xcms_reference_features":len(ref),"openms_nonblank_defined_features":len(alt),
      "reciprocal_unique_matches":matched,
      "fraction_xcms_features_matched":matched/len(ref),
      "fraction_openms_features_matched":matched/len(alt),
      "match_status_counts":status_counts,
      "openms_qc":{"rsd_evaluable":evaluable_rsd,"rsd_gt_30pct":high_rsd,
        "blank_ratio_evaluable":evaluable_blank,"blank_ratio_ge_20pct":high_blank},
      "paired_screen_agreement":agreement,
      "chemical_identity_validated":False,"drift_correction_validated":False,
      "carryover_validated":False,"interlaboratory_accuracy_validated":False,
      "biological_validation":False,"omics_pk_linkage_validated":False,
      "scientific_certification":False,
      "interpretation":"Independent software preprocessing comparison only. Coordinate overlap and QC-screen agreement are descriptive; unmatched/ambiguous features are retained as negative evidence."
    }
    (args.output/"openms_xcms_comparison_summary.json").write_text(
      json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("PASS: independent OpenMS comparison executed; reciprocal matches",
          matched,"/",len(ref),"xcms and",len(alt),"OpenMS; NO chemical/biological/PK certification")

if __name__=="__main__":
    main()
