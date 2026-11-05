#!/usr/bin/env python3
"""Original CYP2C9 known missense constructs: independent VAMP-seq and Click-seq.

Original Amorosi et al., AJHG 2021, DOI 10.1016/j.ajhg.2021.07.001.
Actual assayed *engineered protein constructs*, not patients or germline PGx.
This script source-authenticates the original GitHub git blob by SHA1 and
reuses the MultiOmics biological evidence reporter for every original
paired missense variant. It does not train a PK model and never exports
variant, original molecular score or sample-level rows.
"""
from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from scripts.variant_omic_evidence_sidecar import reconcile_variant_omic_evidence

SHA="1a3041077f399dba6b4115877d31bd9def115787"
URL="https://api.github.com/repos/dunhamlab/CYP2C9/git/blobs/"+SHA
DOI="10.1016/j.ajhg.2021.07.001"
ORIGINAL_VARIANT_ROWS=8701
PAIRED_MISSENSE=4421
TARGETED_INCONSISTENCY_COUNT=400
HEADERS={"variant","class","position","start","end","activity_class",
         "activity_score","abundance_class","abundance_score"}
AMINO="ACDEFGHIKLMNPQRSTVWY"


class OriginalCYP2C9Unavailable(ValueError):
    pass


def _source_label(value):
    return {"wt-like":"reference_like",
            "decreased":"decreased",
            "nonsense-like":"decreased"}.get(value,"unknown")


def _parse_original(content,*,strict=True):
    if not isinstance(content,(bytes,bytearray)) or not 100<len(content)<5_000_001:
        raise OriginalCYP2C9Unavailable("Original public CYP2C9 source CSV bytes missing")
    actual=hashlib.sha1(b"blob "+str(len(content)).encode()+b"\0"+content).hexdigest()
    if strict and actual!=SHA:
        raise OriginalCYP2C9Unavailable("Authenticated original Git blob SHA1 mismatch")
    try:reader=csv.DictReader(io.StringIO(content.decode("utf-8-sig")))
    except UnicodeDecodeError as exc:
        raise OriginalCYP2C9Unavailable("Source CSV encoding invalid") from exc
    if not reader.fieldnames or not HEADERS.issubset(reader.fieldnames):
        raise OriginalCYP2C9Unavailable("Original paired readout headers absent")
    if len(set(reader.fieldnames))!=len(reader.fieldnames):
        raise OriginalCYP2C9Unavailable("Duplicate original assay header")
    seen=set();out=[];all_rows=0
    for rec in reader:
        all_rows+=1
        v=rec.get("variant")
        if not isinstance(v,str) or not v or v in seen:
            raise OriginalCYP2C9Unavailable("Duplicated/invalid original variant construct")
        seen.add(v)
        if rec["class"]!="missense" or any(rec[k]=="NA" for k in (
                "activity_score","abundance_score","activity_class","abundance_class")):
            continue
        try:
            position=int(rec["position"])
            activity=float(rec["activity_score"])
            abundance=float(rec["abundance_score"])
        except (ValueError,TypeError) as exc:
            raise OriginalCYP2C9Unavailable("Original measured paired score is not numeric") from exc
        if (not 1<=position<=500 or rec["start"] not in AMINO
                or rec["end"] not in AMINO or rec["start"]==rec["end"]
                or v!=f'{rec["start"]}{position}{rec["end"]}'
                or any(not(-1e6<x<1e6) for x in (activity,abundance))):
            raise OriginalCYP2C9Unavailable("Original protein variant and numeric readout incongruent")
        out.append((v,rec["abundance_class"],rec["activity_class"]))
    if strict and (all_rows!=ORIGINAL_VARIANT_ROWS or len(out)!=PAIRED_MISSENSE):
        raise OriginalCYP2C9Unavailable("Original known variant source count or paired set shifted")
    if len(out)<30:
        raise OriginalCYP2C9Unavailable("Insufficient original paired protein variants")
    return out,all_rows


def _observation(variant,*,kind,omic,source_call,assay,host,landmark):
    return {
        "assay_id":kind,
        "readout":kind,
        "omic":omic,
        "source_qualified_relative_to_reference":_source_label(source_call),
        "original_assay_matched_reference_validated":True,
        "assay_method":assay,
        "original_feature_id":variant,
        "original_assay_source_url":URL,
        "source_variant_construct_id":variant,
        "biological_unit_id":host,
        "tissue":host,
        "collection_landmark":landmark,
    }


def replay_from_csv(data,*,strict=True):
    variants,n_total=_parse_original(data,strict=strict)
    n_contradiction=0
    aggregate={}
    causal=0
    for vid,prot,activity in variants:
        r=reconcile_variant_omic_evidence({
            "gene":"CYP2C9",
            "genome_assay_status":"engineered_protein_variant_construct",
            "source_variant_id":vid,
            "reference_protein_accession":"UniProt:P11712",
            "variant_consequence":"missense",
            "genetic_evidence_url":URL,
            "original_pairing_scope":"same_construct_variant_independent_assays",
            "observations":[
                _observation(vid,kind="protein_total",omic="proteomics",
                             source_call=prot,assay="VAMPseq_variant_protein_abundance",
                             host="source_cultured_human_cells",landmark="source_untreated_abundance_screen"),
                _observation(vid,kind="measured_substrate_activity",omic="functional_activity",
                             source_call=activity,assay="Clickseq_activity_based_probe",
                             host="source_engineered_yeast",landmark="source_independent_activity_screen"),
            ],
        })
        if (r["source_pairing_scope"]!="same_construct_variant_independent_assays"
                or r["assays_share_source_biological_unit_tissue_and_time"]
                or r["original_person_germline_genotype_measured"]
                or not r["source_variant_is_engineered_construct_not_human_genotype"]
                or r["protein_decay_observed_directly"]
                or r["in_vivo_drug_PK_or_PD_improvement_estimated"]
                or r["causal_variant_mechanism_proven_by_this_report"]):
            raise OriginalCYP2C9Unavailable("Invalid biological evidence promoted to patient functional PK inference")
        for name in r["hypotheses_to_test"]:
            aggregate[name]=aggregate.get(name,0)+1
        if activity=="decreased" and prot=="wt-like":
            n_contradiction+=1
            if "HYPOTHESIS_PRESERVED_ABUNDANCE_LOST_MISSENSE_FUNCTION" not in r["hypotheses_to_test"]:
                raise OriginalCYP2C9Unavailable("Original source direct CYP protein-function discordance not detected")
        causal+=int(r["causal_variant_mechanism_proven_by_this_report"])
    if strict and n_contradiction!=TARGETED_INCONSISTENCY_COUNT:
        raise OriginalCYP2C9Unavailable("Source originally reported sharp phenotype contradiction count shifted")
    return {
        "status":"ORIGINAL_CYP2C9_VARIANT_VAMPSEQ_CLICKSEQ_MOLECULAR_REPLAY_NOT_CLINICAL_PK",
        "original_source":URL,
        "source_DOI":DOI,
        "original_git_blob_sha":SHA if strict else "SYNTHETIC",
        "n_original_source_variants_all_classes":n_total,
        "n_original_source_paired_missense_constructs":len(variants),
        "n_original_decreased_activity_but_wtlike_abundance":n_contradiction,
        "original_evidence_hypothesis_class_counts":dict(sorted(aggregate.items())),
        "n_automatically_proven_functional_or_decay_causal_mechanisms":causal,
        "yeast_activity_and_human_cultured_cell_protein_not_same_original_biospecimen":True,
        "original_substrate_activity_probe_not_specific_clinical_drug_clearance":True,
        "no_original_patient_germline_PGx_or_drug_PK":True,
        "all_original_individual_variant_rows_and_scores_stay_in_memory":True,
    }


def live_original():
    try:
        request=Request(URL,headers={"User-Agent":"MultiOmics-original-public-CYP2C9-source/20261105",
                                     "Accept":"application/vnd.github+json"})
        with urlopen(request,timeout=40) as response:
            blob=json.loads(response.read(6_000_000))
    except (HTTPError,URLError,OSError,ValueError,TimeoutError) as exc:
        raise OriginalCYP2C9Unavailable("Public source GitHub API unavailable") from exc
    if blob.get("sha")!=SHA or blob.get("encoding")!="base64":
        raise OriginalCYP2C9Unavailable("Original source remote SHA / encoding mismatch")
    try:content=base64.b64decode(blob["content"])
    except (ValueError,KeyError,TypeError) as exc:
        raise OriginalCYP2C9Unavailable("Source Git blob base64 invalid") from exc
    return replay_from_csv(content,strict=True)


if __name__=="__main__":
    try:print(json.dumps(live_original(),indent=2,sort_keys=True))
    except OriginalCYP2C9Unavailable as exc:
        print(json.dumps({"status":"SOURCE_ORIGINAL_CYP2C9_UNQUALIFIED_ABSTAIN",
                          "source_replay_validated":False,
                          "reason_category":str(exc)[:160]}))
        raise SystemExit(2)
