#!/usr/bin/env python3
"""Live original MaveDB NUDT15 engineered variant screens → MultiOmics hypotheses.

This script does NOT use patient DNA/genotypes or original clinical PK.
It downloads the immutable SHA256-pinned Suiter et al. MaveDB meta-score set,
matches *protein construct variant identity* across independently performed
VAMP-seq abundance and post-thioguanine cell-tolerance screens, and calls
the exact production sidecar interpreter. Original scores and variant IDs
are transient in memory; only aggregate source counts are printed.

The source score thresholds 0.8 and 0.3 are deliberately exploratory
reference-scale descriptive categories, NOT the source's pathogenicity,
clinical PGx star-allele or patient treatment thresholds.
"""
import csv
import gzip
import hashlib
import io
import json
import re
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

from scripts.variant_omic_evidence_sidecar import reconcile_variant_omic_evidence

URN="urn:mavedb:00000055-0-1"
SOURCE_URL=f"https://api.mavedb.org/api/v1/score-sets/{URN}/scores"
SOURCE_SHA256="2e248907c0dcf54372d6ebf79c4d8a91b9ec2f2cd50d6c8aa365abec39d85a59"
PROTEIN_SOURCE="UniProt:Q9NV35-1"
THREE={"Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
       "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"}


class NUOriginalSourceError(ValueError):
    pass


def _source_valid_amino_acid_substitution(value):
    if not isinstance(value,str):
        return False
    m=re.fullmatch(r"p\.\(?([A-Z][a-z]{2})([1-9][0-9]*)([A-Z][a-z]{2})\)?",value)
    return bool(m and m[1] in THREE and m[3] in THREE
                and m[1]!=m[3] and 1<=int(m[2])<=164)


def _raw_score(value):
    if not isinstance(value,str) or value in {"","NA","NaN","nan"}:
        return None
    try:
        score=float(value)
    except (TypeError,OverflowError,ValueError) as exc:
        raise NUOriginalSourceError("Original published score not numeric") from exc
    if score!=score or score in (float("inf"),float("-inf")):
        raise NUOriginalSourceError("Source phenotype nonfinite")
    return score


def _category(score):
    if score>=.8:
        return "reference_like"
    if score<.3:
        return "decreased"
    return "unknown"


def original_source_parse(data, *, strict=True):
    if not isinstance(data,(bytes,bytearray)) or not 100<len(data)<5_000_001:
        raise NUOriginalSourceError("Original raw MaveDB score set bytes absent or oversized")
    if strict and hashlib.sha256(data).hexdigest()!=SOURCE_SHA256:
        raise NUOriginalSourceError("Original source SHA256 does not match official published variant scores")
    try:
        reader=csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    except UnicodeDecodeError as exc:
        raise NUOriginalSourceError("Original MaveDB score matrix encoding invalid") from exc
    expected={"accession","hgvs_nt","hgvs_splice","hgvs_pro",
              "score","abundance_score","sensitivity_score"}
    if not reader.fieldnames or not expected.issubset(reader.fieldnames):
        raise NUOriginalSourceError("Original public paired MaveDB fields missing")
    if len(set(reader.fieldnames))!=len(reader.fieldnames):
        raise NUOriginalSourceError("Original duplicate CSV field")
    seen=set();valid=[];missing=0;unsupported=0;total=0
    for r in reader:
        total+=1
        accession=r.get("accession")
        if not isinstance(accession,str) or not accession.startswith(URN+"#") or accession in seen:
            raise NUOriginalSourceError("Source variant accession absent or duplicated")
        seen.add(accession)
        abundance=_raw_score(r["abundance_score"])
        tolerance=_raw_score(r["sensitivity_score"])
        if abundance is None or tolerance is None:
            missing+=1
            continue
        if not _source_valid_amino_acid_substitution(r.get("hgvs_pro")):
            unsupported+=1
            continue
        valid.append((r["hgvs_pro"],abundance,tolerance))
    if strict and (total!=3015 or len(valid)!=2844 or missing!=170 or unsupported!=1):
        raise NUOriginalSourceError("Original NUDT15 complete same-variant count changed: ABSTAIN")
    if len(valid)<60:
        raise NUOriginalSourceError("Insufficient biological variant constructs for an evidence audit")
    return valid,{"source_rows":total,"paired_supported_missense":len(valid),
                  "paired_score_missing":missing,"unparsed_or_reference_rows":unsupported}


def _observation(*,readout,omic,signal,variant,source_method,
                 original_biological_system,original_landmark):
    return {
        "assay_id":readout,
        "readout":readout,
        "omic":omic,
        "source_qualified_relative_to_reference":signal,
        # Scores have original WT and stop-variant normalizations, but
        # cutoffs are *exploratory* and are not a statistical p-value.
        "original_assay_matched_reference_validated":True,
        "assay_method":source_method,
        "original_feature_id":variant,
        "original_assay_source_url":SOURCE_URL,
        "source_variant_construct_id":variant,
        "biological_unit_id":original_biological_system,
        "tissue":"engineered_HEK293T_cells",
        "collection_landmark":original_landmark,
    }


def replay_from_bytes(b,*,strict_original=True):
    data,source=original_source_parse(b,strict=strict_original)
    total_stable_sensitive=0
    total_unstable_tolerant=0
    causal_claims=0
    selected={}
    for variant,abundance,tolerance in data:
        molecular=_category(abundance)
        pd=_category(tolerance)
        status=reconcile_variant_omic_evidence({
            "gene":"NUDT15",
            "genome_assay_status":"engineered_protein_variant_construct",
            "variant_consequence":"missense",
            "source_variant_id":variant,
            "reference_protein_accession":PROTEIN_SOURCE,
            "genetic_evidence_url":SOURCE_URL,
            "original_pairing_scope":"same_construct_variant_independent_assays",
            "observations":[
                _observation(readout="protein_total",omic="proteomics",
                             signal=molecular,variant=variant,
                             source_method="source_VAMPseq_protein_abundance",
                             original_biological_system="published_stability_screen",
                             original_landmark="untreated_protein_stability"),
                _observation(readout="drug_cell_tolerance",omic="cellular_PD",
                             signal=pd,variant=variant,
                             source_method="source_HEK293T_thioguanine_6day_PD",
                             original_biological_system="published_thioguanine_screen",
                             original_landmark="post_thioguanine_PD"),
            ],
        })
        if status["status"]!="QUALIFIED_ORIGINAL_ASSAY_DESCRIPTIVE_HYPOTHESES_NOT_CAUSAL_PROOF":
            raise NUOriginalSourceError("Source-variant interpretation refused exact pairing")
        if (status["assays_share_source_biological_unit_tissue_and_time"] or
                status["original_person_germline_genotype_measured"] or
                not status["source_variant_is_engineered_construct_not_human_genotype"] or
                status["causal_variant_mechanism_proven_by_this_report"] or
                status["protein_decay_observed_directly"] or
                status["in_vivo_drug_PK_or_PD_improvement_estimated"]):
            raise NUOriginalSourceError("DMS construct promoted to a patient molecular mechanism or clinical PK")
        hypotheses=status["hypotheses_to_test"]
        if molecular=="reference_like" and pd=="decreased":
            total_stable_sensitive+=1
            if "HYPOTHESIS_PRESERVED_PROTEIN_BUT_REDUCED_CELL_DRUG_TOLERANCE" not in hypotheses:
                raise NUOriginalSourceError("Actual source stability-versus-PD mechanism failed")
        if molecular=="decreased" and pd=="reference_like":
            total_unstable_tolerant+=1
            if "HYPOTHESIS_LOW_PROTEIN_BUT_PRESERVED_CELL_DRUG_TOLERANCE" not in hypotheses:
                raise NUOriginalSourceError("Actual source low-protein/high-PD mechanism failed")
        for h in hypotheses:
            selected[h]=selected.get(h,0)+1
        causal_claims+=int(status["causal_variant_mechanism_proven_by_this_report"])
    if strict_original and (total_stable_sensitive!=131 or total_unstable_tolerant!=222):
        raise NUOriginalSourceError("Original matched assay discordance counts changed")
    return {
        "status":"ORIGINAL_NUDT15_ENGINEERED_VARIANT_OMICS_VS_CELL_PD_HYPOTHESIS_REPLAY_NOT_CLINICAL",
        "source":SOURCE_URL,
        "source_DOI":"10.1073/pnas.1915680117",
        "original_source_SHA256":SOURCE_SHA256 if strict_original else "SYNTHETIC_NOT_ORIGINAL",
        "n_source_rows":source["source_rows"],
        "n_paired_original_missense_constructs":source["paired_supported_missense"],
        "n_missing_one_score":source["paired_score_missing"],
        "n_unparsed_or_reference_rows":source["unparsed_or_reference_rows"],
        "source_reference_qualitative_thresholds":{
            "reference_like_score_at_least":0.8,
            "decreased_score_less_than":0.3,
            "cutoff_type":"exploratory_not_clinical_or_diagnostic",
        },
        "n_preserved_protein_abundance_but_low_cell_drug_tolerance":total_stable_sensitive,
        "n_low_protein_abundance_but_preserved_cell_drug_tolerance":total_unstable_tolerant,
        "hypothesis_counts_aggregated_only":dict(sorted(selected.items())),
        "n_causal_mechanism_claims":causal_claims,
        "source_same_construct_not_same_original_patient":True,
        "no_genetic_PGx_diplotypes_from_engineered_DMS_constructs":True,
        "postdrug_cell_PD_not_eligible_as_predose_PK_feature":True,
        "no_variant_row_score_or_original_person_id_output":True,
        "no_human_PK_or_patient_toxicity_measured":True,
    }


def live_original():
    request=Request(SOURCE_URL,headers={
        "User-Agent":"MultiOmics-source-NuDT15-variant-hypotheses/2026-11-05",
        "Accept":"text/csv"})
    try:
        with urlopen(request,timeout=35) as resp:
            data=resp.read(5_000_001)
    except (HTTPError,URLError,TimeoutError,OSError) as exc:
        raise NUOriginalSourceError("PUBLIC_SOURCE_UNAVAILABLE_NO_CLINICAL_OR_SYNTHETIC_SUBSTITUTE") from exc
    return replay_from_bytes(data,strict_original=True)


if __name__=="__main__":
    try:
        print(json.dumps(live_original(),indent=2,sort_keys=True))
    except NUOriginalSourceError as exc:
        print(json.dumps({"status":"SOURCE_UNAVAILABLE_OR_CHANGED_ABSTAIN",
                          "reason_category":str(exc)[:200],
                          "original_data_reproduced":False}))
        raise SystemExit(2)
