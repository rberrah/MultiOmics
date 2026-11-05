#!/usr/bin/env python3
"""MultiOmics public original NUDT15 variant-omic discordance *aggregation*.

Suiter et al 2020 PNAS DOI 10.1073/pnas.1915680117;
original MaveDB meta-score-set urn:mavedb:00000055-0-1, SHA256 pinned
against original source live fetched. It contains paired VAMP-seq
protein abundance and independently screened HEK293T thiopurine
tolerance for the same variant *identity*, not original same patient.

Produces aggregate assay-label comparisons and testable biochemical
hypotheses, NO measured protein turnover, genomic DNA phasing, mRNA,
direct enzyme activity, patient clinical PD, PK, original per-variant
rows or fitted model. Clinical prediction belongs downstream.
"""
import csv
import hashlib
import io
import json
import math
import re
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

SOURCE_DOI="10.1073/pnas.1915680117"
SOURCE_URN="urn:mavedb:00000055-0-1"
SOURCE_SHA256="2e248907c0dcf54372d6ebf79c4d8a91b9ec2f2cd50d6c8aa365abec39d85a59"
URL="https://api.mavedb.org/api/v1/score-sets/"+SOURCE_URN+"/scores"
AA3={"Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
     "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"}
AMINO_PATTERN=re.compile(r"^(?:p\.)?\(?([A-Z][a-z]{2})([1-9][0-9]*)([A-Z][a-z]{2})\)?$")
SHORT_PATTERN=re.compile(r"^(?:p\.)?([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])$")


class NUDT15OriginalEvidenceError(ValueError):
    pass


def _single_aa_variant(raw):
    if not isinstance(raw,str):
        return False
    m=AMINO_PATTERN.fullmatch(raw.strip())
    if m:
        a,pos,b=m[1],int(m[2]),m[3]
        return a in AA3 and b in AA3 and a!=b and pos<=163
    m=SHORT_PATTERN.fullmatch(raw.strip())
    if m:
        return m[1]!=m[3] and int(m[2])<=163
    return False


def original_discordance_report(raw,*,check_original=True):
    if not isinstance(raw,(bytes,bytearray)) or len(raw)>5_000_000:
        raise NUDT15OriginalEvidenceError("Bounded authentic MaveDB source bytes required")
    try:
        reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    except (UnicodeDecodeError,AttributeError) as exc:
        raise NUDT15OriginalEvidenceError("Invalid original NUDT15 UTF8 source") from exc
    required={"accession","hgvs_nt","hgvs_splice","hgvs_pro",
              "score","abundance_score","sensitivity_score"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise NUDT15OriginalEvidenceError("Original NUDT15 MaveDB score-set headers incomplete")
    if len(reader.fieldnames)!=len(set(reader.fieldnames)):
        raise NUDT15OriginalEvidenceError("Original NUDT15 MaveDB duplicate headers")
    n=0;missing=0;nonstandard=0;missense=0
    stable_sensitive=0;unstable_tolerant=0
    seen=set()
    for row in reader:
        n+=1
        accession=row["accession"]
        if not isinstance(accession,str) or not accession.startswith(SOURCE_URN+"#") or accession in seen:
            raise NUDT15OriginalEvidenceError("Unqualified source variant accession")
        seen.add(accession)
        a,b=row["abundance_score"],row["sensitivity_score"]
        if a in {"","NA","NaN","nan"} or b in {"","NA","NaN","nan"}:
            missing+=1
            continue
        try:
            abundance,tolerance=float(a),float(b)
        except (ValueError,TypeError) as exc:
            raise NUDT15OriginalEvidenceError("Nonnumeric original measured assays") from exc
        if not math.isfinite(abundance) or not math.isfinite(tolerance):
            raise NUDT15OriginalEvidenceError("Nonfinite original measured values")
        if not _single_aa_variant(row["hgvs_pro"]):
            nonstandard+=1
            continue
        missense+=1
        if abundance>=.8 and tolerance<.3:
            stable_sensitive+=1
        if abundance<.3 and tolerance>=.8:
            unstable_tolerant+=1
    if check_original and (n,missing,nonstandard,missense,stable_sensitive,unstable_tolerant)!=(
            3015,170,20,2825,131,222):
        raise NUDT15OriginalEvidenceError(
            "Source original variant/assay or comparison-condition qualification differs")
    if missense<40:
        raise NUDT15OriginalEvidenceError("Not enough original paired variants")
    return {
        "status":"ORIGINAL_PUBLIC_VARIANT_LEVEL_OMIC_PD_DISCORDANCE_AGGREGATES_NOT_CLINICAL",
        "source_urn":SOURCE_URN,
        "source_doi":SOURCE_DOI,
        "number_original_variant_rows":n,
        "number_original_paired_parseable_missense":missense,
        "number_source_rows_missing_either_assay":missing,
        "number_original_nonstandard_or_unqualified_identity":nonstandard,
        "observed_variant_assay_discordances":{
            "protein_abundance_ge_0_8_with_thiopurine_tolerance_lt_0_3":stable_sensitive,
            "protein_abundance_lt_0_3_with_thiopurine_tolerance_ge_0_8":unstable_tolerant,
            "comparison_cutoffs_are_exploratory_not_clinical":True,
        },
        "qualified_original_observation_types":[
            "NUDT15_variant_specific_protein_abundance_VAMPseq",
            "same_variant_thiopurine_cellular_tolerance_screen",
        ],
        "biological_interpretation_hypotheses_not_causal_labels":[
            "PRESERVED_PROTEIN_ABUNDANCE_CAN_COEXIST_WITH_HIGH_DRUG_CELL_SENSITIVITY",
            "DECREASED_PROTEIN_ABUNDANCE_CAN_COEXIST_WITH_DRUG_CELL_TOLERANCE",
            "CHECK_ORTHOGONAL_SUBSTRATE_SPECIFIC_NUDT15_ACTIVITY_AND_DRUG_UPTAKE",
        ],
        "single_aa_variant_identity_not_same_patient_biospecimen":True,
        "no_RNA_expression_or_protein_decay_rate_measured":True,
        "no_direct_NUDT15_catalytic_activity_assay":True,
        "no_human_clinical_PK_or_PD":True,
        "original_variant_rows_predictions_and_identities_not_exported":True,
        "individual_PGx_to_PopPK_model_not_claimed":True,
        "known_effect_is_source_measured_cell_response_not_a_clinical_variant_phenotype":True,
    }


def replay_original():
    try:
        with urlopen(Request(URL,headers={"User-Agent":"MultiOmics-original-NUDT15-functional-20261105/1.0",
                                          "Accept":"text/csv"}),timeout=35) as r:
            raw=r.read(5_000_001)
    except (HTTPError,URLError,TimeoutError,OSError) as exc:
        raise NUDT15OriginalEvidenceError("OFFICIAL_ORIGINAL_MAVEDB_SOURCE_UNAVAILABLE") from exc
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA256:
        raise NUDT15OriginalEvidenceError("PINNED_ORIGINAL_MAVEDB_SCORE_SET_SHA256_MISMATCH")
    out=original_discordance_report(raw,check_original=True)
    out["original_data_fetched_live_source_SHA256"]=SOURCE_SHA256
    return out


if __name__=="__main__":
    try:
        print(json.dumps(replay_original(),indent=2,sort_keys=True))
    except NUDT15OriginalEvidenceError as exc:
        print(json.dumps({"status":"ABSTAIN_SOURCE_ORIGINAL_NOT_QUALIFIED",
                          "reason_category":str(exc)[:140]}))
        raise SystemExit(2)
