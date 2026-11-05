#!/usr/bin/env python3
"""Check the immutable *actual* EATRIS-Plus v4 ISA RDF metadata.

No donor IDs or sample terms are printed: output only controlled,
method/biological-fraction TERM COUNTS as negative/positive source
provenance probes, not a verified per-person assay-material crosswalk.
Published original @ Zenodo 19136013 / eatris_isa_filled_JH.ttl.
"""
import hashlib
import json
import re
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

SOURCE="https://zenodo.org/records/19136013/files/eatris_isa_filled_JH.ttl?download=1"
EXPECTED_MD5="e724062d5f16e2eccc0cdeba2de15405"
MAX_BYTES=8_000_000
KEYWORDS={
    "explicit_plasma":r"\bplasma\b",
    "explicit_whole_blood":r"\bwhole[\s_+.-]+blood\b",
    "explicit_white_blood_cell":r"\b(?:white[\s_+.-]+blood[\s_+.-]+cells?|WBCs?|leukocytes?|buffy[\s_+.-]+coat)\b",
    "explicit_PBMC":r"\bPBMCs?\b",
    "serum":r"\bserum\b",
    "transcriptomics":r"\b(?:RNA-seq|RNAseq|transcriptom\w*|mRNA)\b",
    "proteomics":r"\b(?:proteom\w*|protein)\b",
    "metabolomics":r"\b(?:metabolom\w*|acylcarnitines?)\b",
    "lipidomics":r"\b(?:lipidom\w*|lipid)\b",
    "methylation":r"\b(?:methylat\w*|EM-seq)\b",
    "mass_spectrometry":r"\b(?:mass[\s_+.-]+spectrom\w*|LC-MS|GC-MS)\b",
}
class SourceISAMaterialAbstention(ValueError):
    pass

def audit_original_text(raw,*,strict=True):
    if not isinstance(raw,(bytes,bytearray)) or not 100<len(raw)<MAX_BYTES:
        raise SourceISAMaterialAbstention("Original ISA metadata empty or too large")
    digest=hashlib.md5(raw).hexdigest()
    if strict and digest!=EXPECTED_MD5:
        raise SourceISAMaterialAbstention("Official original Zenodo ISA RDF checksum mismatch")
    try: string=raw.decode("utf-8-sig")
    except UnicodeError as exc:
        raise SourceISAMaterialAbstention("Original ISA UTF8 metadata not readable") from exc
    vals={key:len(re.findall(pattern,string,flags=re.I)) for key,pattern in KEYWORDS.items()}
    return {
        "status":"ORIGINAL_EATRIS_V4_ISA_MATERIAL_TERM_AUDIT_NOT_PER_SUBJECT_PROVENANCE",
        "source_original_ISA":SOURCE,
        "original_ISA_MD5_verified":strict,
        "actual_original_ISA_MD5":digest if strict else "SYNTHETIC_NOT_ORIGINAL",
        "n_original_metadata_characters":len(string),
        "aggregate_material_and_platform_term_mentions":vals,
        "source_ISA_contains_specific_plasma_word":vals["explicit_plasma"]>0,
        "source_ISA_contains_specific_WBC_PBMC_wholeblood_word":any(
            vals[k]>0 for k in ("explicit_whole_blood","explicit_white_blood_cell","explicit_PBMC")),
        "literal_keyword_occurrence_not_original_assay_sample_link_proof":True,
        "individual_source_person_IDs_or_RDF_triples_never_printed":True,
        "source_provenance_only_not_PGx_patient_PK_or_causal_inference":True,
    }

def live_original():
    try:
        r=Request(SOURCE,headers={
            "Accept":"text/turtle",
            "User-Agent":"MultiOmics-public-EATRIS-v4-ISA-original-provenance/20261105"})
        with urlopen(r,timeout=45) as stream:
            raw=stream.read(MAX_BYTES+1)
    except (URLError,HTTPError,OSError,TimeoutError) as exc:
        raise SourceISAMaterialAbstention("ORIGINAL_PUBLIC_ZENODO_ISA_MATERIAL_UNAVAILABLE") from exc
    return audit_original_text(raw)

if __name__=="__main__":
    try:print(json.dumps(live_original(),sort_keys=True,indent=2))
    except SourceISAMaterialAbstention as exc:
        print(json.dumps({
            "status":"ABSTAIN_ORIGINAL_ISA_SOURCE_NOT_QUALIFIED",
            "reason_category":str(exc)[:160],
            "original_provenance_verified":False,
        }))
        raise SystemExit(2)
