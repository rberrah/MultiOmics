#!/usr/bin/env python3
"""Original public EATRIS-Plus 2026 v4 MuData physical-sample/source audit.

Runs on the ACTUAL publisher's original .h5mu, never on guessed arrays or
reconstructed per-person rows. Source SHA is checked against the official
Zenodo v4 published MD5. No original donor identifiers or phenotype/omics
values enter logs, GitHub, reports or artifacts. Each modality is examined
by its exact source-provided obs index, never by row order.

Public white blood cell RNA/miRNA and plasma proteome/metabolome/lipidome
are NOT separate measurements of hepatic CYP enzyme function. EATRIS does
not offer public patient drug PK or complete germline PGx.
"""
from __future__ import annotations
from collections.abc import Mapping
from hashlib import md5
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import json
import re

RECORD="https://zenodo.org/records/19136013"
SOURCE_URL=RECORD+"/files/EATRIS_PLUS_mudata.h5mu?download=1"
MD5="bbf0342797adcbcd6a3c5b242ca64573"
SOURCE_VERSION="EATRIS-PLUS_Zenodo_19136013_2026-03-20_v4"
MAX_SOURCE_BYTES=180_000_000
MAX_ORIGINAL_SUBJECTS=180

class EATRISSourceAbstention(ValueError):
    pass


def _strings(ds):
    """Return source binary/text sample labels in memory, NEVER in output."""
    if not hasattr(ds,"shape") or len(ds.shape)!=1:
        raise EATRISSourceAbstention("Original MuData donor/feature index not one-dimensional")
    vals=ds[()]
    if len(vals)>200_000:
        raise EATRISSourceAbstention("Unqualified source index too large")
    labels=[]
    for raw in vals:
        if isinstance(raw,bytes):
            try: s=raw.decode("utf-8")
            except UnicodeError as exc:
                raise EATRISSourceAbstention("Source index invalid encoding") from exc
        elif isinstance(raw,str):
            s=raw
        else:
            raise EATRISSourceAbstention("Original sample labels are not source strings")
        if not s.strip():
            raise EATRISSourceAbstention("Missing source sample/feature identity")
        labels.append(s)
    if len(set(labels))!=len(labels):
        raise EATRISSourceAbstention("Original source sample/feature identity is duplicated")
    return labels


def _index(node,which):
    if which not in node or not hasattr(node[which],"keys"):
        raise EATRISSourceAbstention("Original MuData "+which+" annotations absent")
    col=node[which]
    # AnnData writes obs/var as a dataframe whose column carrying
    # index is called _index, or by its _index attribute.
    idx=col.attrs.get("_index","_index")
    if isinstance(idx,bytes):
        idx=idx.decode("utf8")
    if not isinstance(idx,str) or idx not in col:
        if "_index" in col:idx="_index"
        elif "index" in col:idx="index"
        else:
            raise EATRISSourceAbstention(
                "No explicit original source sample/feature index; never join by row order")
    return _strings(col[idx])


def _modality_origin_guess(name):
    """Advisory only: method and specimen origin never proven by name."""
    key=re.sub("[^a-z0-9]+","",name.lower())
    if "mirna" in key or "microRNA" in key.lower():
        if "qpcr" in key or "pcr" in key:
            return "POTENTIAL_PLASMA_MIRNA_QPCR_REQUIRES_SOURCE_METHOD_MATCH"
        return "POTENTIAL_WHITE_BLOOD_CELL_SMALL_RNA_REQUIRES_SOURCE_METHOD_MATCH"
    if ("mrna" in key or "rnaseq" in key or
            key.startswith("transcript")):
        return "POTENTIAL_WHITE_BLOOD_CELL_MRNA_REQUIRES_SOURCE_METHOD_MATCH"
    if "protein" in key or "proteom" in key:
        return "POTENTIAL_PLASMA_PROTEOME_REQUIRES_SOURCE_METHOD_MATCH"
    if "metabo" in key or "acylcarn" in key:
        return "POTENTIAL_PLASMA_METABOLOME_REQUIRES_SOURCE_METHOD_MATCH"
    if "lipid" in key:
        return "POTENTIAL_PLASMA_LIPIDOME_REQUIRES_SOURCE_METHOD_MATCH"
    if "methyl" in key or key.startswith("emseq"):
        return "POTENTIAL_WHITE_BLOOD_CELL_METHYLOME_REQUIRES_SOURCE_METHOD_MATCH"
    return "UNRESOLVED_ORIGIN_DO_NOT_INFER_TISSUE_OR_OMICS_FROM_ASSAY_KEY"


def audit_h5mu(path, *, verify_original_checksum=True):
    try:
        import h5py
    except ImportError as exc:
        raise EATRISSourceAbstention("h5py required to parse original HDF5") from exc
    path=Path(path)
    if not path.is_file() or not 100 < path.stat().st_size <= MAX_SOURCE_BYTES:
        raise EATRISSourceAbstention("Original public MuData asset absent/incorrect size")
    d=md5()
    with path.open("rb") as stream:
        for part in iter(lambda:stream.read(2**20),b""):
            d.update(part)
    digest=d.hexdigest()
    if verify_original_checksum and digest!=MD5:
        raise EATRISSourceAbstention("Official Zenodo published md5 does not match original MuData bytes")
    try:
        with h5py.File(path,"r") as h:
            if "mod" not in h or not isinstance(h["mod"],h5py.Group):
                raise EATRISSourceAbstention("Source MuData /mod group absent")
            rawmods=h["mod"]
            if not 2<=len(rawmods)<=20:
                raise EATRISSourceAbstention("Expected at least two actual omic modality matrices")
            global_obs=_index(h,"obs")
            if not 5<=len(global_obs)<=MAX_ORIGINAL_SUBJECTS:
                raise EATRISSourceAbstention("Original global donor count inconsistent with source n125")
            gs=set(global_obs)
            sets={}
            info={}
            for name in sorted(rawmods):
                m=rawmods[name]
                if not isinstance(m,h5py.Group):
                    raise EATRISSourceAbstention("Actual modality not a MuData group")
                ids=_index(m,"obs")
                features=_index(m,"var")
                if not ids or not features or len(ids)>MAX_ORIGINAL_SUBJECTS:
                    raise EATRISSourceAbstention("Original modality donor/feature matrix unavailable")
                if not set(ids).issubset(gs):
                    raise EATRISSourceAbstention(
                        "Modality contains sample IDs absent from true MuData global index; no row-order matching")
                # A source-level qualified MuData view must include a
                # stored feature matrix (dense or compressed sparse).
                if "X" not in m or not isinstance(m["X"],(h5py.Dataset,h5py.Group)):
                    raise EATRISSourceAbstention("Original measured omic expression matrix missing")
                X=m["X"]
                if isinstance(X,h5py.Dataset):
                    if X.ndim!=2 or X.shape!=(len(ids),len(features)):
                        raise EATRISSourceAbstention("Original assay X axes disagree with original ID/feature keys")
                    format="HDF5_dense"
                else:
                    shape=X.attrs.get("shape")
                    if shape is None or tuple(int(v) for v in shape)!=(len(ids),len(features)):
                        raise EATRISSourceAbstention("Compressed assay axes disagreement")
                    if not {"data","indices","indptr"}.issubset(set(X.keys())):
                        raise EATRISSourceAbstention("Incomplete source sparse assay matrix")
                    format="HDF5_sparse"
                sets[name]=set(ids)
                info[name]={
                    "n_original_distinct_source_sample_keys":len(ids),
                    "n_original_measured_features":len(features),
                    "n_missing_from_global_source_index":len(gs-set(ids)),
                    "source_matrix_representation":format,
                    "potential_compartment_from_name_UNVERIFIED":_modality_origin_guess(name),
                    "assay_method_and_blood_fraction_not_proven_from_HDF5_modality_name":True,
                    "source_subject_IDs_and_numeric_measurements_never_output":True,
                }
            pairs={}
            for a,b in combinations(sorted(sets),2):
                pair=sets[a]&sets[b]
                pairs[a+" + "+b]={
                    "n_exact_original_identical_sample_keys":len(pair),
                    "n_assay_A_only":len(sets[a]-sets[b]),
                    "n_assay_B_only":len(sets[b]-sets[a]),
                    "matched_original_subject_identity_is_not_same_physical_aliquot":True,
                }
            allshared=set.intersection(*sets.values())
            return {
                "status":"ORIGINAL_ZENODO_V4_MUDATA_IDENTITY_METADATA_AUDIT_NOT_CLINICAL_PK",
                "official_record":RECORD,
                "official_record_version":SOURCE_VERSION if verify_original_checksum else "SYNTHETIC_FIXTURE",
                "original_MuData_MD5_verified":verify_original_checksum,
                "original_source_md5":digest if verify_original_checksum else "SYNTHETIC_NOT_ORIGINAL",
                "n_original_global_observation_keys":len(global_obs),
                "n_omics_assay_modalities":len(sets),
                "per_original_modality":info,
                "pairwise_exact_original_index_overlap":pairs,
                "n_exact_original_keys_shared_across_all_modalities":len(allshared),
                "public_RNA_source_material":"WHITE_BLOOD_CELLS_QUALIFY_AGAINST_ORIGINAL_ISA",
                "public_proteome_metabolome_lipidome_source":"PLASMA_QUALIFY_AGAINST_ORIGINAL_ISA",
                "original_sample_IDs_not_released":True,
                "original_individual_molecular_measurements_not_released":True,
                "original_complete_germline_PGx_unrestricted_public":False,
                "original_human_drug_PK_endpoint_present":False,
                "joint_PGx_plus_omics_superiority_can_be_tested":False,
                "mechanism_or_cross_fraction_prediction_strength_inferred":False,
            }
    except OSError as exc:
        raise EATRISSourceAbstention("Official source HDF5 cannot be parsed") from exc


def download_original(dst):
    req=Request(SOURCE_URL,headers={
        "User-Agent":"MultiOmics-EATRIS-original-source-metadata-audit/2026-11-05",
        "Accept":"application/octet-stream"})
    size=0
    try:
        with urlopen(req,timeout=90) as r, open(dst,"wb") as out:
            for part in iter(lambda:r.read(2**20),b""):
                size+=len(part)
                if size>MAX_SOURCE_BYTES:
                    raise EATRISSourceAbstention("Source asset exceeds original expected size")
                out.write(part)
    except (HTTPError,URLError,OSError,TimeoutError) as exc:
        raise EATRISSourceAbstention("ORIGINAL_PUBLIC_ZENODO_ASSET_NOT_REACHABLE") from exc
    return size


def live_original():
    with TemporaryDirectory(prefix="original_eatris_v4_") as folder:
        f=Path(folder)/"original_private_processing.h5mu"
        download_original(f)
        return audit_h5mu(f)


if __name__=="__main__":
    try:print(json.dumps(live_original(),indent=2,sort_keys=True))
    except EATRISSourceAbstention as exc:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_EATRIS_SOURCE_UNQUALIFIED",
                          "original_asset_verified":False,
                          "reason_category":str(exc)[:180]}))
        raise SystemExit(2)
