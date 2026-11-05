#!/usr/bin/env python3
"""Actual EATRIS-Plus v4 WBC mRNA + plasma metabolites -> plasma proteome.

Genuine source-MD5-verified, SAME exact original donor IDs in all assays.
Predicts a TRAIN-ONLY selected *plasma protein profile*, not clinical PK.
Fixed ridge+train-only PCA; same donor folds, same heldout outcomes and
available source feature/sample intersections across all arms.

No genotype, drug, experimental disease, hepatic CYP activity or patient
PK in the original public data. This does not validate Omics2PK clinical
dose adjustment, a causal blood-to-plasma mechanism, or omic superiority
over actual pharmacogenetics.
"""
from __future__ import annotations
import json
from math import isfinite
from random import Random
from statistics import mean
from tempfile import TemporaryDirectory
from pathlib import Path
import numpy as np
from scripts.audit_eatrisplus_v4_original_h5mu import (
    audit_h5mu,download_original,EATRISSourceAbstention,_index,
    SOURCE_VERSION,RECORD
)

PREDECLARED_SOURCE_ASSAYS={
    "target_proteome":"proteins",
    "white_blood_cell_mRNA":"mRNA",
    "plasma_acylcarnitines":"acylcarnitines",
    "plasma_aminoacids":"aminoacids",
}
SOURCE_ORIGIN={
    "mRNA":"source_white_blood_cells",
    "acylcarnitines":"source_plasma",
    "aminoacids":"source_plasma",
    "proteins":"source_plasma",
}
ARMS=(
    "train_mean_only",
    "WBC_mRNA",
    "plasma_targeted_metabolome",
    "WBC_mRNA_plus_plasma_targeted_metabolome",
    "WBC_mRNA_permuted_train_negative_control",
)

class EATRISBenchmarkAbstention(ValueError):
    pass


def _read_original_dense_view(h,modname,subset):
    """Only exact published sample indices, no row-order/culture joins."""
    g=h["mod"][modname]
    index=_index(g,"obs")
    if "X" not in g or not hasattr(g["X"],"shape"):
        raise EATRISBenchmarkAbstention("Source preprocessed assay matrix must be qualified dense X")
    X=g["X"]
    if tuple(X.shape)[0]!=len(index) or len(X.shape)!=2:
        raise EATRISBenchmarkAbstention("Original source assay axes inconsistent")
    # Public v4 is only 141 MB total. Individual matrices are read to
    # RAM transiently; never write per-donor output.
    arr=np.asarray(X[()],dtype=np.float64)
    if arr.shape[1]<3 or arr.shape[0]<25:
        raise EATRISBenchmarkAbstention("Source quantitative view too small")
    if np.isinf(arr).any():
        raise EATRISBenchmarkAbstention("Original public assays contain infinity")
    lookup={s:i for i,s in enumerate(index)}
    if len(lookup)!=len(index) or not set(subset).issubset(lookup):
        raise EATRISBenchmarkAbstention("Strict original donor identity join failed")
    return arr[[lookup[s] for s in subset],:]


def _preprocess_train_only(A, tr, te, *,cap, components):
    """Use source numeric scale unchanged, no assumed log/count units."""
    from sklearn.decomposition import PCA
    Xtr=A[tr].copy()
    Xte=A[te].copy()
    enough=(np.isfinite(Xtr).sum(axis=0)>=int(.80*len(tr)))
    Xtr[:,~enough]=np.nan
    Xte[:,~enough]=np.nan
    if not enough.any():
        raise EATRISBenchmarkAbstention("No source assay features have 80% training coverage")
    Xm=Xtr[:,enough];Z=Xte[:,enough]
    med=np.nanmedian(Xm,axis=0)
    Xm=np.where(np.isfinite(Xm),Xm,med)
    Z=np.where(np.isfinite(Z),Z,med)
    variance=np.var(Xm,axis=0)
    eligible=np.flatnonzero(np.isfinite(variance)&(variance>1e-14))
    if len(eligible)<3:
        raise EATRISBenchmarkAbstention("Near-constant original assay layer after train-only QC")
    # Exactly train-only feature variance selection, deterministic ties.
    selected=eligible[np.argsort(-variance[eligible],kind="stable")[:cap]]
    train=Xm[:,selected]
    test=Z[:,selected]
    mu=train.mean(axis=0)
    sd=train.std(axis=0)
    sd=np.maximum(sd,1e-8)
    train=(train-mu)/sd
    test=(test-mu)/sd
    k=min(components,train.shape[0]-1,train.shape[1])
    # full SVD deterministic within training, no test covariance leakage.
    pc=PCA(n_components=k,svd_solver="full")
    return pc.fit_transform(train),pc.transform(test)


def _target_train_only(Y,tr,te,cap=80):
    raw=Y[tr]
    good=np.isfinite(raw).sum(axis=0)>=int(.8*len(tr))
    if not good.any():
        raise EATRISBenchmarkAbstention("Original source measured plasma proteome lacks qualified numeric coverage")
    source=raw[:,good]
    meds=np.nanmedian(source,axis=0)
    source=np.where(np.isfinite(source),source,meds)
    variance=np.var(source,axis=0)
    candidates=np.flatnonzero((variance>1e-14)&np.isfinite(variance))
    if len(candidates)<10:
        raise EATRISBenchmarkAbstention("Too few original variable plasma protein targets")
    select=candidates[np.argsort(-variance[candidates],kind="stable")[:cap]]
    # The training subset alone defines the standardized target basis.
    tmean=np.mean(source[:,select],axis=0)
    tstd=np.maximum(np.std(source[:,select],axis=0),1e-8)
    Ys=(source[:,select]-tmean)/tstd
    true=Y[te][:,good][:,select]
    # Missing genuine test outcome is ignored, NOT zero imputed.
    ok=np.isfinite(true)
    if not ok.any(axis=1).all():
        raise EATRISBenchmarkAbstention("At least one heldout subject without measured protein target")
    Ytest=np.where(ok,(true-tmean)/tstd,0.0)
    return Ys,Ytest,ok,len(select)


def _report_bootstrap(errors,ref,*,seed=20261105,boots=999):
    rng=Random(seed)
    ids=list(range(len(errors)))
    obs=mean(errors[i]-ref[i] for i in ids)
    vals=[]
    for _ in range(boots):
        sample=[rng.choice(ids) for _ in ids]
        vals.append(mean(errors[i]-ref[i] for i in sample))
    vals.sort()
    lo=vals[int(.025*(boots-1))]
    hi=vals[int(.975*(boots-1))]
    return {
        "delta_standardized_protein_MAE_negative_better":round(obs,6),
        "original_person_paired_residual_bootstrap_95pct_interval":
            [round(lo,6),round(hi,6)],
        "interval_includes_zero":lo<=0<=hi,
        "bootstrap_refits_models_or_adjusts_model_choice":False,
    }


def benchmark_h5mu(path, *,strict_original=True,repeats=3):
    try:
        import h5py
        from sklearn.model_selection import KFold
        from sklearn.linear_model import Ridge
    except ImportError as exc:
        raise EATRISBenchmarkAbstention("NumPy, scikit-learn and h5py are required") from exc
    if repeats not in (1,2,3):
        raise EATRISBenchmarkAbstention("Only fixed 1–3 repetition benchmarking allowed")
    source=audit_h5mu(path,verify_original_checksum=strict_original)
    if not set(PREDECLARED_SOURCE_ASSAYS.values()).issubset(
            source["per_original_modality"]):
        raise EATRISBenchmarkAbstention("Expected original public assay names are absent")
    with h5py.File(path,"r") as h:
        nsets={s:set(_index(h["mod"][s],"obs"))
               for s in PREDECLARED_SOURCE_ASSAYS.values()}
        # Original source-index intersection, not 125-wide zero filling.
        exact=sorted(set.intersection(*nsets.values()))
        if len(exact)<30:
            raise EATRISBenchmarkAbstention("Insufficient original cross-compartment matched donors")
        arrays={s:_read_original_dense_view(h,s,exact)
                for s in PREDECLARED_SOURCE_ASSAYS.values()}
    n=len(exact)
    Y=arrays["proteins"]
    if n>180 or len(arrays)!=4:
        raise EATRISBenchmarkAbstention("Unexpected original cohort or assay count")
    measurements={arm:np.full((repeats,n),np.nan,dtype=float)
                  for arm in ARMS}
    ytarget_counts=[]
    for rep in range(repeats):
        cv=KFold(n_splits=5,shuffle=True,random_state=20261105+rep)
        for fold,(tr,te) in enumerate(cv.split(np.arange(n))):
            Ytrain,Ytest,mask,k=_target_train_only(Y,tr,te)
            ytarget_counts.append(k)
            blocks={}
            for s in ("mRNA","acylcarnitines","aminoacids"):
                blocks[s]=_preprocess_train_only(
                    arrays[s],tr,te,cap=200 if s=="mRNA" else 80,
                    components=3)
            rtr,rte=blocks["mRNA"]
            a_tr,a_te=blocks["acylcarnitines"]
            m_tr,m_te=blocks["aminoacids"]
            metab_tr=np.hstack([a_tr,m_tr])
            metab_te=np.hstack([a_te,m_te])
            randidx=np.random.default_rng(20261105+rep*13+fold).permutation(len(tr))
            xarms={
                "train_mean_only":None,
                "WBC_mRNA":(rtr,rte),
                "plasma_targeted_metabolome":(metab_tr,metab_te),
                "WBC_mRNA_plus_plasma_targeted_metabolome":
                    (np.hstack([rtr,metab_tr]),np.hstack([rte,metab_te])),
                "WBC_mRNA_permuted_train_negative_control":
                    (rtr[randidx],rte),
            }
            for name,inputs in xarms.items():
                if inputs is None: pred=np.zeros_like(Ytest)
                else:
                    Xtr,Xte=inputs
                    mdl=Ridge(alpha=10.0,fit_intercept=True)
                    mdl.fit(Xtr,Ytrain)
                    pred=mdl.predict(Xte)
                valid=np.where(mask,np.abs(Ytest-pred),0.)
                subject_mae=valid.sum(axis=1)/mask.sum(axis=1)
                measurements[name][rep,te]=subject_mae
    if any(np.isnan(mat).any() for mat in measurements.values()):
        raise EATRISBenchmarkAbstention("Some original heldout donors lack predictions")
    meanline={name:np.mean(a,axis=0).tolist()
              for name,a in measurements.items()}
    summary={
        name:round(mean(e),6) for name,e in meanline.items()
    }
    contrasts={
        name:{
            "vs_train_mean_only":_report_bootstrap(
                meanline[name],meanline["train_mean_only"]),
            **({"vs_plasma_metabolites_only":_report_bootstrap(
                meanline[name],meanline["plasma_targeted_metabolome"],
                seed=20261106)} if name=="WBC_mRNA_plus_plasma_targeted_metabolome" else {}),
        }
        for name in ARMS if name!="train_mean_only"
    }
    return {
        "status":"ORIGINAL_EATRIS_V4_H5MU_DONOR_HELDOUT_CROSS_COMPARTMENT_PROTEOME_RECONSTRUCTION_NOT_PK",
        "original_Zenodo_record":RECORD,
        "original_source_hash_verified":strict_original,
        "source_MuData_modality_count":source["n_omics_assay_modalities"],
        "original_global_person_count":source["n_original_global_observation_keys"],
        "n_exact_same_original_subjects_in_all_four_model_input_output_assays":n,
        "source_assay_to_biological_material":SOURCE_ORIGIN,
        "original_WBC_mRNA_feature_count":arrays["mRNA"].shape[1],
        "original_plasma_protein_target_feature_count":Y.shape[1],
        "original_plasma_acylcarnitine_feature_count":arrays["acylcarnitines"].shape[1],
        "original_plasma_amino_acid_feature_count":arrays["aminoacids"].shape[1],
        "fixed_outer_cross_validation":"%d_repetitions_x_5_person_disjoint_folds"%repeats,
        "train_only_preprocessing_target_and_molecular_components":True,
        "n_train_selected_target_proteins_per_fold_min":min(ytarget_counts),
        "n_train_selected_target_proteins_per_fold_max":max(ytarget_counts),
        "endpoint":"observed_plasma_protein_profile_standardized_on_training_fold",
        "measured_endpoint_is_not_clinical_drug_PK_or_PD":True,
        "cross_omic_regression_loss_is_not_genotype_or_PopPK_added_value":True,
        "heldout_standardized_plasma_protein_profile_MAE_by_arm":summary,
        "exploratory_original_biological_person_paired_contrasts":contrasts,
        "per_donor_source_rows_and_values_emitted":False,
        "original_germline_CYP3A_PGx_publicly_confirmed":False,
        "no_causal_blood_RNA_to_plasma_protein_claim":True,
    }


def live_original():
    with TemporaryDirectory(prefix="public_eatris_v4_paper_") as work:
        f=Path(work)/"private_original_mudata.h5mu"
        download_original(f)
        return benchmark_h5mu(f)


if __name__=="__main__":
    try:print(json.dumps(live_original(),indent=2,sort_keys=True))
    except (EATRISSourceAbstention,EATRISBenchmarkAbstention) as exc:
        print(json.dumps({"status":"ABSTAIN_UNQUALIFIED_OR_UNAVAILABLE_EATRIS_SOURCE_OR_PROFILES",
                          "reason_category":str(exc)[:190],
                          "actual_source_benchmark_reproduced":False}))
        raise SystemExit(2)
