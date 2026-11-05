#!/usr/bin/env python3
"""EATRIS v4 source-real outcome-free-PK capacity-matched cross-omic check.

On the SAME original n90 white-blood-cell RNA+plasma metabolome+plasma
protein subjects as the fixed-budget pilot, fit source-verified
plasma-proteome reconstruction. Restrict every omic arm to 3 total
PCs and select Ridge alpha by *inner original biological donor CV*
only. Allow one WBC and two plasma PCs, or vice versa; test a
train-permuted RNA null. Outer original persons never select model,
assay, target feature or inner hyperparameter.

This is a model sensitivity check, not a patient PK or PGx benefit
experiment. No original subject rows or source numeric matrices logged.
"""
from __future__ import annotations
import json
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scripts.audit_eatrisplus_v4_original_h5mu import (
    audit_h5mu,download_original,_index,EATRISSourceAbstention,RECORD,
)
from scripts.benchmark_eatrisplus_v4_cross_fraction import (
    _read_original_dense_view,_preprocess_train_only,_target_train_only,
    _report_bootstrap,EATRISBenchmarkAbstention,
)
from math import isfinite
from statistics import mean

ALPHAS=(1.,10.,100.,1000.)
ARMS=(
    "train_mean_only",
    "WBC_mRNA_3PC",
    "plasma_metabolome_3PC",
    "joint_WBC1_plasma2_3PC",
    "joint_WBC2_plasma1_3PC",
    "WBC_mRNA_3PC_train_permuted",
)
MODALITIES=("mRNA","proteins","acylcarnitines","aminoacids")


def _features(rna,meta,tr,te):
    rtr,rte=_preprocess_train_only(rna,tr,te,cap=200,components=3)
    mtr,mte=_preprocess_train_only(meta,tr,te,cap=80,components=3)
    if min(rtr.shape[1],mtr.shape[1])<3:
        raise EATRISBenchmarkAbstention("Fixed three omic PCs not present in train-only source")
    return rtr,rte,mtr,mte


def _arms_from_source(rtr,rte,mtr,mte,*,seed):
    rng=np.random.default_rng(seed)
    ix=rng.permutation(rtr.shape[0])
    return {
        "WBC_mRNA_3PC":(rtr,rte),
        "plasma_metabolome_3PC":(mtr,mte),
        "joint_WBC1_plasma2_3PC":
            (np.hstack([rtr[:,:1],mtr[:,:2]]),
             np.hstack([rte[:,:1],mte[:,:2]])),
        "joint_WBC2_plasma1_3PC":
            (np.hstack([rtr[:,:2],mtr[:,:1]]),
             np.hstack([rte[:,:2],mte[:,:1]])),
        "WBC_mRNA_3PC_train_permuted":(rtr[ix],rte),
    }


def _fit_predict(xt,yt,xv,alpha):
    from sklearn.linear_model import Ridge
    mdl=Ridge(alpha=alpha,fit_intercept=True)
    mdl.fit(xt,yt)
    return mdl.predict(xv)


def _per_subject_error(y,p,mask):
    return (np.where(mask,np.abs(y-p),0.).sum(axis=1) /
            mask.sum(axis=1))


def benchmark_h5mu_matched(path,*,strict_original=True,repeats=3):
    from sklearn.model_selection import KFold
    import h5py
    if repeats not in (1,2,3):
        raise EATRISBenchmarkAbstention("Original donor CV repetitions must be predeclared 1-3")
    source=audit_h5mu(path,verify_original_checksum=strict_original)
    if not set(MODALITIES).issubset(source["per_original_modality"]):
        raise EATRISBenchmarkAbstention("Exact four original source omic matrices absent")
    with h5py.File(path,"r") as h:
        original_keys=sorted(set.intersection(*(
            set(_index(h["mod"][name],"obs")) for name in MODALITIES)))
        if len(original_keys)<30:
            raise EATRISBenchmarkAbstention("Original source donor overlap too small")
        matrices={n:_read_original_dense_view(h,n,original_keys)
                  for n in MODALITIES}
    y=matrices["proteins"]
    rna=matrices["mRNA"]
    meta=np.hstack((matrices["acylcarnitines"],matrices["aminoacids"]))
    n=len(original_keys)
    out={name:np.full((repeats,n),np.nan) for name in ARMS}
    alpha_counts={name:Counter() for name in ARMS if name!="train_mean_only"}
    per_fold_target=[]
    for rep in range(repeats):
        for fold,(tr,te) in enumerate(KFold(
                n_splits=5,shuffle=True,
                random_state=20261105+rep).split(np.arange(n))):
            ty,vy,valid,nproteins=_target_train_only(y,tr,te)
            per_fold_target.append(nproteins)
            # Inner folds train-only: target variables and PCA views
            # refit from inner training donor subset. No use of
            # outer test donor info in hyperparameter or method choice.
            in_cv=KFold(n_splits=3,shuffle=True,
                        random_state=20301105+rep*10+fold)
            inner=[]
            for il,(itr,ival) in enumerate(in_cv.split(tr)):
                itr=tr[itr];ival=tr[ival]
                ytr,yval,mask,_=_target_train_only(y,itr,ival)
                rtr,rval,mtr,mval=_features(rna,meta,itr,ival)
                x=_arms_from_source(rtr,rval,mtr,mval,
                                    seed=20261105+rep*100+fold*10+il)
                inner.append((ytr,yval,mask,x))
            choices={}
            for arm in ARMS:
                if arm=="train_mean_only":continue
                scored=[]
                for alpha in ALPHAS:
                    losses=[]
                    for ytr,yval,mask,inputs in inner:
                        xt,xv=inputs[arm]
                        pred=_fit_predict(xt,ytr,xv,alpha)
                        losses.extend(_per_subject_error(yval,pred,mask))
                    scored.append((mean(losses),alpha))
                _,chosen=min(scored)
                choices[arm]=chosen
                alpha_counts[arm][str(int(chosen))]+=1
            rtr,rte,mtr,mte=_features(rna,meta,tr,te)
            inputs=_arms_from_source(
                rtr,rte,mtr,mte,seed=20261105+rep*100+fold*10+9)
            out["train_mean_only"][rep,te]=_per_subject_error(
                vy,np.zeros_like(vy),valid)
            for arm,alpha in choices.items():
                xt,xv=inputs[arm]
                pred=_fit_predict(xt,ty,xv,alpha)
                out[arm][rep,te]=_per_subject_error(vy,pred,valid)
    if any(np.isnan(v).any() for v in out.values()):
        raise EATRISBenchmarkAbstention("One original heldout donor not predicted")
    donor_means={name:np.mean(value,axis=0).tolist()
                 for name,value in out.items()}
    summary={name:round(mean(y),6) for name,y in donor_means.items()}
    ref=donor_means["train_mean_only"]
    return {
        "status":"EATRIS_V4_SOURCE_REAL_MATCHED_CAPACITY_NESTED_ALPHA_WBC_PLASMA_PROTEOME_NOT_PK",
        "original_Zenodo_record":RECORD,
        "original_h5mu_integrity_checked":strict_original,
        "n_exact_matched_original_blood_and_plasma_donors":n,
        "n_source_mRNA_features":rna.shape[1],
        "n_source_plasma_target_protein_features":y.shape[1],
        "n_source_targeted_plasma_metabolites":meta.shape[1],
        "original_individual_source_IDs_or_molecular_values_emitted":False,
        "original_drug_exposure_or_germline_PGx_present":False,
        "outer_biological_subject_split":
            "%dx5 stratification_free donor-disjoint crossvalidation"%repeats,
        "inner_biological_subject_split":"3 fold solely within each outer training sample",
        "all_prediction_arms_same_original_subjects_and_assay_outcome":True,
        "all_omic_arms_have_exactly_three_total_train_only_PCs":True,
        "outer_training_selected_plasma_protein_targets_min":min(per_fold_target),
        "outer_training_selected_plasma_protein_targets_max":max(per_fold_target),
        "inner_only_ridge_alpha_candidates":list(ALPHAS),
        "inner_selected_alpha_counts_by_arm":{
            name:dict(sorted(d.items())) for name,d in alpha_counts.items()},
        "target":"original_plasma_protein_profile_outer_train_standardized_MAE",
        "heldout_plasma_protein_target_standardized_MAE":summary,
        "paired_person_bootstrap_delta_MAE_vs_source_mean": {
            name:_report_bootstrap(vals,ref)
            for name,vals in donor_means.items() if name!="train_mean_only"},
        "matched_PC_budget_joint_minus_plasma_metab_only":{
            name:_report_bootstrap(donor_means[name],
                                   donor_means["plasma_metabolome_3PC"],
                                   seed=20261106)
            for name in (
                "joint_WBC1_plasma2_3PC",
                "joint_WBC2_plasma1_3PC",
            )},
        "no_causal_association_or_clinical_pk_efficacy_claim":True,
        "no_omic_winner_selected_on_outer_heldout_errors":True,
    }


def live_original():
    with TemporaryDirectory(prefix="eatris_real_capacity_") as d:
        p=Path(d)/"original_study_mudata.h5mu"
        download_original(p)
        return benchmark_h5mu_matched(p)


if __name__=="__main__":
    try:print(json.dumps(live_original(),sort_keys=True,indent=2))
    except (EATRISSourceAbstention,EATRISBenchmarkAbstention) as exc:
        print(json.dumps({"status":"ABSTAIN_PUBLIC_ORIGINAL_SOURCE_OR_ASSAY",
                          "source_real_benchmark_computed":False,
                          "reason_category":str(exc)[:190]}))
        raise SystemExit(2)
