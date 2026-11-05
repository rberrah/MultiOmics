#!/usr/bin/env python3
"""Exploratory 2×2 EATRIS v4 PC capacity × Ridge shrinkage sensitivity.

Original public MD5-pinned MuData from EATRIS 125 healthy persons; 90
exactly joined WBC mRNA, plasma metabolites and measured plasma proteins.
This follows previous original-source results; it is an EXPLORATORY
follow-up, not a prospectively registered study. Orthogonally vary:
 - total molecular PCs exactly 3 versus 9, every omic arm;
 - same prespecified fixed alpha 10 versus 1_000_000.
Mean-only reference is held identical across all conditions.
Outcome: 80 fold-training-selected plasma source protein signals,
standardized using training donors, *not PK*. No site/genotype effects
or clinical dosing inference. All participant-level data remain in
memory in CI and never enter returned JSON.
"""
from __future__ import annotations
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from math import isfinite
from statistics import mean
import numpy as np

from scripts.audit_eatrisplus_v4_original_h5mu import (
    audit_h5mu,download_original,_index,RECORD,
)
from scripts.benchmark_eatrisplus_v4_cross_fraction import (
    EATRISBenchmarkAbstention,_read_original_dense_view,
    _preprocess_train_only,_target_train_only,_report_bootstrap,
)

ORIGINAL_MODS=("mRNA","proteins","acylcarnitines","aminoacids")
CAPACITIES=(3,9)
ALPHAS=(10.0,1000000.0)
ARMS=("WBC_mRNA","plasma_targeted_metabolome","WBC_plus_plasma")
SEED=20261105


def _factorial_design(source_rna,source_metabolites,tr,te):
    """Identical train-only PCA basis; slice 3/9 and 1+2/3+6 PCs."""
    rx,rv=_preprocess_train_only(source_rna,tr,te,cap=200,components=9)
    mx,mv=_preprocess_train_only(source_metabolites,tr,te,cap=80,components=9)
    if min(rx.shape[1],mx.shape[1])!=9:
        raise EATRISBenchmarkAbstention("Original assay cannot support exactly 9 train-only PCs")
    views={}
    for capacity in CAPACITIES:
        if capacity==3:ra,ma=1,2
        elif capacity==9:ra,ma=3,6
        else:raise EATRISBenchmarkAbstention("Unregistered capacity")
        views[capacity]={
            "WBC_mRNA":(rx[:,:capacity],rv[:,:capacity]),
            "plasma_targeted_metabolome":(mx[:,:capacity],mv[:,:capacity]),
            "WBC_plus_plasma":(
                np.concatenate([rx[:,:ra],mx[:,:ma]],axis=1),
                np.concatenate([rv[:,:ra],mv[:,:ma]],axis=1)),
        }
        assert all(t[0].shape[1]==capacity and t[1].shape[1]==capacity
                   for t in views[capacity].values())
    return views


def _subject_loss(true,pred,valid):
    if true.shape!=pred.shape or valid.shape!=true.shape:
        raise EATRISBenchmarkAbstention("Unqualified source target dimensions")
    n=valid.sum(axis=1)
    if (n<1).any():
        raise EATRISBenchmarkAbstention("At least one original donor has no measurable target protein")
    return np.where(valid,np.abs(true-pred),0.).sum(axis=1)/n


def _arm_key(omic,cap,alpha):
    return f"{omic}|{cap}_total_PCs|alpha_{int(alpha)}"


def benchmark_original(path,*,verify_source=True,repeats=3):
    import h5py
    from sklearn.model_selection import KFold
    from sklearn.linear_model import Ridge
    if repeats not in (1,2,3):
        raise EATRISBenchmarkAbstention("Predeclared original person repeated fold count only 1–3")
    source=audit_h5mu(path,verify_original_checksum=verify_source)
    if not set(ORIGINAL_MODS).issubset(source["per_original_modality"]):
        raise EATRISBenchmarkAbstention("Original expected white-blood-cell and plasma source assays missing")
    with h5py.File(path,"r") as f:
        ids={m:set(_index(f["mod"][m],"obs")) for m in ORIGINAL_MODS}
        shared=sorted(set.intersection(*ids.values()))
        if len(shared)!=90 and verify_source:
            raise EATRISBenchmarkAbstention("Actual source matched study participants drifted from original n90")
        if len(shared)<30:
            raise EATRISBenchmarkAbstention("Too few original same-person WBC/plasma donors")
        arrays={m:_read_original_dense_view(f,m,shared) for m in ORIGINAL_MODS}
    rna=arrays["mRNA"]
    meta=np.concatenate([arrays["acylcarnitines"],arrays["aminoacids"]],axis=1)
    target=arrays["proteins"]
    n=len(shared)
    names=[_arm_key(omic,cap,alpha) for omic in ARMS
           for cap in CAPACITIES for alpha in ALPHAS]
    output={k:np.full((repeats,n),np.nan) for k in names}
    output["mean_only"]=np.full((repeats,n),np.nan)
    target_counts=[]
    for rep in range(repeats):
        outer=KFold(n_splits=5,shuffle=True,random_state=SEED+rep)
        for tr,te in outer.split(np.arange(n)):
            train_y,test_y,valid,count=_target_train_only(target,tr,te,cap=80)
            target_counts.append(count)
            output["mean_only"][rep,te]=_subject_loss(
                test_y,np.zeros_like(test_y),valid)
            views=_factorial_design(rna,meta,tr,te)
            for cap,per_omic in views.items():
                for omic,(train_x,test_x) in per_omic.items():
                    for alpha in ALPHAS:
                        p=Ridge(alpha=alpha,fit_intercept=True).fit(
                            train_x,train_y).predict(test_x)
                        output[_arm_key(omic,cap,alpha)][rep,te]=_subject_loss(
                            test_y,p,valid)
    if any(not np.isfinite(v).all() for v in output.values()):
        raise EATRISBenchmarkAbstention("Unobserved original subject heldout model")
    individual={k:values.mean(axis=0).tolist() for k,values in output.items()}
    ref=individual["mean_only"]
    scores={k:round(mean(y),6) for k,y in individual.items()}
    delta={k:_report_bootstrap(y,ref,seed=20261105,boots=999)
           for k,y in individual.items() if k!="mean_only"}
    # Factor contrasts hold the omic arm and 3-PC/9-PC basis fixed
    # when varying shrinkage, or hold alpha fixed when changing PCs.
    contrasts={}
    for omic in ARMS:
        for cap in CAPACITIES:
            b=individual[_arm_key(omic,cap,ALPHAS[0])]
            a=individual[_arm_key(omic,cap,ALPHAS[1])]
            contrasts[f"alpha_change_{omic}_{cap}PC"]=_report_bootstrap(
                a,b,seed=20261106,boots=999)
        for alpha in ALPHAS:
            b=individual[_arm_key(omic,3,alpha)]
            a=individual[_arm_key(omic,9,alpha)]
            contrasts[f"PC_change_{omic}_alpha_{int(alpha)}"]=_report_bootstrap(
                a,b,seed=20261107,boots=999)
    return {
        "status":"SOURCE_ORIGINAL_EATRIS_V4_ORTHOGONAL_2BY2_CAPACITY_SHRINKAGE_EXPLORATORY_NOT_PK",
        "official_original_Zenodo_record":RECORD,
        "source_original_MD5_verified":verify_source,
        "original_matched_unique_persons":n,
        "original_WBC_mRNA_features":rna.shape[1],
        "original_targeted_plasma_metabolite_features":meta.shape[1],
        "original_plasma_protein_features":target.shape[1],
        "source_proteome_target_selected_outer_train_min":min(target_counts),
        "source_proteome_target_selected_outer_train_max":max(target_counts),
        "n_outer_fold_repetitions":repeats,
        "n_outer_folds_each_repetition":5,
        "independent_biological_unit":"person",
        "all_arms_exactly_same_original_donors_and_plasma_protein_outcomes":True,
        "exact_total_PC_budgets":list(CAPACITIES),
        "fixed_no_search_alpha_values":list(ALPHAS),
        "model_setting_chosen_after_previous_original_experiments":True,
        "this_is_exploratory_post_hoc_sensitivity_not_prospective_confirmation":True,
        "plasma_protein_outcome_scale":"outer_train_standardized_assay_signal",
        "heldout_protein_profile_MAE_by_factorial_arm":scores,
        "paired_original_person_bootstrap_vs_no_omics":delta,
        "orthogonal_fixed_factor_95pct_subject_bootstrap":contrasts,
        "no_repeated_CV_pseudoreplicated_as_independent_subjects":True,
        "never_retained_original_participant_IDs_or_values":True,
        "does_not_identify_relative_causal_contribution_due_to_posthoc_alphas":True,
        "no_patient_dose_AUC_CL_PGx_or_PD_claim":True,
    }

def live_original():
    with TemporaryDirectory(prefix="eatris_v4_orthogonal_") as d:
        p=Path(d)/"study.h5mu"
        download_original(p)
        return benchmark_original(p,verify_source=True)

if __name__=="__main__":
    try:print(json.dumps(live_original(),sort_keys=True,indent=2))
    except (EATRISBenchmarkAbstention,OSError,ValueError) as exc:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_SOURCE_OR_MODEL",
                          "reason_category":str(exc)[:190],
                          "real_original_benchmark_reproduced":False}))
        raise SystemExit(2)
