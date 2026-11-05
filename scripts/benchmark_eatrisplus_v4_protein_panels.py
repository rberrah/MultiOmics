#!/usr/bin/env python3
"""Public EATRIS v4 source assay sensitivity: proteins vs proteomics-ipp.

Exploratory post-hoc comparison, not PK. Both protein target panels are
measured on exactly the SAME original donors, using identical folds, source
input features, fixed PC budgets and Ridge penalties. The protein panels
are different measured outcomes: raw and standardized MAEs are NOT
directly comparable as interchangeable targets or analytical platforms.
No subject IDs, sample rows, protein identifiers, or per-subject errors
are emitted. Source original MD5 check fails closed.
"""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from scripts.audit_eatrisplus_v4_original_h5mu import (
    RECORD, EATRISSourceAbstention, _index, audit_h5mu, download_original,
)
from scripts.benchmark_eatrisplus_v4_cross_fraction import (
    EATRISBenchmarkAbstention, _read_original_dense_view, _preprocess_train_only,
    _target_train_only, _report_bootstrap,
)
from scripts.benchmark_eatrisplus_v4_orthogonal_capacity_shrinkage import (
    _factorial_design, _subject_loss,
)

SOURCE_MODS = ("mRNA", "acylcarnitines", "aminoacids", "proteins", "proteomics-ipp")
TARGET_PANELS = ("proteins", "proteomics-ipp")
ARMS = ("WBC_mRNA", "plasma_targeted_metabolome", "WBC_plus_plasma")
CAPACITIES = (3, 9)
ALPHAS = (10.0, 1_000_000.0)
SEED = 20261105


def _key(arm, capacity, alpha):
    return f"{arm}|{capacity}PC|alpha_{int(alpha)}"


def _scale_contrast(errors, reference):
    """Person-level relative skill vs the same panel's mean-only baseline."""
    ref = float(np.mean(reference))
    if not np.isfinite(ref) or ref <= 0:
        raise EATRISBenchmarkAbstention("Unqualified target-specific baseline")
    return round(float((ref - np.mean(errors)) / ref), 8)


def benchmark_proteome_panels(path, *, verify_source=True, repeats=3):
    import h5py
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold

    if repeats not in (1, 2, 3):
        raise EATRISBenchmarkAbstention("Only 1 to 3 fixed donor-disjoint repetitions")
    provenance = audit_h5mu(path, verify_original_checksum=verify_source)
    if not set(SOURCE_MODS).issubset(provenance["per_original_modality"]):
        raise EATRISBenchmarkAbstention("Required original proteomic panels not present")

    with h5py.File(path, "r") as h:
        sample_sets = {m: set(_index(h["mod"][m], "obs")) for m in SOURCE_MODS}
        matched = sorted(set.intersection(*sample_sets.values()))
        # Public n90 mRNA is the limiting assay. Never invent missing protein rows.
        if len(matched) < 30 or len(matched) > 90:
            raise EATRISBenchmarkAbstention("Unqualified five-assay same-person overlap")
        arrays = {m: _read_original_dense_view(h, m, matched) for m in SOURCE_MODS}
        # Exact source feature names are diagnostic only: do not assume
        # equivalent protein identity, isoform, physical assay or aliquot.
        original_feature_keys = {p: _index(h["mod"][p], "var")
                                 for p in TARGET_PANELS}
        exact_key_overlap = len(set(original_feature_keys["proteins"])
                                & set(original_feature_keys["proteomics-ipp"]))
        p_normalized = [k.casefold().strip() for k in original_feature_keys["proteins"]]
        i_normalized = [k.casefold().strip() for k in original_feature_keys["proteomics-ipp"]]
        norm_collision = (len(set(p_normalized)) != len(p_normalized)
                          or len(set(i_normalized)) != len(i_normalized))
        casefold_overlap = (None if norm_collision else
                            len(set(p_normalized) & set(i_normalized)))

    rna = arrays["mRNA"]
    metabolites = np.hstack([arrays["acylcarnitines"], arrays["aminoacids"]])
    n = len(matched)
    keys = ["mean_only"] + [_key(a, k, alpha) for a in ARMS
                             for k in CAPACITIES for alpha in ALPHAS]
    keys.append("WBC_mRNA_train_permuted|3PC|alpha_10")
    errors = {p: {name: np.full((repeats, n), np.nan) for name in keys}
              for p in TARGET_PANELS}
    target_selection_counts = {p: [] for p in TARGET_PANELS}
    available_targets = {p: [] for p in TARGET_PANELS}

    for rep in range(repeats):
        cv = KFold(n_splits=5, shuffle=True, random_state=SEED + rep)
        for fold, (tr, te) in enumerate(cv.split(np.arange(n))):
            # Predictor preprocessing is identical for both output matrices.
            views = _factorial_design(rna, metabolites, tr, te)
            perm = np.random.default_rng(SEED + rep * 13 + fold).permutation(len(tr))
            for panel in TARGET_PANELS:
                train_y, test_y, observed, n_targets = _target_train_only(
                    arrays[panel], tr, te, cap=80)
                target_selection_counts[panel].append(n_targets)
                available_targets[panel].append(int(observed.sum()))
                errors[panel]["mean_only"][rep, te] = _subject_loss(
                    test_y, np.zeros_like(test_y), observed)
                for capacity, arm_views in views.items():
                    for arm, (train_x, test_x) in arm_views.items():
                        for alpha in ALPHAS:
                            predicted = Ridge(alpha=alpha, fit_intercept=True).fit(
                                train_x, train_y).predict(test_x)
                            errors[panel][_key(arm, capacity, alpha)][rep, te] = (
                                _subject_loss(test_y, predicted, observed))
                xtr, xte = views[3]["WBC_mRNA"]
                permuted = Ridge(alpha=10., fit_intercept=True).fit(
                    xtr[perm], train_y).predict(xte)
                errors[panel]["WBC_mRNA_train_permuted|3PC|alpha_10"][rep, te] = (
                    _subject_loss(test_y, permuted, observed))

    if any(not np.isfinite(v).all() for panel in errors.values() for v in panel.values()):
        raise EATRISBenchmarkAbstention("Original assay heldout predictions incomplete")
    # Repeated predictions of one biological person are averaged FIRST; n is
    # the number of unique matched persons, not 3*n or 15*n.
    person_errors = {p: {key: val.mean(axis=0) for key, val in errors[p].items()}
                     for p in TARGET_PANELS}
    summaries = {}
    for panel in TARGET_PANELS:
        reference = person_errors[panel]["mean_only"]
        models = person_errors[panel]
        summaries[panel] = {
            "n_measured_original_features": int(arrays[panel].shape[1]),
            "original_training_selected_targets_per_fold_min":
                min(target_selection_counts[panel]),
            "original_training_selected_targets_per_fold_max":
                max(target_selection_counts[panel]),
            "heldout_observed_target_cells_per_fold_min": min(available_targets[panel]),
            "heldout_observed_target_cells_per_fold_max": max(available_targets[panel]),
            "mean_only_standardized_target_MAE": round(float(np.mean(reference)), 6),
            "standardized_target_MAE_by_omic_arm": {
                name: round(float(np.mean(values)), 6)
                for name, values in models.items() if name != "mean_only"
            },
            "within_target_panel_relative_skill_vs_mean_only": {
                name: _scale_contrast(values, reference)
                for name, values in models.items() if name != "mean_only"
            },
            "exploratory_paired_person_delta_MAE_vs_same_panel_mean_only": {
                name: _report_bootstrap(
                    values.tolist(), reference.tolist(), seed=SEED, boots=999)
                for name, values in models.items() if name != "mean_only"
            },
        }

    # Across panels, compare within-panel RELATIVE changes from each
    # panel's own benchmark. This remains non-causal: different analytes,
    # missingness, signal distributions, and measurement technologies.
    relative_sensitivity = {}
    for name in keys:
        if name == "mean_only":
            continue
        pair = {}
        for panel in TARGET_PANELS:
            baseline = person_errors[panel]["mean_only"]
            score = person_errors[panel][name]
            pair[panel] = (baseline - score) / float(np.mean(baseline))
        # Same-donor resampling of the paired *relative gain* difference.
        # Difference is IPP minus proteins. Positive = higher relative gain
        # for IPP, NOT proof that the platform itself causes this difference.
        rel_delta = (pair["proteomics-ipp"] - pair["proteins"]).tolist()
        relative_sensitivity[name] = _report_bootstrap(
            rel_delta, [0.] * n, seed=SEED + 1, boots=999)

    return {
        "status": "SOURCE_ORIGINAL_EATRIS_V4_TWO_PROTEIN_PANELS_MATCHED_PERSONS_EXPLORATORY_NOT_PK",
        "official_original_Zenodo": RECORD,
        "original_h5mu_MD5_verified": verify_source,
        "n_original_independent_matched_people_in_all_five_assays": n,
        "original_assay_features": {
            "white_blood_cell_mRNA": int(rna.shape[1]),
            "plasma_acylcarnitines_plus_amino_acids": int(metabolites.shape[1]),
        },
        "target_panels": list(TARGET_PANELS),
        "original_assay_source_fraction": {
            "mRNA": "white_blood_cells_from_study_description",
            "acylcarnitines": "plasma_from_study_description",
            "aminoacids": "plasma_from_study_description",
            "proteins": "plasma_related_method_not_fully_qualified",
            "proteomics-ipp": "plasma_related_distinct_method_not_fully_qualified",
        },
        "source_ISA_assay_level_platform_and_same_aliquot_qualification": "PENDING",
        "source_feature_identity_diagnostic": {
            "exact_identical_original_feature_keys_both_protein_panels":
                exact_key_overlap,
            "casefold_identical_feature_keys_both_panels_if_no_collisions":
                casefold_overlap,
            "casefold_identifier_collisions_present": norm_collision,
            "source_identical_key_is_not_verified_same_protein_analyte": True,
            "cross_platform_analyte_and_proteoform_mapping_proven": False,
        },
        "identical_donor_split_and_input_representation_all_panels": True,
        "outer_folds": f"{repeats}x5_person_disjoint",
        "target_and_predictor_preprocessing": "outer_training_only",
        "fixed_total_PC_budgets": list(CAPACITIES),
        "fixed_ridge_alpha_values": list(ALPHAS),
        "protein_profile_target_selection": "top_80_train_only_variance_per_panel_per_fold",
        "panel_results": summaries,
        "original_same_person_relative_skill_sensitivity_IPP_minus_proteins":
            relative_sensitivity,
        "different_protein_panels_MAE_not_directly_interchangeable": True,
        "posthoc_source_replay_exploratory_not_preregistered": True,
        "bootstrap_is_person_residual_only_not_model_refit_or_multiplicity_adjusted": True,
        "raw_individual_donor_ids_and_molecular_rows_emitted": False,
        "original_germline_full_PGx_or_human_drug_PK_measured": False,
        "no_platform_caused_difference_or_functional_CYP_activity_inference": True,
    }


def live_original():
    with TemporaryDirectory(prefix="eatris_public_protein_panel_") as td:
        path = Path(td) / "temporary_source.h5mu"
        download_original(path)
        return benchmark_proteome_panels(path, verify_source=True)


if __name__ == "__main__":
    try:
        print(json.dumps(live_original(), sort_keys=True, indent=2))
    except (EATRISSourceAbstention, EATRISBenchmarkAbstention, ValueError, OSError) as exc:
        print(json.dumps({
            "status": "ABSTAIN_UNQUALIFIED_OR_UNAVAILABLE_PUBLIC_EATRIS_PROTEIN_PANELS",
            "reason_category": str(exc)[:190],
            "original_cross_panel_result_available": False,
        }))
        raise SystemExit(2)
