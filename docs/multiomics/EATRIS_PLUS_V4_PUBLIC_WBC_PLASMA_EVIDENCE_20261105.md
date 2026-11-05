# EATRIS-Plus v4: original white-blood-cell vs plasma multiomics

**Date:** 2026-11-05. **Scope:** strictly original public
Zenodo [record 19136013, version 4, published 2026-03-20](https://zenodo.org/records/19136013);
actual individual matrices are accessed transiently in source-replay CI
without committing or printing individual participants.

**Official, primary original source-data reproduction:**
[GitHub Actions #42036020160 SUCCESS](https://github.com/rberrah/MultiOmics/actions/runs/42036020160).
The original **141.7 MB** `EATRIS_PLUS_mudata.h5mu`
download matched its official MD5
`bbf0342797adcbcd6a3c5b242ca64573`.
The h5mu file was opened directly with h5py and its 12
`/mod/<assay>/obs/_index` arrays matched by original exact
sample key, not row order, guessed specimen labels or
source names. Only aggregate counts were produced.

## 1. Original public assay sizes (actual source HDF5)

| Original MuData modality name | Source-identified physical fraction (original Zenodo description; detailed method still requires ISA audit) | n source distinct subjects | n source quantified feature columns |
|---|---|---:|---:|
| `mRNA` | **White blood cells**; mRNA sequencing | **90** | **16,386** |
| `miRNA-mature` | White blood cells; miRNA sequencing | **89** | **2,130** |
| `EM-seq` | White blood cells; DNA methylation assay | **117** | **100,000** |
| `miRNA-qRT-PCR` | **Plasma**; miRNA qPCR | **86** | **179** |
| `proteins` | Plasma-associated protein assay (do not assume same platform as `proteomics-ipp`) | **125** | **1,393** |
| `proteomics-ipp` | Plasma-associated proteomic assay (distinct original matrix) | **123** | **6,831** |
| `acylcarnitines` | Plasma; targeted metabolomics | **125** | **33** |
| `aminoacids` | Plasma; targeted metabolomics | **125** | **46** |
| `fattyacids` | Plasma; targeted metabolomics | **125** | **5** |
| `lipidomics-neg` | Plasma; lipidomics negative mode | **125** | **164** |
| `lipidomics-pos` | Plasma; lipidomics positive mode | **125** | **196** |
| `lipidomics-spline` | Plasma; lipidomics postprocessed view, actual processing yet to source-qualify | **124** | **257** |

Important caveat: source type assignments above are the repository
deposit's **study-level method descriptions**, not a demonstration
that all individual HDF5 modality labels contain exact ISA
sample-source crosswalk. In particular `proteins` and
`proteomics-ipp` require **separate** source-platform
qualification. We neither call circulating proteins intracellular
hepatocyte CYP, nor replace WBC mRNA with liver RNA.

## 2. Complete-case loss is already observed, not assumed

The original HDF5 contains **125 global original subject keys**.
Using all **12** public assays would restrict the exact
same-subject intersection to **55 subjects**. This may exacerbate
selection bias and variance, and is **not** a scientifically
justified default just because 12 assays are available.

Selected fully source-key matched important intersections:
- **90** original persons share WBC `mRNA`, plasma
  `proteins`, `acylcarnitines`, and `aminoacids`.
- **89** original persons share WBC `mRNA` and
  `miRNA-mature` (white blood cell gene regulatory study).
- **63** original persons share WBC `mRNA` and
  plasma `miRNA-qRT-PCR` (not tissue-matched).
- **123** original persons share both `proteins` and
  `proteomics-ipp` (distinct assay platforms).

The HDF5 IDs establish the same source participant key, not
necessarily the same **physical plasma aliquot** or same
collection moment. Detailed reference metadata/collection timing
must still be checked in the ISA export.

## 3. Prospective predeclared *multiomic* benchmark (not PK)

The purpose is to **falsify the hypothesis that more omics
automatically improves biological prediction**, using the
original same-participant matrices.

**Outcome:** the genuinely source-measured plasma
`proteins` profile. This is **a laboratory cross-omic
prediction**, not a measured AUC, clearance, active CYP
protein, drug toxicity, or causal treatment effect.

**Source-matched persons:** n90 at intersection.
**Predictors on exactly those original people:**
WBC mRNA, plasma acylcarnitines + amino acids, or both,
with constant target-protein prediction and train-permuted
WBC RNA controls. **Never introduce post-treatment drug
outcomes as baseline covariates**.

Predeclare 3×5 **individual-disjoint** CV, train-only
missingness screening, variance-based feature selection,
standardization/PCA (max 3 PCs per predictor layer),
target-proteome standardized on outer-train donors,
fixed Ridge(alpha=10). Compare per-subject heldout macro
**MAE on the same selected plasma proteins** (original
protein features may vary between training folds; all arms
within a fold share identical targets). Use **paired
original-person bootstrap** for exploratory error
differences, without refitting or post-selection correction.
Test controls:
- no molecular source, train-mean proteome prediction;
- plasma targeted metabolites only;
- WBC mRNA only;
- WBC mRNA + plasma targeted metabolites;
- WBC mRNA donor-label correspondence **permuted in
  training only** (negative control).

The analysis does **not** claim improved clinical PK,
genuine patient pharmacogenetics, hepatic CYP activity,
or mechanistic secretion from a statistical correlation.
A negative result is scientifically admissible.

Reproducible original software:
`scripts/audit_eatrisplus_v4_original_h5mu.py` and
`scripts/benchmark_eatrisplus_v4_cross_fraction.py`.
Corresponding synthetic index/privacy/anti-leakage tests
are separate from source-real data in GitHub Actions.
Do not commit any source HDF5, original subject matrices
or direct original participant IDs.

## 4. Actually completed original public-source predictive validation

Two independent GitHub Actions reproductions downloaded the **exact
same published MD5-pinned original 141.7 MB MuData**, with all
numeric source donor-level measurements processed transiently.
Results below are genuinely source-executed, **not synthetic
examples or imputed values from a paper**.

1. [Original-source **fixed Ridge/PCA** (SUCCESS
   #42037048717)](https://github.com/rberrah/MultiOmics/actions/runs/42037048717):
   three train-only PCs **per predictor layer**,
   fixed Ridge `alpha=10`, 3×5 outer subject-disjoint
   folds. Plasma acylcarnitines and amino acids supply
   **six total PCs** versus three from WBC RNA.
2. [Original-source **capacity-matched/nested** model
   (SUCCESS #42038211111)](https://github.com/rberrah/MultiOmics/actions/runs/42038211111):
   exactly **three total PCs per molecular arm**,
   tested joint budgets RNA 1 + metabolites 2 and RNA 2 +
   metabolites 1. Alpha chosen by **3-fold inner original
   biological donor CV** using a fixed search grid
   `(1, 10, 100, 1000, 10000, 100000, 1000000)`,
   separately in each outer training partition.
   Train-only target protein feature screening, sample
   imputation, predictor scaling and PCA refit inside
   the inner folds. 3×5 outer person-disjoint folds
   and the exact same n90 original subjects.

**Observed outcome: outer-training-standardized plasma
protein-profile MAE**, lower is better, macro averaged
across **80 train-selected source protein targets** per
outer fold. **Not** drug clearance, patient dosing,
germline genotype or clinical PD. Protein target feature
sets may vary by outer fold but are **identical across
model arms within that fold**.

| Original n90 WBC/plasma arm | Fixed alpha=10, PCs per layer | Matched 3 total PCs + inner-only tuned alpha |
|---|---:|---:|
| Outer-training mean protein profile (no omics) | **0.760615** | **0.760615** |
| WBC mRNA only | 0.783413 | 0.760756 |
| Plasma targeted acylcarnitines+amino acids | 0.792970 | 0.761520 |
| WBC mRNA + plasma metabolites | **0.819299** | 0.761429 (RNA1 + plasma2) |
| WBC mRNA + plasma metabolites, alternate fixed view allocation | Not tested | 0.760720 (RNA2 + plasma1) |
| WBC mRNA, independently train-permuted negative control | 0.784368 | 0.760747 |

**Actual original subject-level paired exploratory contrasts:**

- Fixed alpha=10 joint view versus no omic:
  **ΔMAE +0.058684** [95% descriptive CI **+0.045340,
  +0.072923**]. Joint versus plasma metabolites only
  **+0.026329** [**+0.019070,+0.033997**].
- Matched total 3 PCs, inner-tuned:
  **WBC mRNA alone** versus no omic Δ **+0.000141**
  [**+0.000057,+0.000221**];
  **plasma metabolite alone** Δ **+0.000905**
  [**−0.000148,+0.001984**];
  **joint RNA1+plasma2** Δ **+0.000814**
  [**+0.000311,+0.001388**];
  **joint RNA2+plasma1** Δ **+0.000105**
  [**+0.000022,+0.000186**].
- Matched joint minus matched plasma metabolites alone:
  RNA1/plasma2 Δ **−0.000091**
  [**−0.001086,+0.000764**];
  RNA2/plasma1 Δ **−0.000800**
  [**−0.001862,+0.000176**]. Neither interval
  excludes zero; no robust WBC added value demonstrated
  over the plasma molecular layer.
- WBC RNA and the **train-permuted** WBC RNA control
  perform very similarly under both designs;
  their similarity is a *negative control* on stable
  cross-compartment predictive gain in this narrow
  modeling setup, **not a demonstration that there are
  no true WBC-to-plasma physiological pathways**.

**Key finite-data/software insight:** the initial fixed
model produces substantial error degradation as more
views are added. With shared dimensionality and inner-fold
regularization, the models nearly revert to the
no-omics baseline. The inner-fold best alpha is often
**100,000–1,000,000** for WBC-containing arms:
the fitting procedure effectively shrinks the
molecular predictions almost away.

This is exactly the distinction that Omics2PK should
demonstrate carefully:

> Additional information cannot intrinsically harm an
> optimal model permitted to ignore it. Its *estimation
> and exploitation* in a finite cohort can degrade
> out-of-sample prediction; sufficiently regularized
> learning may neutralize the harm rather than
> reveal a new useful physiological biomarker.

**Very important limitation:** the two columns change
**both** the total number/allocation of PCs **and**
the method of choosing the shrinkage; therefore this
experiment does **not separately identify the causal
contribution of regularization versus dimensionality**.
No standalone superiority, mechanism or true causal
null follows. The reported 95% intervals are
**subject-level residual bootstraps with original
frozen heldout predictions**, **not** refitted nested
model uncertainty, not multiplicity corrected and not
external cohort evidence. The 15 outer test instances
per subject from repeated CV are *correlated*; inference
uses the **90 unique original individuals** after
averaging repeats. Actual tissue-specific sample
origin remains as in source methods; a systemic protein
profile is not measured liver CYP activity.

**Scientific novelty boundary:** predicting cross-omic
reference profiles is a sensitivity control for
source-aware MultiOmics integration and a warning
against automatic all-omics fusion. It does **not**
answer the core Omics2PK clinical+full-PGx
drug-AUC/CL improvement question.

The exact predeclared artifacts and tests are:
`scripts/benchmark_eatrisplus_v4_cross_fraction.py`,
`scripts/test_eatrisplus_v4_cross_fraction_benchmark.py`,
`scripts/benchmark_eatrisplus_v4_capacity_matched.py`,
`scripts/test_eatrisplus_v4_capacity_matched.py`.
Source-real CI analyses output aggregate values and
original source hashes, **never participants' source
assay rows**.

## 5. EATRIS public vs inaccessible components

The current open v4 original source is substantial and
allows true transcriptome/proteome/metabolome integration
analyses **now**; it must not be described as blocked on
other future cohorts. However **complete WGS/clinical
phenotype and drug-specific longitudinal human PK are not
in the public version**; the separately access-controlled
ClinData records cannot be presumed available. This
precludes the final Omics2PK research question *"does
omics add to full germline PGx for clinical tacrolimus
AUC/CL/C0?"* from being answered here.

This public-only study can validate physical source
annotations, model robustness and false discoveries,
without importing any unavailable patient dataset.
