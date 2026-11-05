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

The analysis will **not** claim improved clinical PK,
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

## 4. EATRIS public vs inaccessible components

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
