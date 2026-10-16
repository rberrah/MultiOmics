# Genuine Nutrimouse per-animal frozen-score study — narrow research qualification

**Status: candidate scientific CI gate for issue #20, not independent biological or clinical validation.**

## Measured source and experimental design

Official public mixOmics Nutrimouse source at immutable Git commit `ef3e760526623d9e91e945e4b50be48d764efc40`, `data/nutrimouse.rda`: 40 animals measured for 120 preselected liver macroarray gene-expression signals and 21 liver fatty-acid percentages by gas chromatography. These are neither RNA-seq raw counts nor LC-MS peak areas. It is the **same published single cohort** used in the existing DIABLO/MOFA2 benchmarks; cross-cohort replication is not established.

One animal = one biological subject. RNA and lipid assays share the original animal `subject_id`, but receive **distinct specimen identifiers**, because a shared physical aliquot is not documented. They are NOT two independent subjects. Genotype and diet are real experimental metadata; sampling/dose time and any PK concentration remain unknown.

## Locked fit and independent check

1. Predeclare the **exact** original feature `Cyp3a11` (case-insensitive lookup only), without screening for p-values or choosing a substitute. Fail if absent/ambiguous.
2. Predeclare a phenotype-blind split: test indices 4,8,12,…,40 (10 animals), and remaining 30 animals for training. Repeated assays always stay within their own animal.
3. Fit only the **training** mean and sample SD of measured single-gene expression, then freeze. Score = `(individual measured expression − training mean) / training sample SD`. No genotype, diet, test values, imputation, model-selection, or PK outcomes enter fitting.
4. Save real per-assay heldout scores only, marked `training_frozen_test_projection`. Archive the 30 training IDs/values, ten heldout IDs, exact R mean/SD, feature specification, source and model SHA-256.
5. Independently recompute training mean/SD and the ten heldout scores in Node from the full original matrices; check train/test subject nonoverlap and file SHA-256; generate the normal production exchange package containing the **full real 120-feature transcriptome and 21-feature lipidome**, paired sample sheet, and optional ten scores. Validate independently with the strict Python interchange reader.
6. Perform a heldout-only perturbation: altering test expression must not alter frozen training fit. Retain 199 exploratory genotype-label permutations *after scoring*, including negative results; these do not provide confirmatory FDR or external validation.

## Claim boundaries

A successful run establishes **a real per-assay score, fitted with training-only preprocessing and projected to genuine untouched animals**, plus technical output reproducibility. It is **not** a clinical CYP3A4/5 phenotype, patient PK covariate, molecular biomarker, cross-study validation, proof of multi-omic superiority, or improvement in PopPK/PD. The 10-animal holdout has low precision and diet/genotype may be associated. The original 120 genes were preselected in the published study, and original platform measurement QA was not recreated.

The generated exchange intentionally retains `validation.status=exploratory` and `individual_signature_export.status=schema_checked_only_not_biologically_validated` even when independent numerical checks pass. The last scientific step in issue #20 remains an independent cohort and meaningful multifeature biological signature and ultimately *actual PK* benchmarking downstream.

## Running the pinned public study

From the pinned Bioconductor R4.5.2/Bioc3.22 container already used in the Nutrimouse DIABLO workflow:

```bash
Rscript multiomics-engine/tests/test_real_nutrimouse_frozen_individual_score.R \
  external/mixOmics tmp/multiomics/real-nutrimouse-frozen-score
node scripts/export_real_nutrimouse_frozen_score.mjs \
  tmp/multiomics/real-nutrimouse-frozen-score
python3 scripts/validate_multiomics_exchange.py \
  tmp/multiomics/real-nutrimouse-frozen-score/exchange
```

R source and frozen original mixOmics checkout are controlled by `.github/workflows/nutrimouse-diablo-validation.yml`. CI archives genuine model hashes, independent matrices, heldout scores, and all adverse results. No data from PMx Explain, MIPD-Engine, or Omics2PK are modified.
