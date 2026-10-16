# MultiOmics — Public GSE52778 end-to-end replay evidence

**Date:** 2026-10-15  
**Source dataset:** Bioconductor `airway`, GEO GSE52778, Himes et al., PLOS ONE (2014).  
**Scientific boundary:** software/data-path reproducibility, *not* independent clinical validation, global FDR certification, or validation of all six MultiOmics analysis families.

## Study design and actual result

Four independent airway smooth-muscle cell lines were measured untreated and dexamethasone treated (eight paired assays). The native R backend uses a subject-blocked design (`~ subject_id + condition`) and an explicit treatment-vs-untreated contrast. Its prespecified RNA count-expression filter retained **17,199 genes**; all feature rows are preserved.

The checked native R model was independently compared with a direct Bioconductor DESeq2 invocation using the **same design, preprocessing rule, contrast and analysis cohort**. All six native DESeq2 numerical columns (`baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`) and missingness positions are checked by `multiomics-engine/tests/test_public_airway_paired_backend.R`.

The integration harness `scripts/test_public_airway_website_r_numeric_figures.mjs` then:

1. Projects every one of the **17,199 complete native R fitted rows** to the site's primary numerical result structure; checks feature IDs, log2FC, p and q for every row and checks a known glucocorticoid-responsive marker.
2. Builds the actual JavaScript display-data heatmap from the original public source count matrix and sample sheet, but using **native R model effects**, not independently inferred JavaScript statistical tests.
3. Generates a publication-style heatmap with base R from the archived display-ready JS visualization object. Compares **30 features × 8 samples** point by point (including missing values) against the values plotted by JavaScript.
4. Checks the R figure's declared z-score units, colors, missing-value treatment and the JS/R heatmap color legend (small RGB tolerance for interpolation).
5. Builds the actual browser ZIP export with the original public source inputs, the original reference R request and response, and byte-for-byte matching R source SHA-256. Extracts the ZIP without changing its contents and runs **`Rscript run_all.R`**.
6. Requires `PASS_R_REFERENCE_ONLY` in the replay report, figures written from the R source, and all full-feature numerical outputs to be available in the extracted research folder.

## Observed CI result

The initial real-data integration run completed successfully: [GitHub Actions #38792865652](https://github.com/rberrah/rberrah.github.io/actions/runs/38792865652).

It printed:

> PUBLIC AIRWAY PAIRED BACKEND: PASS; 17199 genes; four paired cell lines; CRISPLD2 positive; complete CSV and HTTP JSON; native DESeq2 matched

> PASS: Airway/GSE52778 17199 native R values reproduced feature-by-feature; browser & R figure plot data are identical for 30 features × 8 samples

The figure-palette/unit/NA assertions were subsequently strengthened and require the **newest** CI run to pass before being called confirmed. Generated SVG/PDF figures are **similar in scientific content**, not pixel-identical renderings.

## What researchers must know

- The browser's basic inference models can be **simplified or approximate relative to the appropriate native R reference methods**, depending on the study and feature type. A difference in numbers does not imply a rounding error; it can mean a different model or preprocessing pipeline. Browser-only numbers are exploratory.
- A genuine native R result should be projected to the site *as-is*, with its complete table downloadable; the R code replay must match the producing backend sources, runtime and model settings.
- Native R is more **methodologically appropriate** for supported inference (e.g. DESeq2 for raw RNA-seq counts), but R itself does not guarantee valid findings. Study design, independence, assumptions, batch effects, calibration and multiplicity must still be independently assessed.
- A heatmap built from the same data values is scientifically comparable, but interactive Svelte geometry and R plotting geometry are different. Exact visual pixel identity is not claimed.
- Airway validates a single, paired, public RNA-seq workflow. Separate public-gold tests remain needed for edgeR, limma, fgsea, MOFA2, DIABLO, real multi-omic figures, and other sample types.

## Beginner acceptance criterion

Extract the downloaded research ZIP, use the built-in `COMMENCER_ICI_FR_EN.md` guide and run from the extracted folder:

```sh
Rscript run_all.R
```

The standard replay and numerical verification are written in R (plus the `jsonlite` R package); no Python installation is required for the one-command workflow. An independent Python verifier remains available for an optional second opinion. Inspect `research_reproduction_status.json`. **PASS_R_REFERENCE_ONLY** means the native R calculation passed the original numerical replay comparison. **GENERATED_SIMILAR_PDF_SVG** means R generated graphical equivalents; it does not mean the site plots are identical pixel by pixel. Any FAIL or NOT_CHECKED is a limitation, not a silent pass.

Future omics–PK/PD analyses can reuse subject, sample, assay, tissue and time keys, but no PK/PD inference is validated by this benchmark.
