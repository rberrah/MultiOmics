# MultiOmics — exact reproducibility acceptance contract (15 October 2026)

**Scientific requirement:** The user downloads **the code actually used** to produce the numerical tables and the visible figures. A statistically similar independent implementation does **not** satisfy this requirement.

## Observed architecture

- `runDeterministicAnalysis()` generates browser-level statistical results in JavaScript.
- `buildMultiomicsVisualizationData()` plus Svelte components generate interactive figures in JavaScript/SVG, including heatmaps, metabolic network and metabolograms.
- `run_backend_analysis()` executes the actual R reference analyses for supported research designs. These are distinct results, not a replacement for browser numbers.
- The old `run_exchange_analysis.R` was removed because it silently substituted a different set of preprocessing/filtering/statistical choices.
- The new `reproduction/` bundle captures the **actual** JSON request, original response, R sources and a native replay script; a numerical comparator fails on changed results. Browser SVG snapshots are clearly labeled as not R-generated.

## Definition of done: researcher-facing claims

For each displayed result or figure, store an immutable provenance object:

`{artifact_id, producing_engine, producing_function, code_commit, input_sha256, processed_data_sha256, parameters, seed, dependencies_lock, external_database_release, output_sha256, verification_status}`

Distinguish `computed_in_browser`, `computed_in_r`, `rendered_in_browser`, and `rendered_in_r`. Never imply a JavaScript computation is an R fit. The claims **R replay verified**, **numerical results match**, **figure replay verified** must be independently checked and are false by default.

### Acceptance gate A — numbers

1. Preserve original raw inputs, metadata, contrasts, exclusions, units, covariates, model formulae, sample-level/block-level identities and all method options. Use the same preprocessing and algorithm source, not a look-alike recipe.
2. Replay in a clean environment restored from an R/Bioconductor and Python/JS lock (plus public dataset SHA and database releases). Record the actual runtime.
3. Compare each full R table (all features, IDs, NA/q, coefficient sign, sorting), QC metrics, counts and fitted statistics. A method failing due to a missing package is **not** counted as reproduced.
4. Report tolerances *before* running. Prefer exact categorical/ID/NA matches and tight numeric tolerances; document any BLAS/parallel nondeterminism. Do not relabel approximate results as exact.

### Acceptance gate B — figures

1. Choose **one canonical figure-generation engine** for publication figures. If the requirement is reproducible R code, R must generate each publication figure that the site displays as an R result. The website must display that exact returned SVG/PNG/PDF, not independently rebuild it with Svelte.
2. Include the exact figure source R, resolved input tables, themes/palettes, fonts policy, labels, units, axis limits, grouping, plot dimensions and file outputs. Never deliver an R heatmap that only resembles a different browser heatmap.
3. Compare exported SVG structural/text/data attributes or raster image pixels using declared rendering tolerances. Validate legends, label contrast, captions, missing-value colors and missing measurements.
4. Interactive browser-only exploration may remain for beginners, but it must be tagged **interactive preview / non-R figure** until a canonical export is generated.
5. Export all *publication figures* individually with source and metadata. A snapshot of the currently selected tab does **not** establish reproducibility of figures in hidden tabs.

### Acceptance gate C — supported methods

Validate DESeq2, edgeR QL, edgeR/limma-voom, limma, fgsea, MOFA2, DIABLO individually on independently sourced public data. Pair independent native reference calls and source-exact replay on each dataset. Retain failures and missing packages. MOFA2/DIABLO require repeatable seeds, training-only transformations and versioned factors/feature selections; fgsea requires **frozen pathways and identifier mapping**.

### Beginner workflow

1. Run an obvious demo; choose one omic or several; hide advanced settings by default.
2. See the result with method/figure provenance in plain language.
3. Download one ZIP containing the **actual** producing code and the data needed to replay it, plus a `run_all.R` or equivalent one-command entry point.
4. A report immediately states: `identical numeric replay: PASS/FAIL/NOT CHECKED` and `identical figure replay: PASS/FAIL/NOT CHECKED`. Do not hide a FAIL inside a green build.

### Future omics → PK/PD

The v1 subject/specimen/assay/time/tissue/biological-namespace contract is intended for future exposure-response modeling and longitudinal cohort joins. Keep it model-agnostic now. All network extensions are opt-in and separated from the validated primary analysis.

## Frozen R dependency closure for Research v1

The reference R environment is pinned at three levels: R 4.5.2 / Bioconductor 3.22, an immutable Bioconductor container digest, and the complete **171-package** installed-library inventory in `multiomics-engine/environment/reference-runtime-inventory.json`. The frozen-reference workflow restores CRAN packages from the 2026-10-16 snapshot and exact-version Bioconductor source archives where required, then rejects any changed, missing or extra package before accepting the reference runtime. This closes the previous gap where a fixed container could still install a moving dependency set.

## Current scope decision

**The stronger claim "identical R replay of numerical and visual site results"
remains BLOCKED and is explicitly outside Research v1.** Research v1 does not
claim that independent browser/Svelte graphics have byte-identical or
geometry-identical R renderings. It freezes exact provenance of browser
artifacts and complete native-R numerical tables where native inference is
actually supported, while keeping the producing engine visible.

This distinction is deliberate: a reproducible research release can preserve
the exact producing code and numbers without falsely asserting that two
different rendering engines produce the same SVG geometry. If a future release
wants the stronger identical-cross-engine figure claim, it must satisfy the
additional gates below under a new qualification scope.

## R-generated figure equivalents added in PR #110

The archive now includes `run_all.R` and `figures_r/generate_figures.R`.
The single-command entry point attempts the exact R reference replay where
the backend actually ran, invokes an independent numerical verification script,
then plots descriptive **R-native PDF/SVG figure equivalents** from the archived
browser display dataset. It records separate statuses in
`research_reproduction_status.json`.

The figure script uses the **same input values**, heatmap row z-scores,
sample order, effects, effect scales, and curated biochemical map annotations
already placed in `provenance/analysis_complete.json`. It additionally writes
the plotted values to `*_plot_values.csv`, permitting point-by-point audits.
It produces (when present) heatmaps, a descriptive circular metabologram,
and a central carbon schematic. This is scientifically similar visualization,
**not** a claim of identical SVG geometry, font metrics, all interactive tabs,
pathway enrichment, or a new inference.

### Remaining blockers for the stronger identical cross-engine replay claim

- Compare real public cohort outputs from the website and from a clean R
  re-execution (native DESeq2/edgeR/limma/fgsea/MOFA2/DIABLO separately).
- Freeze R/Bioconductor versions, dependencies, external identifier mappings
  and pathways; account for runtime stochasticity.
- If an R result is the *published* result, show that same R table in the UI
  rather than its independently calculated JavaScript exploratory counterpart.
- Validate R-generated figures on a public cohort, including captions, sample
  annotations, units and missing values; avoid silently claiming that a saved
  browser SVG was generated in R.
- Validate novice ZIP extraction, one-command replay and error messages on
  Windows and Linux.
