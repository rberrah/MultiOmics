# Multi-omics advanced R engine

This directory contains reference-package adapters for analyses that should use their actual R implementations rather than browser approximations carrying the same names.

The design goal is conservative automation: the backend may automate a method only when its required study structure is identifiable. It does not turn a statistically impossible design into an analysable one, infer statistical power from sample size alone, combine p/q-values from different packages as if they were independent evidence, or describe internal cross-validation as external validation.

## Browser engine versus reference implementations

The public browser tool accepts **one, two or three omics layers**, and performs deterministic QC, adjusted feature-wise models, repeated-measures modelling, exploratory PCA, nested cross-validated prediction and Reactome interpretation. Balanced multi-block PCA and cross-omics correlations run only when at least two layers are loaded; with one layer the PCA is explicitly single-omic. The R backend accepts one matrix for supported per-layer reference methods. MOFA2 and DIABLO require multiple layers and are never silently used for one-layer work. A separate browser-side feasibility audit records blockers, cautions and advanced-method eligibility before strong interpretation.

`advanced_methods.R` contains reference implementations for:

- raw RNA-seq differential analysis through **DESeq2**, **edgeR quasi-likelihood**, and **edgeR + limma-voom**;
- normalized expression/proteomics-style linear modelling through `limma`;
- repeated-measures mixed models through `lmerTest`, with an automatically attempted random time slope only when repeated observations can support it and deterministic fallback to a random-intercept model when necessary;
- ranked enrichment through `fgsea`;
- MOFA2 through the `MOFA2` package;
- DIABLO through `mixOmics::block.splsda`, with repeated internal cross-validation for sparsity/component tuning when feasible.

## Input contract

`blocks` must be a named list of numeric matrices in samples × features orientation. Sample IDs must be row names. Constant and insufficiently observed variables are removed before advanced integration rather than silently carried into the factor/discriminant model.

Example:

    blocks <- list(
      transcriptomics = rna_matrix,
      proteomics = protein_matrix,
      metabolomics = metabolite_matrix
    )

## MOFA2

    source("multiomics-engine/advanced_methods.R")
    run_mofa2_blocks(
      blocks = blocks,
      output_dir = "results/mofa2",
      factors = 5,
      seed = 20260924
    )

The adapter deliberately requires at least **16 subjects shared across the analysed blocks**. This is an eligibility guardrail, not a claim that 16 subjects guarantee adequate power. It removes non-informative features, records view scaling, saves factor/weight tables and, when available, variance explained.

Known technical effects must be addressed before interpreting MOFA factors. A latent factor is a covariance pattern, **not** evidence that a mechanism is causal or that a biomarker is validated. If technical structure is confounded with biology, the browser feasibility audit blocks strong interpretation instead of relabelling the factor as biological.

## DIABLO

The outcome must be categorical. Prefer a named vector keyed by sample ID.

    outcome <- c(S01="control", S02="control", S03="treated", S04="treated", ...)
    run_diablo_blocks(
      blocks = blocks,
      outcome = outcome,
      output_dir = "results/diablo",
      ncomp = 2,
      seed = 20260924,
      tune = TRUE,
      nrepeat = 5
    )

Methodological safeguards:

- at least three subjects in every class are required for the automated tuning route;
- folds are bounded by the smallest class;
- the adapter enforces at least **five repeated M-fold CV runs**, including when the browser bridge requests fewer repeats;
- `tune.block.splsda()` is used with repeated M-fold cross-validation;
- BER is the tuning criterion so a larger class does not dominate the error metric;
- the candidate `keepX` grid, final `keepX` sparsity settings and component count are recorded;
- tuning failures and `perf()` failures are retained explicitly rather than silently discarded;
- `mixOmics::perf()` is run with a deterministic seed distinct from the tuning seed;
- the compact performance summary reports median, mean, minimum and maximum extracted BER rather than highlighting only the best observed value;
- selected variables are reported as a supervised multivariate signature, not as independently validated biomarkers;
- **external validation remains required** for claims of generalisable prediction or biomarker performance.

The fitted model, tuning object, performance object, sample variates, selected variables and compact summary are saved when available.

Important interpretation boundary: tuning and the `perf()` diagnostic are performed on the same cohort. Repeated CV reduces sensitivity to one arbitrary split but does **not** create an independent validation cohort and should not be presented as unbiased external performance.

## Raw RNA-seq counts

DESeq2 remains available as a full count-model implementation:

    run_deseq2_counts(
      counts = counts_matrix,
      metadata = metadata,
      output_dir = "results/deseq2",
      design_formula = ~ batch + condition,
      contrast = c("condition", "treated", "control")
    )

A quasi-likelihood edgeR route is also available:

    run_edger_ql_counts(
      counts = counts_matrix,
      metadata = metadata,
      output_dir = "results/edger_ql",
      design_formula = ~ batch + condition
    )

The same design can be analysed independently with edgeR filtering/TMM normalization followed by limma-voom:

    run_voom_counts(
      counts = counts_matrix,
      metadata = metadata,
      output_dir = "results/voom",
      design_formula = ~ batch + condition
    )

Rank-deficient design matrices are rejected by the edgeR/voom/limma adapters rather than fitted with silently non-estimable coefficients.

The automatic browser bridge currently uses DESeq2 as the primary raw-count analysis and edgeR/limma-voom as an independent sensitivity analysis when available. `run_edger_ql_counts()` is exposed for reference workflows and future routing. Method agreement is summarised with `compare_rnaseq_methods()` using:

- Spearman correlation between signed test statistics;
- number of BH/FDR-significant genes in each method;
- overlap and Jaccard index of significant genes;
- sign agreement among genes declared significant by either method.

The purpose is sensitivity analysis. **q-values are never averaged, pooled or converted into an artificial consensus q-value.** A disagreement remains visible.

The request protocol may set `rnaCountMethod` to `deseq2`, `voom`, `both` or `auto` for the currently routed methods.

## Already normalized / log-scale matrices

    run_limma_matrix(
      matrix = expression_matrix,
      metadata = metadata,
      output_dir = "results/limma",
      design_formula = ~ batch + condition
    )

The model matrix is checked for full rank before fitting.

## Reference repeated-measures model

    run_lmer_matrix(
      matrix = expression_matrix,
      metadata = metadata,
      output_dir = "results/lmer",
      fixed_formula = "condition * time + batch",
      subject_column = "subject_id",
      random_slope = "auto"
    )

With `random_slope="auto"`, a subject-specific time slope is attempted only when enough subjects have at least three distinct observations. Singular or failed random-slope fits fall back to a random-intercept model. The structure actually used is recorded per result. This is preferable to forcing a complex random-effects structure that the data cannot identify.

## Ranked enrichment

    run_fgsea_ranked(
      ranks = signed_statistics,
      pathways = reactome_pathways,
      output_dir = "results/fgsea",
      seed = 20260924
    )

Pathway enrichment organises molecular evidence; it is not by itself proof that a pathway is activated, causal or clinically relevant. The measured/mappable feature universe and identifier coverage should be retained in the study report.

### Reactome release provenance

Whenever Reactome analysis is used in the browser, the tool also queries the Reactome database-version endpoint and attaches the returned release identifier to the final analysis object. The provenance is stored in three places so downstream exporters cannot easily drop it:

- `externalDatabaseProvenance.reactome`;
- `reactome.provenance`;
- `reproducibility.externalDatabases.reactome`.

The stored record includes the Reactome release, version endpoint, AnalysisService base URL, captured analysis token(s) and timestamp. If the version endpoint is unavailable, the pathway analysis remains usable but the provenance record is marked `unavailable` with the retrieval error. The failure does not silently invent a database version.

The downloadable analysis JSON therefore carries the external knowledge-base release used at analysis time. Reactome content may change independently of PMx Explain, so this information should be retained with manuscript or regulated-workflow archives.

## Automatic browser bridge

Install the backend dependencies once:

    Rscript multiomics-engine/install_backend_dependencies.R

Start the local backend:

    Rscript multiomics-engine/run_backend.R

The default address is:

    http://127.0.0.1:8787

In **Auto** mode the browser calls `/health`. If the backend is available, it sends the current sample-sheet contract and loaded matrices to `/run` and executes applicable reference methods. If the local backend is absent, browser-native deterministic analysis remains available. In **Require R** mode, analysis fails instead of silently falling back.

Current automatic routing includes (for one or more loaded layers, when applicable):

- raw RNA-seq counts + independent groups -> DESeq2 and edgeR/limma-voom when available, with explicit concordance summary;
- normalized/log-scale independent group comparisons -> limma;
- repeated/longitudinal designs -> lmerTest;
- unsupervised integration -> MOFA2, subject to adapter/browser eligibility guardrails;
- categorical supervised integration -> tuned/validated mixOmics DIABLO when the class structure supports it;
- signed reference statistics with compatible Ensembl/UniProt IDs -> fgsea using current Reactome mappings.

### Local security policy

`run_backend.R` now applies a second origin guard on top of the annotated Plumber routes. Browser requests are accepted only from the default PMx origin and common local Vite/preview origins, plus any explicitly declared origins in `PMX_MULTIOMICS_ALLOWED_ORIGINS` (comma-separated). Rejected origins receive HTTP 403 before an analysis endpoint executes. Responses are marked `Cache-Control: no-store`.

The launcher also refuses a non-loopback bind by default. A deliberate remote bind requires:

    PMX_MULTIOMICS_ALLOW_REMOTE=yes

This flag is **not** a security mechanism. It is only an acknowledgement that a remote deployment must still be placed behind TLS, authentication, request-size limits, rate limits and network/origin controls. The recommended default remains `127.0.0.1`.

GitHub Pages does not run R itself. The bridge must run locally or on a controlled server.

### Environment provenance

The launcher exposes a read-only `/environment` endpoint containing:

- exact R runtime/version and platform;
- installed versions of `plumber`, `jsonlite`, `DESeq2`, `edgeR`, `limma`, `lmerTest`, `fgsea`, `MOFA2` and `mixOmics`;
- R library search paths.

No submitted study data are returned by this endpoint. For a manuscript or regulated workflow, archive this environment manifest together with the analysis result, original inputs and their cryptographic checksums.

## Raw LC-MS / GC-MS preprocessing

The reference raw-MS entry point is deliberately restricted to open, standardised spectra readable by the Bioconductor stack. Vendor formats should be converted upstream to mzML, for example with ProteoWizard/msconvert when the appropriate vendor reader is available.

Use:

    Rscript multiomics-engine/run_raw_ms.R raw_ms_manifest.csv raw_ms_output raw_ms_parameters.json

Templates are shipped in:

- `static/multiomics/templates/raw_ms_manifest.csv`;
- `static/multiomics/templates/raw_ms_parameters.json`.

The reference pipeline performs centWave chromatographic peak detection, optional obiwarp retention-time alignment, peak-density correspondence/grouping and chromatographic gap filling, then exports a peak-area feature matrix plus feature definitions, sample sheet and processing manifest.

This is **feature detection, not compound identification**. An m/z-retention-time feature must not be presented as a named metabolite until a separate identification/annotation workflow supports that claim.

## Frozen independent prediction validation

External predictive validation is implemented as a separate evaluator rather than as another resampling mode of the development model:

    Rscript multiomics-engine/run_external_validation.R external_validation_predictions_binary.csv external_validation_binary.json external_validation_output

Templates are shipped in:

- `static/multiomics/templates/external_validation_predictions_binary.csv`;
- `static/multiomics/templates/external_validation_binary.json`.

The evaluator accepts already-frozen predictions and never performs feature selection, hyperparameter tuning or refitting on the validation cohort. It supports binary, continuous, count, survival and multiclass outcomes. When estimable, deterministic percentile bootstrap intervals are added without replacing failed resamples by model refitting.

The status `external_validation` is emitted only when cohort independence is explicitly asserted. The software cannot infer historical independence from the numbers themselves; provenance must establish that the validation cohort was not used for model development.

## Auditable browser preprocessing

Browser preprocessing is intentionally conservative and is recorded per omics layer in `qc.preprocessingAudit`. The audit reports the declared measurement type, observed value characteristics, input range, chosen output scale, pseudocount when applicable, median-centering status and any compatibility warnings.

Current browser rules are explicit:

- transcriptomic raw counts: library-size CPM followed by log2(CPM + 0.5) for deterministic **screening** only; publication-grade count inference is delegated to DESeq2/edgeR/voom;
- proteomic spectral counts: library-size normalisation followed by log2 transformation;
- TPM, LFQ intensity, peak area and absolute concentration: positive-scale log2 transformation with a recorded data-derived pseudocount;
- LFQ intensity and metabolomic peak area: sample-wise median centering on the log scale;
- declared log-scale matrices: kept on their supplied log scale;
- declared normalized or unknown matrices: kept as supplied rather than silently imposing another normalization.

Negative values declared as raw/spectral counts are rejected. Fractional count-like values are explicitly flagged for review. Positive abundance types containing negative values are retained as supplied and flagged as likely already transformed/centred.

The browser does **not** silently infer or apply ComBat, quantile normalisation, VST or another study-specific preprocessing choice. Batch is instead handled as an explicit design/nuisance variable when identifiable.

## Single-omics and advanced MS sample annotations

One-layer RNA, proteomics and metabolomics analyses share the same sample-sheet contract and scientific design checks as multi-omics. If the sheet includes more modalities than the uploaded data, summary denominators and statistical contrasts are scoped to the loaded modalities. For the exploratory-only beginner shortcut, a one-to-one independent biological subject/assay relationship must be explicitly confirmed by the user; a group assignment is never invented. Browser fallback is always identified in the result report, including when local R is unavailable.

For CSV/TSV peak-area input, native m/z, retention time, adduct and identification metadata are exported as a separate technical annotation CSV, never interpreted as measured sample intensities. SHA-256 hashes of the original vendor file and converted quantitative matrix can be retained with the report; a user-provided feature-to-ChEBI mapping is reported as unverified rather than a validated MS identification.

## Advanced MS sample annotations

For LC-MS/GC-MS workflows the shared sample sheet additionally supports:

- `sample_type`: `biological`, `pooled_qc`/`qc`, or `blank`;
- `injection_order`: numeric sequence position.

The browser and R QC routes can exclude technical injections from biological inference, flag blank-associated features, estimate pooled-QC drift along injection order, apply a pooled-QC RSD threshold, and optionally perform explicitly requested low-tail imputation for left-censored MNAR values.

### Blank handling policy

Blank/sample ratios are QC evidence, not biological truth. The default mode therefore **flags** features whose biological median is below the declared biological/blank ratio without deleting them. An explicit **remove** mode is available for laboratories whose validated SOP requires exclusion. The selected mode and threshold are recorded.

Left-censored imputation is likewise an explicit sensitivity assumption. Confirmatory work should retain a no-imputation sensitivity analysis.

## Reproducibility and reporting

Random seeds are explicit where stochastic algorithms are used. The browser exposes two portable exports:

- `multiomics-methods-report.md`: human-readable methods/interpretation record;
- `multiomics-reproducibility-manifest.json`: machine-readable snapshot of the scientific question, design/endpoint settings, non-file interface parameters, local input-file metadata, visible corrected evidence, warnings and interpretation limits.

The browser now computes **SHA-256 locally** over the exact uploaded file bytes and over a canonical representation of the sample table. These hashes are stored in `reproducibility.cryptographicInputs` and exported with the analysis. A deterministic FNV-1a fingerprint may still be retained in `inputManifest` for fast accidental-change detection, but it is not treated as a cryptographic substitute for SHA-256.

The public truth benchmark goes further: its external data repositories are pinned to immutable Git commit SHAs, each exported benchmark input is hashed with SHA-256, and the resulting source/file provenance is attached to `benchmark-report.json` and `benchmark-provenance.json`.

Reference-package outputs should retain the `/environment` package/runtime manifest when used for a manuscript or regulated workflow.

## Multiple-testing interpretation policy

The user interface distinguishes two evidence levels rather than using an ambiguous single notion of “significant”:

- **primary evidence:** `q ≤ 0.05`;
- **exploratory/suggestive evidence:** `0.05 < q ≤ 0.10`.

Raw p-values are not used alone to select features from a high-dimensional family of tests. Neither threshold replaces effect size, direction, uncertainty, study design or external validation. A q-value is not interpreted as the probability that one individual biological hypothesis is false.

## What automation cannot certify

The pipeline can detect many deterministic incompatibilities, prevent several common statistical errors and make method choices reproducible. It cannot automatically certify:

- that the biological sampling frame is unbiased;
- that all confounders were measured;
- that the chosen effect size is clinically meaningful;
- adequate statistical power without a prespecified effect/variance model;
- causal identification without an appropriate causal design;
- external validity without an independent target population/cohort.

Those limitations remain explicit rather than being hidden behind a successful software run.
