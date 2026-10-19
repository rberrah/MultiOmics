# MultiOmics → downstream methods: exchange contract v1 (2026-10-15)

**Status: implemented and independently consumer-tested.** The contract is generated inside the research ZIP and can also be downloaded directly from the results page as `multiomics_omics2pk_exchange_v1.zip`. The producer normalizes compatible CSV/TSV/semicolon inputs, mapped metadata columns and transposed matrices into one strict canonical exchange; it does not change the original archived source files. This interchange layer does not imply that network analyses, biomarkers or PK effects are validated.

## Directory layout

```text
exchange/
  manifest.json
  samples.csv
  features.csv
  matrices/
    transcriptomics.csv
    proteomics.csv
    metabolomics.csv
```

Only include matrices actually supplied. One layer is valid. No network inference is run implicitly.

### samples.csv (one row per assay/measurement)

Required columns: `assay_id,sample_id,subject_id,omic`.
Recommended: `tissue,collection_time,collection_time_unit,condition,batch,technical_replicate_id,source_dataset_id`.
The MultiOmics #20 interoperability extension **appends** `timepoint,timepoint_unit,time_origin` at the end of the existing V1 sample sheet. These optional columns may be blank in older V1 archives and are not a breaking schema revision.

- `timepoint` preserves the **source's categorical original label** (e.g. `T0`, `T12`); it is **not** an elapsed time or PK sampling clock.
- `timepoint_unit` and `time_origin` are exported only if explicitly present in the original sample sheet. Missing means unknown, not zero hours or time since dose.
- `collection_time` and `collection_time_unit` retain the original quantitative measurement and unit when supplied; MultiOmics **never** derives 12 h or any unit from `T12`. PK temporal alignment must be qualified downstream.
- Source `technical_replicate` is copied into the V1 contract's `technical_replicate_id` unless the explicit V1 name is already provided. Multiple assays from one specimen must not be interpreted as additional subjects.
- Cross-omic assays that share `sample_id` must not disagree on subject, tissue or any supplied time fields. The JS producer and offline Python validator both reject conflicting original metadata.
- The browser manifest records the actual source Git checkout SHA **when available at build time**; a build from a source archive without Git metadata records `unknown` rather than a fabricated hash. This is the **producer** commit, not the revision of a model fitted from omic signatures.
- **No individual signatures by default**: `exchange/matrices/*.csv` are individual original measurements, not independently qualified frozen-model scores. Group logFC/q-values cannot be promoted to individual PK covariates. A future optional `individual_signatures.csv` requires true per-assay scores, frozen held-out projection or out-of-fold training provenance, with model and feature SHA-256 verified by the Omics2PK consumer.

- `subject_id` is the biological unit (patient, animal, cell line); `sample_id` is the specimen; `assay_id` is the **unique matrix column key**. Never treat repeated assays from one subject as independent subjects.
- Pairing and temporal order must be explicit; unknown tissue/time are **empty with missingness documented**, never guessed.
- Use pseudonymous IDs, not names or direct patient identifiers.

### Optional individual biological signatures: `individual_signatures.csv`

This file is **absent by default**. It is only included when an actual per-assay scoring result is explicitly supplied to `buildResearchBundle` as `result.individualSignatures`; the exporter never computes scores from group differences, `q_value`, differential `logFC`, or an unnormalized matrix. Missing or invalid explicit scores cause **export failure**, not a silently partial ZIP.

Required fields per score: `subject_id,sample_id,assay_id,omic,signature_id,score,score_scale,fit_scope,reference_model_sha256,feature_definition_sha256,source_exchange_analysis_id,score_provenance`.

- `subject_id/sample_id/assay_id/omic` must exactly match an existing matrix's biological identity.
- `score` must be a finite, genuinely per-assay numeric scoring result; `score_scale` is explicit and cannot be a group statistic or raw counts.
- `fit_scope` must be `training_frozen_test_projection` or `out_of_fold_training`; in-sample and transductive scores are rejected.
- Both SHA-256 values are required and syntactically checked. **The exporter cannot independently attest that the referenced model was really fitted exclusively on training subjects**. The original model artifact, training IDs/split and frozen preprocessing still need independent custody and audit by Omics2PK.
- `source_exchange_analysis_id` must match the exact source exchange manifest; `score_provenance` is `native_r_individual_score` or `external_frozen_projection` only.
- The file has a verified SHA-256 in `manifest.files.individual_signatures` and is checked by the offline Python validator and the pinned consumer Omics2PK integration gate.
- The presence of a structurally valid signature is **not** evidence of a biological biomarker, cross-cohort generalization or improvement in PopPK/PD. The manifest retains exploratory validation status.

Existing ZIP consumers that ignore the optional `files.individual_signatures` property remain supported. Data remain local in the research ZIP; public CI only uses synthetic fixtures.

### features.csv (one row per matrix row)

Required: `omic,feature_id,id_namespace`. Recommended: `display_name,organism,taxon_id,annotation_database,annotation_version,tissue_scope,mapping_status`.
`(omic,feature_id)` must be unique. Keep source IDs without auto-conversion: Ensembl, UniProt and ChEBI are different namespaces. Mapping to another namespace belongs in a separate provenance-bearing crosswalk; one-to-many mappings cannot be silently collapsed.

### matrices/*.csv

Header `feature_id,<assay_id>,...`; one row per feature, **all** analyzed features (not top 100), and measurement values on the declared scale. CSV NA is literally `NA` for missing values, not zero. No pre-filtered significance-only matrix, no imputed values mislabelled as raw.
Every matrix column must have exactly one matching `samples.csv` row for the same omic; every matrix feature must have a matching features row for the same omic. The converse may include explicitly excluded features, but exclusions require a QC record.

### manifest.json

See [JSON Schema](./DOWNSTREAM_EXCHANGE_V1.schema.json). The manifest enumerates all files with relative paths and SHA-256, the original dataset access conditions, sample-unit definition, values scale/normalization, producer commit, exact runtime/dependency versions **or explicit nulls when unavailable**, and per-method qualification. Neither a green pipeline nor `passed_technical` means inferential or biological confirmation.

**No inferred tissue/time** and **no invented lock versions**. For controlled/private data, local archiving and opt-in transfer only; clinical data must not enter public CI artifacts.

## Validation and compatibility

Run `python3 scripts/validate_multiomics_exchange.py exchange/` (Python standard library only). Checks file hashes, required columns, referential integrity, duplicate IDs, non-finite values, row counts and missingness. This is an **interchange validation**, not a biological benchmark.

The browser dossier attempts an optional `exchange/` directory from the archived **original** source files and writes `downstream_exchange.status` plus an aggregate Omics2PK-readiness summary in its outer manifest. The dedicated Omics2PK download exposes the same contract at the ZIP root for direct consumption. Source metadata may use researcher-mapped column names and CSV/TSV/semicolon delimiters; compatible transposed matrices are normalized to `feature_rows_sample_columns`. Empty numeric source cells become explicit `NA`, never zero. For incompatible layouts the research dossier records a reason rather than misrepresenting the data.

The independent Omics2PK consumer gate verifies hashes, subject/specimen/assay identity, matrix values, categorical time preservation and optional frozen/out-of-fold individual-score syntax. This remains a **technical source-exchange qualification**, not evidence that an exported molecular feature improves measured PK.

## Scientific limits retained after Research v1

Research v1 freezes R 4.5.2 / Bioconductor 3.22, the immutable container digest, the full R package inventory, Node/package locks and the qualification workflows. The exchange contract itself is therefore no longer waiting for a generic runtime freeze.

Three scientific boundaries remain deliberately open:

- **#20 — individual signatures → Omics2PK:** schema/provenance checks and consumer compatibility are qualified; biological transferability and incremental PK value still require genuine subject-matched PK/PGx data and held-out/external validation.
- **#44 — DESeq2 calibration:** native DESeq2 execution is reproducible, but null simulations showed excess false discoveries in several tested generators. A q-value is not a universal 5% FDR certificate.
- **#47 — LC-MS analytical validation:** authentic raw Sacurine xcms execution and QC are qualified, including unfavourable QC observations; chemical identity, standards/interlaboratory accuracy and a universally valid feature-exclusion policy remain outside the claim.

No claim of universally controlled FDR, cross-cohort score transfer, metabolite identity or clinical PK benefit is supported merely by this contract.

## Exact source-based replay of R results (when R actually ran)

The previous `run_exchange_analysis.R` was intentionally **removed**: it fitted a different statistical workflow and would not reproduce the reported results. Do not use ad-hoc recipes as substitutes for the original computational pipeline.

For a completed **local R reference backend** analysis, the researcher ZIP contains:

- `reproduction/reference_request.json`: the exact JSON payload POSTed to the R backend, including experimental design, contrasts, original files and options.
- `reproduction/reference_expected.json`: the original R API response, as returned in the analysis.
- `reproduction/multiomics-engine/*.R`: snapshots of the R functions that actually executed the calculations.
- `reproduction/run_reference_backend.R`: a replay runner that sources these same R functions and calls `run_backend_analysis()` on the exact request.
- `reproduction/README.md`: explicit setup, limitations and reproducibility instructions.

```bash
# Run from the extracted ZIP root with the original R/Bioconductor package versions:
Rscript reproduction/run_reference_backend.R reproduction/
```

This is a **source-exact replay candidate**, not yet a tested guarantee of matching numeric output. The original R environment must be restored and replayed, with numerical comparisons and frozen Reactome mappings when those were used. If R was not used successfully, no R replay is included.

### Figures and exploratory JavaScript

The Svelte figures are drawn in JavaScript, but the **primary numerical source is explicit per layer**. When a complete native-R table passes the strict completeness and identity contract, the displayed feature-wise effects/p/q values are projected from that table and the figures consume those same projected effects without refitting. Otherwise the layer remains browser-exploratory. Mixed-source analyses are labelled as such.

The ZIP can contain `figures_browser/*.svg` for the currently displayed vector figures. These are saved browser output images, not figures regenerated by the R runner. Therefore `browser_svg_snapshots.r_replayable = false` and `exact_r_numeric_and_figure_equivalence_certified = false` remain explicit even when the underlying primary numerical table is native R.

### Beginner workflow

1. Select your research question and upload an omic matrix, then run the guided workflow.
2. Download the researcher ZIP, or use the dedicated **Omics2PK exchange** download when the goal is downstream pharmacometrics. Each primary table records whether it is a complete native-R table or a browser-exploratory table.
3. If the local R backend was used, locate `reproduction/` and run its `run_reference_backend.R` with the exact package versions. Compare the replayed results to `reference_expected.json` and the full R CSVs before citing reproducibility.
4. Retain `figures_browser/` as vector snapshots; do not assert R recreated them.
5. For future omics-to-PK/PD work, join biological subject, specimen, assay, collection time and tissue by explicit identifiers. Never join by row position or mistake repeated measurements for independent subjects.

No PK/PD fit or network analysis is automatically run by this export.

## Reproducing figures in R, without re-estimating the statistics

In the downloaded research dossier, run **one command from its extracted root**:

```bash
Rscript run_all.R
```

This script (1) attempts replay of the original native R calculations if they actually ran and the exact versions/source SHA-256 match, (2) verifies the numerical comparison using the bundled independent Python script, and (3) creates similar publication figures directly with R based on the archived display-ready values.

Generated files include `figures_r/output/*.pdf`, SVGs when the R build supports Cairo, the underlying `*_plot_values.csv` and `figure_manifest.json`. A status report at `research_reproduction_status.json` separates numerical R replay from graph generation. **A browser-only run does not become a statistically equivalent R analysis.**

The R-generated heatmaps, metabologram and central-carbon map are intended to convey the same measured information as the web versions. Font metrics, interactivity and SVG pixel layout may differ, and any non-displayed interactive network state remains outside the exported figure set. R generates these figures using the actual archived values, without refitting statistical models.

For future PK/PD integration, the validated subject/sample/assay/tissue/time identifiers remain preserved and networks are still optional.
