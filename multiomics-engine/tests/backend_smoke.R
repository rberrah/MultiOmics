#!/usr/bin/env Rscript
source(file.path("multiomics-engine", "server.R"))

stopifnot(exists("apply_ms_qc_reference"))
stopifnot(exists("reference_preprocess_matrix"))
stopifnot(exists("run_backend_analysis"))
stopifnot(exists("run_deseq2_counts"))
stopifnot(exists("run_edger_ql_counts"))
stopifnot(exists("run_voom_counts"))
stopifnot(exists("compare_rnaseq_methods"))
stopifnot(exists("filter_informative_features"))
stopifnot(exists("make_keepx_grid"))
stopifnot(exists("normalise_keepx"))
stopifnot(exists("safe_perf_summary"))
stopifnot(exists("run_diablo_blocks"))
stopifnot(exists("prepare_multiblock_integration"))
stopifnot(exists("batch_completely_confounded"))
stopifnot(exists("validate_external_predictions"))

# -------------------------------------------------------------------------
# 0) Package-free guardrail helpers are deterministic
# -------------------------------------------------------------------------
helper_x <- cbind(
  constant = rep(1, 6),
  varying = 1:6,
  sparse = c(NA, NA, NA, NA, 1, 2)
)
rownames(helper_x) <- paste0("H", 1:6)
filtered <- filter_informative_features(helper_x, min_observed=3L)
stopifnot(identical(colnames(filtered), "varying"))
stopifnot(identical(make_keepx_grid(1L), 1L))
stopifnot(all(make_keepx_grid(12L) >= 1L))
stopifnot(all(make_keepx_grid(12L) <= 12L))
normalised_keepx <- normalise_keepx(NULL, list(a=matrix(1, 4, 7), b=matrix(1, 4, 3)), 2L)
stopifnot(length(normalised_keepx$a) == 2L)
stopifnot(length(normalised_keepx$b) == 2L)
stopifnot(all(normalised_keepx$a <= 7L), all(normalised_keepx$b <= 3L))

# DIABLO must keep repeated validation conservative and transparent even when
# mixOmics itself is not installed in the lightweight CI job.
stopifnot(as.integer(formals(run_diablo_blocks)$nrepeat) == 5L)
diablo_perf_mock <- list(error.rate=list(BER=c(0.40, 0.20, 0.30, NA_real_)))
diablo_perf_summary <- safe_perf_summary(diablo_perf_mock)
stopifnot(identical(diablo_perf_summary$status, "available"))
stopifnot(diablo_perf_summary$ber_values_extracted == 3L)
stopifnot(abs(diablo_perf_summary$ber_median - 0.30) < 1e-12)
stopifnot(abs(diablo_perf_summary$ber_mean - 0.30) < 1e-12)
stopifnot(abs(diablo_perf_summary$ber_min - 0.20) < 1e-12)
stopifnot(abs(diablo_perf_summary$ber_max - 0.40) < 1e-12)
stopifnot(grepl("same cohort", diablo_perf_summary$note, fixed=TRUE))
diablo_body <- paste(deparse(body(run_diablo_blocks)), collapse="\n")
stopifnot(grepl("max(5L, nrepeat)", diablo_body, fixed=TRUE))
stopifnot(grepl("performance_seed", diablo_body, fixed=TRUE))
stopifnot(grepl("candidate_keepX", diablo_body, fixed=TRUE))
stopifnot(grepl("external validation remains required", diablo_body, fixed=TRUE))

# A completely confounded technical series must block multiblock integration
# rather than being residualised together with the biological contrast.
conf_meta <- data.frame(
  sample_id=paste0("S",1:6),
  condition=c("A","A","A","B","B","B"),
  timepoint="T0",
  batch=c("BA","BA","BA","BB","BB","BB"),
  stringsAsFactors=FALSE
)
rownames(conf_meta) <- conf_meta$sample_id
conf_block <- matrix(seq_len(24), nrow=6, dimnames=list(conf_meta$sample_id, paste0("F",1:4)))
conf_guard <- prepare_multiblock_integration(
  blocks=list(transcriptomics=conf_block, proteomics=conf_block + 1),
  metas=list(transcriptomics=conf_meta, proteomics=conf_meta),
  covariates=character(),
  target_columns="condition"
)
stopifnot(identical(conf_guard$status, "blocked"))
stopifnot(grepl("completely confounded", conf_guard$message, fixed=TRUE))

# A non-confounded technical series may be adjusted deterministically.
bal_meta <- conf_meta
bal_meta$batch <- c("B1","B2","B1","B2","B1","B2")
bal_guard <- prepare_multiblock_integration(
  blocks=list(transcriptomics=conf_block, proteomics=conf_block + 1),
  metas=list(transcriptomics=bal_meta, proteomics=bal_meta),
  covariates=character(),
  target_columns="condition"
)
stopifnot(identical(bal_guard$status, "ok"))
stopifnot(all(vapply(bal_guard$details, function(x) x$status %in% c("adjusted","not_needed","not_estimable"), logical(1))))

# -------------------------------------------------------------------------
# 1) Reference MS QC: blank filter + pooled-QC LOESS + RSD + explicit MNAR
# -------------------------------------------------------------------------
assays <- c("BL1","BL2", paste0("QC",1:5), paste0("S",1:6))
orders <- seq_along(assays)
sample_type <- c("blank","blank", rep("pooled_qc",5), rep("biological",6))
meta <- data.frame(
  subject_id = assays,
  sample_id = assays,
  assay_id = assays,
  omic = "metabolomics",
  condition = c(rep("",7), rep(c("A","B"), each=3)),
  timepoint = "T0",
  batch = "B1",
  technical_replicate = "1",
  outcome = "",
  survival_time = "",
  survival_event = "",
  sample_type = sample_type,
  injection_order = orders,
  stringsAsFactors = FALSE
)

x <- cbind(
  MS_TRUE_SIGNAL = c(
    1,1,
    80,90,100,110,120,
    95,100,105,115,120,125
  ),
  MS_BLANK_CONTAM = c(
    10,11,
    18,19,20,21,22,
    20,21,19,22,20,21
  ),
  MS_UNSTABLE_QC = c(
    0.2,0.3,
    10,100,12,110,9,
    40,42,41,43,39,44
  ),
  MS_MNAR = c(
    0.1,0.1,
    50,51,49,50,52,
    60,NA,62,64,65,66
  )
)
rownames(x) <- assays

protocol <- list(
  msBlankFilter = "yes",
  msBlankFold = 5,
  msQcRsdFilter = "yes",
  msQcRsdThreshold = 0.30,
  msDriftCorrection = "yes",
  msMnarStrategy = "left_censored"
)
qc <- apply_ms_qc_reference(x, meta, protocol)

stopifnot(nrow(qc$matrix) == 6L)
stopifnot(!"MS_BLANK_CONTAM" %in% colnames(qc$matrix))
stopifnot(!"MS_UNSTABLE_QC" %in% colnames(qc$matrix))
stopifnot("MS_TRUE_SIGNAL" %in% colnames(qc$matrix))
stopifnot("MS_MNAR" %in% colnames(qc$matrix))
stopifnot(all(is.finite(qc$matrix[,"MS_MNAR"])))
stopifnot(qc$summary$blank_filtered_features >= 1L)
stopifnot(qc$summary$qc_rsd_filtered_features >= 1L)
stopifnot(qc$summary$drift_corrected_features >= 1L)
stopifnot(qc$summary$mnar_imputed_values >= 1L)

# -------------------------------------------------------------------------
# 2) DESeq2/voom concordance summary is deterministic and package-free
# -------------------------------------------------------------------------
deseq_mock <- data.frame(
  feature=c("G1","G2","G3","G4"),
  stat=c(5,-4,0.5,2),
  padj=c(0.001,0.01,0.8,0.04),
  stringsAsFactors=FALSE
)
voom_mock <- data.frame(
  feature=c("G1","G2","G3","G4"),
  t=c(4.5,-3.5,0.2,1.8),
  adj.P.Val=c(0.002,0.02,0.9,0.08),
  stringsAsFactors=FALSE
)
concordance <- compare_rnaseq_methods(deseq_mock, voom_mock, alpha=0.05)
stopifnot(identical(concordance$status, "ok"))
stopifnot(concordance$common_features == 4L)
stopifnot(concordance$significant_overlap == 2L)
stopifnot(is.finite(concordance$spearman_statistics))
stopifnot(concordance$sign_agreement_among_any_significant == 1)
stopifnot(grepl("not averaged", concordance$interpretation, fixed=TRUE))

# -------------------------------------------------------------------------
# 3) Browser-to-R payload routing, without requiring heavy packages in CI
# -------------------------------------------------------------------------
metadata_csv <- paste(
  "subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate",
  "P1,P1,R1,transcriptomics,A,T0,B1,1",
  "P2,P2,R2,transcriptomics,A,T0,B1,1",
  "P3,P3,R3,transcriptomics,B,T0,B1,1",
  "P4,P4,R4,transcriptomics,B,T0,B1,1",
  "P1,P1,P1A,proteomics,A,T0,B1,1",
  "P2,P2,P2A,proteomics,A,T0,B1,1",
  "P3,P3,P3A,proteomics,B,T0,B1,1",
  "P4,P4,P4A,proteomics,B,T0,B1,1",
  sep="\n"
)

rna_csv <- paste(
  "feature_id,R1,R2,R3,R4",
  "GENE1,5,5.2,7.5,7.8",
  "GENE2,3,3.1,3.2,3.0",
  sep="\n"
)
protein_csv <- paste(
  "feature_id,P1A,P2A,P3A,P4A",
  "P001,4,4.1,5.5,5.7",
  "P002,2,2.2,2.1,2.0",
  sep="\n"
)

payload <- list(
  protocol = list(
    organism="human",
    objective="groups",
    longitudinal=FALSE,
    designType="independent",
    outcomeType="none",
    covariateColumns=list()
  ),
  dataTypes = list(
    transcriptomics="log_expression",
    proteomics="log_intensity"
  ),
  identifierTypes = list(
    transcriptomics="gene_symbol",
    proteomics="uniprot"
  ),
  columnMapping = list(
    subject_id="subject_id",
    sample_id="sample_id",
    assay_id="assay_id",
    omic="omic",
    condition="condition",
    timepoint="timepoint",
    batch="batch",
    technical_replicate="technical_replicate"
  ),
  metadataCsv = metadata_csv,
  matrices = list(
    transcriptomics=rna_csv,
    proteomics=protein_csv
  )
)

result <- run_backend_analysis(payload)
stopifnot(identical(result$status, "ok"))
stopifnot(identical(result$engine$version, "1.4.0"))
stopifnot(all(c("transcriptomics","proteomics") %in% names(result$preprocessing)))
stopifnot("transcriptomics_differential" %in% names(result$methods))
stopifnot("proteomics_differential" %in% names(result$methods))
stopifnot(result$methods$transcriptomics_differential$status %in% c("ok","unavailable"))
stopifnot(result$methods$proteomics_differential$status %in% c("ok","unavailable"))
stopifnot(is.list(result$packages))
stopifnot("edgeR" %in% names(result$packages))

# -------------------------------------------------------------------------
# 3b) One layer is valid for reference differential inference and browser PCA
# -------------------------------------------------------------------------
single_payload <- payload
single_payload$matrices <- list(transcriptomics=rna_csv)
single_payload$dataTypes <- list(transcriptomics="log_expression")
single_result <- run_backend_analysis(single_payload)
stopifnot(identical(single_result$status, "ok"))
stopifnot(identical(single_result$analysisMode, "single_omic"))
stopifnot(identical(single_result$analysedLayers, "transcriptomics"))
stopifnot(!"diablo" %in% names(single_result$methods))
stopifnot(!"mofa2" %in% names(single_result$methods))
stopifnot(single_result$methods$transcriptomics_differential$status %in% c("ok","unavailable"))
single_payload$protocol$objective <- "explore"
single_explore <- run_backend_analysis(single_payload)
stopifnot(identical(single_explore$status, "ok"))
stopifnot(identical(single_explore$methods$single_omic_exploration$status, "not_applicable"))
stopifnot(!"mofa2" %in% names(single_explore$methods))

# -------------------------------------------------------------------------
# 4) Frozen external validation is evaluated without model fitting
# -------------------------------------------------------------------------
external_demo <- data.frame(
  outcome=c(0,0,0,1,1,1,0,1),
  prediction=c(0.05,0.20,0.25,0.72,0.85,0.93,0.30,0.78)
)
external_result <- validate_external_predictions(
  external_demo,
  outcome_type="binary",
  prediction_kind="probability",
  independent_cohort=TRUE,
  cohort_label="CI independent cohort",
  bootstrap_repetitions=50L,
  seed=20260928L
)
stopifnot(identical(external_result$status, "external_validation"))
stopifnot(identical(external_result$evaluation_status, "ok"))
stopifnot(is.finite(external_result$metrics$auc))
stopifnot(external_result$metrics$auc > 0.9)
stopifnot(external_result$bootstrap$valid_repetitions > 0L)

cat("multiomics reference R backend smoke: PASS\n")
