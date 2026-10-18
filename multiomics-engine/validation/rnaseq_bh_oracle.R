#!/usr/bin/env Rscript
# Prespecified BH algorithmic reference: NOT native DESeq2/edgeR/limma FDR.
# All subjects and features are synthetic; no patient's measurements are used.
source("multiomics-engine/reference_inference_audit.R")
args <- commandArgs(trailingOnly=TRUE)
output <- if(length(args)>=1L) args[[1L]]
          else Sys.getenv("MULTIOMICS_VALIDATION_OUTPUT",
                          "tmp/reference-methods")
if(length(args)>1L || !nzchar(output))
  stop("Provide at most one output path.")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
fixed_seed <- 20261018L
oracle <- run_independent_bh_oracle(replicates=2048L,
                                    features=200L,seed=fixed_seed,
                                    alpha=0.05)
utils::write.csv(oracle$per_study,
  file.path(output,"bh_oracle_independent_null_per_study.csv"),
  row.names=FALSE)
utils::write.csv(oracle$summary,
  file.path(output,"bh_oracle_independent_null_summary.csv"),
  row.names=FALSE)
writeLines(c(
  "schema: multiomics.bh-oracle-integrity.v1",
  paste("source_sha:",Sys.getenv("GITHUB_SHA","local_unverified")),
  paste("seed:",fixed_seed),
  "studies_per_mechanism: 2048",
  "features_per_study: 200",
  "null_model: independent Uniform(0,1) p-values; no fitted NB count model",
  "negative_control: independent Beta(0.5,1) anticonservative p-values",
  "adjustment: stats::p.adjust(p, method='BH')",
  "reference_fdr: 0.05; algorithmic regression interval +-0.03 (not a proof of statistical calibration)",
  "negative_control_alarm: empirical FWER/FDR under all-null > 0.20",
  "two_sided_MC_interval: exact binomial, independent synthetic datasets",
  "true_DESeq2_edgeR_limma_FDR_calibration: NOT assessed by this oracle",
  "biological_validation: NOT assessed",
  capture.output(sessionInfo())
),file.path(output,"bh_oracle_independent_null_manifest.txt"))
print(oracle$summary,row.names=FALSE)
cat("PASS: prespecified independent-BH oracle and anticonservative failure control; no native FDR certification\n")
