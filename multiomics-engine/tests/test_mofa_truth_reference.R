#!/usr/bin/env Rscript
# True MOFA2 HDF5/mofapy2 fitting, not a PCA surrogate. Generates synthetic
# causal/truth factors solely for POST-FIT assessment. Factor order/sign
# ambiguity is accommodated by maximal ABS correlation.
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/validation_benchmarks.R")
if(!requireNamespace("MOFA2",quietly=TRUE))stop("Actual Bioconductor MOFA2 is required.")
seeds <- c(20261009L,20261010L)
scenarios <- c("biological","biological_plus_batch","batch_only")
records <- list()
work <- file.path("tmp","multiomics","mofa-truth")
dir.create(work,recursive=TRUE,showWarnings=FALSE)
for(scenario in scenarios)for(seed in seeds) {
  simulated <- simulate_mofa_known_truth(seed,scenario)
  dirname <- file.path(work,paste(scenario,seed,sep="_"))
  fit <- run_mofa2_blocks(
    simulated$blocks,dirname,factors=3L,
    seed=seed,convergence_mode="fast",
    min_samples=60L,max_features_per_block=40L
  )
  stopifnot(
    identical(fit$summary$matrixOrientation,
        "feature_rows_subject_columns_for_MOFA2"),
    identical(fit$summary$implementationContract,
        "mofa2-feature-rows-sample-columns-v2"),
    isTRUE(fit$summary$subjectIdentityPreserved),
    fit$summary$samples==nrow(simulated$truth)
  )
  scores <- MOFA2::get_factors(fit$model,as.data.frame=TRUE)
  audit <- mofa_factor_truth_diagnostics(scores,simulated$truth)
  records[[length(records)+1L]] <- data.frame(
    seed=seed,scenario=scenario,
    maxBiological=unname(audit$maxBiological),
    maxTechnical=unname(audit$maxTechnical),
    biologicalRecovered=audit$biologicalRecovered,
    technicalDetected=audit$batchSignalDetected,
    technicalDominated=audit$batchDominated,
    pass=if(scenario=="biological") {
      audit$biologicalRecovered
    } else if(scenario=="biological_plus_batch") {
      audit$distinctTwoFactorRecovery
    } else {
      !audit$biologicalRecovered && audit$batchSignalDetected
    },stringsAsFactors=FALSE
  )
}
result <- do.call(rbind,records)
utils::write.csv(result,file.path(work,"truth_results.csv"),row.names=FALSE)
writeLines(c(
  "ACTUAL MOFA2 synthetic factor recovery evidence",
  "Two deterministic seeds x biological, mixed bio+batch and technical-only scenarios",
  "Absolute correlations account for sign / factor ordering, not formal factor identifiability.",
  "Failures remain in truth_results.csv: never drop a failed seed/scenario.",
  "No proof of MNAR, clinical generalization, cross-study robustness, or optimal factor count."
),file.path(work,"methodological_scope.txt"))
print(result)
if(nrow(result)!=6L || anyNA(result$pass) || !all(result$pass))
  stop("MOFA2 did not recover prespecified strong factors or misclassified technical-only factors. See truth_results.csv")
cat("MOFA2 real engine synthetic biological/batch/null factor benchmark: 6/6 scenarios PASS\n")
