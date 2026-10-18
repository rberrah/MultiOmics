#!/usr/bin/env Rscript
# Base-R tests: a passing benchmark must not be mistaken for calibrated 5% FDR.
source("multiomics-engine/null_reference_metrics.R")

fixture <- do.call(rbind, lapply(c("balanced_base","unequal_partial_batch"),
  function(scenario) {
    do.call(rbind, lapply(c("DESeq2","edgeR_QL","limma_voom"),
      function(method) {
        data.frame(scenario=scenario,replicate=seq_len(12L),
          method=method,seed=seq_len(12L)+
            if (scenario=="balanced_base") 100L else 1000L,
          features_tested=400L,finite_q=400L,
          discoveries_q05=rep(c(1L,0L),c(1L,11L)),
          false_discovery_fraction=rep(c(1L,0L),c(1L,11L)),
          stringsAsFactors=FALSE)
      }))
  }))
result <- null_reference_summary(fixture)
stopifnot(nrow(result)==6L,
  all(result$replicates==12L),
  all(result$datasets_with_any_false_discovery==1L),
  all(abs(result$empirical_global_null_fdr-1/12)<1e-12),
  all(result$mc_exact_ci_low<0.05),
  all(result$mc_exact_ci_high>0.05),
  all(result$mc_exact_ci_high-result$mc_exact_ci_low>0.20),
  all(result$simultaneous_ci_low<=result$mc_exact_ci_low),
  all(result$simultaneous_ci_high>=result$mc_exact_ci_high),
  all(result$simultaneous_family_cells==6L),
  all(result$simultaneous_family_confidence==0.95),
  all(result$simultaneous_inconclusive),
  all(!result$scientific_certified),
  !nrow(null_reference_gross_failures(result)))
# Even zero out of twelve events does not establish a precise 5% null FDR.
zero <- transform(fixture,
  discoveries_q05=0L,false_discovery_fraction=0L)
zero_summary <- null_reference_summary(zero)
stopifnot(all(zero_summary$mc_exact_ci_low==0),
          all(zero_summary$mc_exact_ci_high>0.20),
          all(zero_summary$empirical_global_null_fdr==0),
          all(zero_summary$simultaneous_ci_high>=zero_summary$mc_exact_ci_high),
          all(!zero_summary$scientific_certified))

reject <- function(bad, min_replicates=8L) {
  inherits(try(null_reference_summary(bad,min_replicates),silent=TRUE),"try-error")
}
bad <- fixture[-1,]; stopifnot(!reject(bad))
bad <- fixture[fixture$replicate<=7L,]; stopifnot(reject(bad))
bad <- fixture; bad$replicate[2L] <- bad$replicate[1L]; stopifnot(reject(bad))
bad <- fixture; bad$seed[2L] <- bad$seed[1L]; stopifnot(reject(bad))
bad <- fixture; bad$discoveries_q05[2L] <- 500L; stopifnot(reject(bad))
bad <- fixture; bad$false_discovery_fraction[1L] <- 0L; stopifnot(reject(bad))
bad <- fixture; bad$finite_q[1L] <- 401L; stopifnot(reject(bad))
bad <- fixture; bad$scenario[1L] <- NA_character_; stopifnot(reject(bad))
bad <- fixture; bad$false_discovery_fraction[1L] <- NA_real_; stopifnot(reject(bad))
bad <- fixture; bad$method[1L] <- ""; stopifnot(reject(bad))
bad <- fixture; bad$replicate[1L] <- 1.5; stopifnot(reject(bad))
# Even with 128 all-zero simulations, a simultaneous upper 95% bound
# does not establish an exact 5% FDR (and cannot be pooled across methods).
many_zero <- do.call(rbind,lapply(c("DESeq2","edgeR_QL","limma_voom"),
  function(method) data.frame(
    scenario="balanced_base",replicate=seq_len(128L),method=method,
    seed=seq_len(128L)+20000L,features_tested=400L,finite_q=400L,
    discoveries_q05=0L,false_discovery_fraction=0L)))
many_summary <- null_reference_summary(many_zero,min_replicates=128L)
stopifnot(nrow(many_summary)==3L,
  all(many_summary$replicates==128L),
  all(many_summary$simultaneous_family_cells==3L),
  all(many_summary$simultaneous_ci_high>0),
  all(!many_summary$scientific_certified))
# Artificial gross anti-calibration must be visible; do not turn this into a
# passing 5% certification because the software workflow ran successfully.
bad_null <- many_zero
bad_null$false_discovery_fraction[1:60] <- 1L
bad_null$discoveries_q05[1:60] <- 1L
bad_summary <- null_reference_summary(bad_null,min_replicates=128L)
stopifnot(any(bad_summary$simultaneous_inflation_detected),
          any(bad_summary$simultaneous_ci_low>0.05))
stopifnot(inherits(try(null_reference_summary(
  fixture,simultaneous_family_alpha=0),silent=TRUE),"try-error"))
inflated <- result
inflated$empirical_global_null_fdr[[1L]] <- 0.375
stopifnot(nrow(null_reference_gross_failures(inflated))==1L)
cat("RNA-seq global-null exact-binomial uncertainty: PASS (no scientific certification)\n")
