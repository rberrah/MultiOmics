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
  !nrow(null_reference_gross_failures(result)))
# Even zero out of twelve events does not establish a precise 5% null FDR.
zero <- transform(fixture,
  discoveries_q05=0L,false_discovery_fraction=0L)
zero_summary <- null_reference_summary(zero)
stopifnot(all(zero_summary$mc_exact_ci_low==0),
          all(zero_summary$mc_exact_ci_high>0.20),
          all(zero_summary$empirical_global_null_fdr==0))

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
inflated <- result
inflated$empirical_global_null_fdr[[1L]] <- 0.375
stopifnot(nrow(null_reference_gross_failures(inflated))==1L)
cat("RNA-seq global-null exact-binomial uncertainty: PASS (no scientific certification)\n")
