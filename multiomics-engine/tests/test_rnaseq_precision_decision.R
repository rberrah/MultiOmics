#!/usr/bin/env Rscript
# Pure-base-R contract test for the fixed-horizon report: synthetic fixtures,
# NEVER scientific evidence. Executes the real report on temporary CSVs.
out <- tempfile("multiomics_precision_contract_")
dir.create(out,recursive=TRUE)
on.exit(unlink(out,recursive=TRUE),add=TRUE)
sc_null <- c("balanced_base","balanced_depth_dispersion",
             "unequal_partial_batch")
sc_mixed <- c("one_pct","five_pct","twenty_pct")
methods <- c("DESeq2","edgeR_QL","limma_voom")
null <- list()
mixed <- list()
for(s in seq_along(sc_null)) {
  for(rep in 1001:1128) {
    seed <- 780000L+s*10000L+rep
    for(method in methods) {
      any_false <- as.integer(rep %% 20L == 0L)
      null[[length(null)+1L]] <- data.frame(
        scenario=sc_null[[s]],replicate=rep,method=method,seed=seed,
        features_tested=400L,finite_q=400L,
        discoveries_q05=any_false,false_discovery_fraction=any_false)
    }
  }
}
for(s in seq_along(sc_mixed)) {
  n_true <- c(4L,20L,80L)[[s]]
  for(rep in 1001:1064) {
    seed <- 620000L+s*10000L+rep
    for(method in methods) {
      discoveries <- as.integer(rep %% 12L == 0L)
      mixed[[length(mixed)+1L]] <- data.frame(
        scenario=sc_mixed[[s]],method=method,replicate=rep,seed=seed,
        false_discovery_proportion=0, power=discoveries/n_true,
        selected_features=discoveries,false_discoveries=0L,
        tested_features=400L,total_planted_features=400L,
        n_true=n_true,mean_absolute_log2fc_error=0.2,
        signed_log2fc_bias=0.01)
    }
  }
}
null <- do.call(rbind,null)
mixed <- do.call(rbind,mixed)
np <- file.path(out,"bioconductor_null_replicates.csv")
mp <- file.path(out,"bioconductor_mixed_truth_per_study.csv")
utils::write.csv(null,np,row.names=FALSE)
utils::write.csv(mixed,mp,row.names=FALSE)
script <- "multiomics-engine/validation/rnaseq_precision_decision.R"
run <- function() {
  suppressWarnings(system2("Rscript",c(script,shQuote(out)),
    stdout=TRUE,stderr=TRUE))
}
success <- run()
stopifnot(is.null(attr(success,"status")) || attr(success,"status")==0L,
  file.exists(file.path(out,"RNASEQ_MONTE_CARLO_PRECISION_REPORT.md")),
  file.exists(file.path(out,"rnaseq_precision_null_simultaneous_95.csv")),
  file.exists(file.path(out,"rnaseq_precision_mixed_95.csv")))
summary <- utils::read.csv(
  file.path(out,"rnaseq_precision_null_simultaneous_95.csv"))
stopifnot(nrow(summary)==9L,
  all(summary$replicates==128L),
  all(summary$simultaneous_family_cells==9L),
  all(!summary$scientific_certified))
# A missing study must cause an explicit failure, NEVER a shorter horizon.
utils::write.csv(null[-1L,,drop=FALSE],np,row.names=FALSE)
failed <- run()
stopifnot(!is.null(attr(failed,"status")) && attr(failed,"status")!=0L)
utils::write.csv(null,np,row.names=FALSE)
# Replicating one study under two methods with conflicting seeds is prohibited.
tampered <- mixed
tampered$seed[[1L]] <- 99L
utils::write.csv(tampered,mp,row.names=FALSE)
failed <- run()
stopifnot(!is.null(attr(failed,"status")) && attr(failed,"status")!=0L)
cat("PASS: fixed-horizon precision report, independent seeds and fail-closed missing-study controls; synthetic fixture only\n")
