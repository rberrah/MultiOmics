#!/usr/bin/env Rscript
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
for(scenario in c("balanced_base","balanced_depth_dispersion",
                  "unequal_partial_batch")) {
  a <- rnaseq_oracle_null_study(scenario,replicate=2001L)
  b <- rnaseq_oracle_null_study(scenario,replicate=2001L)
  stopifnot(identical(a$counts,b$counts),
    identical(a$depth,b$depth),
    identical(a$meta,b$meta),
    a$planted_condition_log2fc==0,
    nrow(a$counts)==40L,ncol(a$counts)==400L,
    identical(rownames(a$counts),rownames(a$meta)),
    !anyDuplicated(a$meta$condition),
    !anyDuplicated(colnames(a$counts)))
  if(scenario=="unequal_partial_batch")
    stopifnot(identical(as.integer(table(a$meta$condition)),c(12L,28L)))
}
stopifnot(inherits(try(rnaseq_oracle_null_study("not_original",2),
  silent=TRUE),"try-error"))
stopifnot(inherits(try(rnaseq_oracle_null_study("balanced_base",0),
  silent=TRUE),"try-error"))
set.seed(999L)
state <- get(".Random.seed",envir=.GlobalEnv)
rnaseq_oracle_null_study("balanced_base",2002L)
stopifnot(identical(state,get(".Random.seed",envir=.GlobalEnv)))
rows <- do.call(rbind,lapply(
  c("balanced_base","balanced_depth_dispersion","unequal_partial_batch"),
  function(s) do.call(rbind,lapply(1:32,function(rep)
    data.frame(scenario=s,replicate=rep,seed=780000L+rep,
      method=c("DESeq2_default","DESeq2_oracle_true_depth"),
      any_false_discovery=as.integer(rep %% 8L==0),
      finite_q=400L,size_factor_log_rmse=c(.08,0))))))
summary <- rnaseq_null_oracle_summary(rows,replicates=32L)
stopifnot(nrow(summary)==6L,
  all(summary$empirical_global_null_fdr==.125),
  all(summary$independent_null_studies==32L),
  all(!summary$scientific_calibration_certified))
rows <- rows[-1L,,drop=FALSE]
stopifnot(inherits(try(rnaseq_null_oracle_summary(rows,32L),
  silent=TRUE),"try-error"))
cat("PASS: deterministic fixed null generator, matched default-vs-oracle complete cohorts and noncertifying Monte Carlo uncertainty\n")
