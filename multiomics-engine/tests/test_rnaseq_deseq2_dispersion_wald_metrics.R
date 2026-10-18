#!/usr/bin/env Rscript
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
source("multiomics-engine/rnaseq_deseq2_dispersion_wald_metrics.R")
scenarios <- c("balanced_base","balanced_depth_dispersion",
               "unequal_partial_batch")
methods <- c("native_default","oracle_true_dispersion",
             "oracle_true_depth_and_dispersion")
rows<-do.call(rbind,lapply(seq_along(scenarios),function(k) {
  s<-scenarios[[k]]
  do.call(rbind,lapply(3001:3032,function(i)
    data.frame(scenario=s,method=methods,replicate=i,
      seed=as.integer(780000L+k*10000L+i),
      fdr_any_discovery=c(i%%8==0,i%%10==0,i%%12==0)*1L,
      bh_discoveries=c(i%%8==0,i%%10==0,i%%12==0)*1L,
      finite_p=400L,finite_q=400L,
      wald_stat_squared_mean=c(1.15,.99,1.01),
      dispersion_log_rmse=c(.2,0,0),
      depth_factor_log_rmse=c(.05,.05,0),
      native_bh_integrity_verified=TRUE)))
}))
x<-rnaseq_dispersion_wald_summary(rows)
stopifnot(nrow(x)==9L,
  all(x$independent_global_null_datasets==32L),
  all(x$exact_binomial_ci95_low<=x$empirical_global_null_fdr),
  all(x$exact_binomial_ci95_high>=x$empirical_global_null_fdr),
  all(!x$mechanism_causally_identified),
  all(!x$universal_FDR_control_certified))
reject<-function(x)inherits(try(rnaseq_dispersion_wald_summary(x),
  silent=TRUE),"try-error")
stopifnot(reject(rows[-1L,,drop=FALSE]))
x<-rows;x$replicate[[1L]]<-3002L;stopifnot(reject(x))
x<-rows;x$seed[[1L]]<-9L;stopifnot(reject(x))
x<-rows;x$finite_q[[1L]]<-100L;stopifnot(reject(x))
x<-rows;x$fdr_any_discovery[[1L]]<-2L;stopifnot(reject(x))
x<-rows;x$native_bh_integrity_verified[[1L]]<-FALSE;stopifnot(reject(x))
x<-rows;x$dispersion_log_rmse[[1L]]<-NA_real_;stopifnot(reject(x))
x<-rows;x$bh_discoveries[[1L]]<-2L;stopifnot(reject(x))
st<-rnaseq_oracle_null_study("balanced_base",3001L)
stopifnot(length(st$dispersion)==ncol(st$counts),
  all(st$dispersion>0),st$planted_condition_log2fc==0)
cat("PASS: fixed 3×3 DESeq2 native true-dispersion whole-study null metrics; forged/skipped studies fail closed\n")
