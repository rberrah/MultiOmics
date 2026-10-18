#!/usr/bin/env Rscript
# Synthetic fixtures are NOT biological evidence.
source("multiomics-engine/rnaseq_deseq2_dispersion_component_metrics.R")
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
scenarios<-c("balanced_base","balanced_depth_dispersion",
             "unequal_partial_batch")
arms<-c("map_final","gene_wise_CR","fitted_trend","true_nb_alpha")
dat<-do.call(rbind,lapply(seq_along(scenarios),function(k) {
  sc<-scenarios[[k]]
  do.call(rbind,lapply(4001:4032,function(id){
    study<-rnaseq_oracle_null_study(sc,id)
    data.frame(scenario=sc,replicate=id,seed=study$seed,arm=arms,
      source_features=400L,modeled_features=400L,
      finite_p=400L,finite_q=400L,
      false_calls_bh05=as.integer(c(id%%8==0L,id%%9==0L,
                                id%%10==0L,id%%11==0L)),
      any_false_call=as.integer(c(id%%8==0L,id%%9==0L,
                                id%%10==0L,id%%11==0L)),
      native_BH_replay_ok=TRUE,
      dispersion_log_rmse=c(.2,.3,.25,0),
      mean_wald_stat_squared=c(1.1,1.2,1.1,1.0),
      raw_p_under_005=c(18L,20L,19L,17L),
      raw_p_under_001=c(2L,3L,2L,1L),
      stringsAsFactors=FALSE)
  })))
}))
s<-rnaseq_deseq2_dispersion_component_summary(dat)
stopifnot(nrow(s)==12L,all(s$whole_null_studies==32L),
          all(s$individual_clopper_pearson_low<=s$empirical_global_null_FDR),
          all(s$individual_clopper_pearson_high>=s$empirical_global_null_FDR),
          all(!s$universal_FDR_calibration_certified))
reject<-function(d)inherits(try(
  rnaseq_deseq2_dispersion_component_summary(d),
  silent=TRUE),"try-error")
stopifnot(reject(dat[-1L,,drop=FALSE]))
x<-dat;x$seed[[1L]]<-999L;stopifnot(reject(x))
x<-dat;x$replicate[[1L]]<-4002L;stopifnot(reject(x))
x<-dat;x$arm[[1L]]<-"fabricated_pipeline";stopifnot(reject(x))
x<-dat;x$modeled_features[[1L]]<-399L;stopifnot(reject(x))
x<-dat;x$dispersion_log_rmse[[1L]]<-NA_real_;stopifnot(reject(x))
x<-dat;x$native_BH_replay_ok[[1L]]<-FALSE;stopifnot(reject(x))
x<-dat;x$raw_p_under_001[[1L]]<-21L;stopifnot(reject(x))
x<-dat;x$false_calls_bh05[[1L]]<-2L;stopifnot(reject(x))
cat("PASS: 12-cell four-arm real DESeq2 dispersion plan rejects missing and false evidence\n")
