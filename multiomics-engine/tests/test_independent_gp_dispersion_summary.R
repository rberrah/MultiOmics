#!/usr/bin/env Rscript
# Package-free falsification evidence contract. Fixtures are synthetic
# test objects, not native DESeq2 output or independent biological evidence.
source("multiomics-engine/rnaseq_independent_gamma_poisson_null.R")
source("multiomics-engine/rnaseq_gp_dispersion_summary.R")
source("multiomics-engine/tests/test_independent_gamma_poisson_null.R")
scenarios<-c("gp_balanced_heterogeneous","gp_depth_heterogeneous",
             "gp_unbalanced_batch")
arms<-c("map_final","gene_wise_CR","fitted_trend","true_nb_alpha")
dat<-do.call(rbind,lapply(seq_along(scenarios),function(k){
  scenario<-scenarios[[k]]
  do.call(rbind,lapply(5001:5032,function(rep){
    src<-independent_gp_null_study(scenario,rep)
    flag<-as.integer(c(rep%%8L==0L,rep%%9L==0L,
                       rep%%10L==0L,rep%%11L==0L))
    data.frame(scenario=scenario,replicate=as.integer(rep),
      seed=as.integer(src$seed),arm=arms,
      source_features=400L,modeled_features=400L,
      finite_p=400L,finite_q=400L,
      false_calls_bh05=flag,any_false_call=flag,
      native_BH_replay_ok=TRUE,
      dispersion_log_rmse=c(.25,.31,.42,0),
      mean_wald_stat_squared=c(1.06,1.08,1.12,1),
      raw_p_under_005=c(18L,20L,21L,16L),
      raw_p_under_001=c(2L,3L,4L,1L),
      stringsAsFactors=FALSE)
  }))
}))
s<-gp_dispersion_independent_summary(dat)
stopifnot(nrow(dat)==384L,nrow(s)==12L,
   setequal(unique(s$scenario),scenarios),
   all(s$whole_null_studies==32L),
   all(s$source_gamma_poisson_independently_sampled),
   all(s$generator_includes_off_trend_per_gene_alpha),
   all(!s$generic_research_FDR_claim_permitted),
   all(s$individual_clopper_pearson_low<=s$empirical_global_null_FDR),
   all(s$individual_clopper_pearson_high>=s$empirical_global_null_FDR))
reject<-function(x)inherits(try(gp_dispersion_independent_summary(x),
                             silent=TRUE),"try-error")
stopifnot(reject(dat[-1L,,drop=FALSE]))
a<-dat;a$seed[[1L]]<-999L;stopifnot(reject(a))
a<-dat;a$replicate[[1L]]<-4001L;stopifnot(reject(a))
a<-dat;a$arm[[1L]]<-"fake_oracle";stopifnot(reject(a))
a<-dat;a$scenario[[1L]]<-"balanced_base";stopifnot(reject(a))
a<-dat;a$native_BH_replay_ok[[1L]]<-FALSE;stopifnot(reject(a))
a<-dat;a$modeled_features[[1L]]<-399L;stopifnot(reject(a))
cat("PASS: independent source 384-study/arm labels and false-calibration negative controls\n")
