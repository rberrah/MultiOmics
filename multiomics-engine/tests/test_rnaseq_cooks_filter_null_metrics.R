#!/usr/bin/env Rscript
# Package-free adversarial schema/control tests, NOT native biological evidence.
source("multiomics-engine/rnaseq_cooks_filter_null_metrics.R")
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
source("multiomics-engine/rnaseq_independent_gamma_poisson_null.R")
original<-c("balanced_base","balanced_depth_dispersion",
            "unequal_partial_batch")
gp<-c("gp_balanced_heterogeneous","gp_depth_heterogeneous",
      "gp_unbalanced_batch")
arms<-c("cooks_on_if_on","cooks_on_if_off",
        "cooks_off_if_on","cooks_off_if_off")
rows<-do.call(rbind,lapply(c(original,gp),function(sc) {
  is_orig<-sc%in%original
  ids<-if(is_orig)6001:6032 else 5001:5032
  do.call(rbind,lapply(ids,function(rep) {
    s<-if(is_orig)rnaseq_oracle_null_study(sc,rep)
      else independent_gp_null_study(sc,rep)
    data.frame(
      scenario=sc,source_generator=if(is_orig)
        "original_smooth_NB" else "independent_heterogeneous_GP",
      replicate=rep,seed=s$seed,arm=arms,
      source_genes=400L,modeled_genes=400L,total_native_fit_count=1L,
      finite_p=c(399L,399L,400L,400L),
      finite_q=c(380L,399L,382L,400L),
      discovery_count=c(0L,0L,as.integer(rep%%8L==0L),
                        as.integer(rep%%8L==0L)),
      any_discovery=c(0L,0L,as.integer(rep%%8L==0L),
                      as.integer(rep%%8L==0L)),
      p_removed_by_cooks=c(1L,1L,0L,0L),
      q_removed_by_filter=c(19L,0L,18L,0L),
      raw_p_005=c(16L,16L,17L,17L),
      native_bh_replay_ok=TRUE,same_model_pvalue_check=TRUE,
      stringsAsFactors=FALSE)
  }))
}))
s<-rnaseq_cook_if_summary(rows)
stopifnot(nrow(rows)==768L,nrow(s)==24L,
    all(s$whole_null_studies==32L),
    all(!s$FDR_5pct_certified),
    all(s$individual_exact_ci95_low<=s$empirical_global_null_FDR),
    all(s$individual_exact_ci95_high>=s$empirical_global_null_FDR))
reject<-function(x)inherits(try(rnaseq_cook_if_summary(x),
  silent=TRUE),"try-error")
stopifnot(reject(rows[-1L,,drop=FALSE]))
x<-rows;x$arm[[1L]]<-"fake";stopifnot(reject(x))
x<-rows;x$seed[[1L]]<-123L;stopifnot(reject(x))
x<-rows;x$replicate[[1L]]<-6002L;stopifnot(reject(x))
x<-rows;x$source_generator[[1L]]<-"independent_heterogeneous_GP";stopifnot(reject(x))
x<-rows;x$finite_p[[1L]]<-401L;stopifnot(reject(x))
x<-rows;x$finite_p[[1L]]<-400L;stopifnot(reject(x))
x<-rows;x$native_bh_replay_ok[[1L]]<-FALSE;stopifnot(reject(x))
x<-rows;x$same_model_pvalue_check[[1L]]<-FALSE;stopifnot(reject(x))
x<-rows;x$total_native_fit_count[[1L]]<-2L;stopifnot(reject(x))
x<-rows;x$discovery_count[[1L]]<-1L;stopifnot(reject(x))
x<-rows;x$finite_q[[1L]]<-400L;stopifnot(reject(x))
cat("PASS: 192 fixed null cohorts and 768 paired Cook/IF results must be complete; no FDR certification\n")
