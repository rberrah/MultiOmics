# Pre-registered complete study-level summary for native DESeq2 shrinkage ablation.
# Four arms use the SAME count matrix, size factors and Wald test.
# Gene-wise estimates are Cox-Reid adjusted; not an independent assay.
rnaseq_deseq2_dispersion_component_summary <- function(dat,
    per_scenario=32L, replicate_start=4001L) {
  scenarios <- c("balanced_base","balanced_depth_dispersion",
                 "unequal_partial_batch")
  arms <- c("map_final","gene_wise_CR","fitted_trend","true_nb_alpha")
  required <- c("scenario","replicate","seed","arm",
     "source_features","modeled_features","finite_p","finite_q",
     "false_calls_bh05","any_false_call","native_BH_replay_ok",
     "dispersion_log_rmse","mean_wald_stat_squared",
     "raw_p_under_005","raw_p_under_001")
  if(!is.data.frame(dat)||!nrow(dat)||
     !all(required %in% names(dat))||
     anyNA(dat[,required])||
     length(per_scenario)!=1L||!is.finite(per_scenario)||
     per_scenario!=floor(per_scenario)||per_scenario<16L||
     length(replicate_start)!=1L||!is.finite(replicate_start)||
     replicate_start<3001L||
     !setequal(unique(dat$scenario),scenarios)||
     !setequal(unique(dat$arm),arms)||
     any(dat$source_features!=400L)||
     any(dat$modeled_features<300L|dat$modeled_features>400L)||
     any(dat$finite_p<300L|dat$finite_q<300L)||
     any(dat$false_calls_bh05<0L)||
     any(!dat$any_false_call%in%c(0L,1L))||
     any((dat$false_calls_bh05>0L)!=(dat$any_false_call==1L))||
     any(!dat$native_BH_replay_ok)||
     any(!is.finite(dat$dispersion_log_rmse))||
     any(dat$dispersion_log_rmse<0)||
     any(!is.finite(dat$mean_wald_stat_squared))||
     any(dat$mean_wald_stat_squared<0)||
     any(dat$raw_p_under_005<0|dat$raw_p_under_001<0)||
     any(dat$raw_p_under_001>dat$raw_p_under_005))
     stop("Invalid / missing native dispersion component ablation evidence.")
  group <- split(dat,list(dat$scenario,dat$arm),drop=TRUE)
  if(length(group)!=12L||
     any(vapply(group,nrow,integer(1))!=per_scenario))
    stop("Expected exactly twelve pre-registered scenario-by-arm cells.")
  result <- list()
  for(sc in scenarios) {
    byarm <- lapply(arms,function(a) {
      x<-dat[dat$scenario==sc & dat$arm==a,,drop=FALSE]
      x<-x[order(x$replicate),,drop=FALSE]
      if(nrow(x)!=per_scenario||
         anyDuplicated(x$replicate)>0L||
         anyDuplicated(x$seed)>0L||
         !identical(as.integer(x$replicate),
           as.integer(seq.int(replicate_start,
                              length.out=per_scenario))))
        stop("Missing, duplicated or unregistered synthetic cohort.")
      x
    })
    ref<-byarm[[1L]]
    for(j in seq_along(arms)){
      d<-byarm[[j]]
      if(!identical(d$seed,ref$seed)||
         !identical(d$source_features,ref$source_features)||
         !identical(d$modeled_features,ref$modeled_features))
        stop("Dispersion models not fitted on same cohort and gene family.")
      n<-nrow(d)
      k<-sum(d$any_false_call)
      ci<-stats::binom.test(k,n,conf.level=.95)$conf.int
      result[[length(result)+1L]]<-data.frame(
        scenario=sc,arm=arms[[j]],
        whole_null_studies=n,
        studies_with_bh_false_calls=k,
        empirical_global_null_FDR=k/n,
        individual_clopper_pearson_low=ci[[1L]],
        individual_clopper_pearson_high=ci[[2L]],
        mean_false_calls_per_study=mean(d$false_calls_bh05),
        mean_dispersion_log_rmse=mean(d$dispersion_log_rmse),
        mean_wald_stat_squared=mean(d$mean_wald_stat_squared),
        mean_raw_p_fraction_under_005=mean(d$raw_p_under_005/d$finite_p),
        mean_raw_p_fraction_under_001=mean(d$raw_p_under_001/d$finite_p),
        no_source_gene_pseudoreplication=TRUE,
        no_post_hoc_optimized_method=TRUE,
        universal_FDR_calibration_certified=FALSE,
        source_params_known_in_real_biology=FALSE,
        stringsAsFactors=FALSE)
    }
  }
  ans<-do.call(rbind,result)
  rownames(ans)<-NULL
  ans
}
