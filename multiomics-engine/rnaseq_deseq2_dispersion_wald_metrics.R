# Fixed-horizon whole-study summary of the true-dispersion DESeq2 diagnostic.
# This is an exploratory mechanism comparison under a *synthetic* global null,
# NOT a proof that real-data nominal BH FDR is controlled.
rnaseq_dispersion_wald_summary <- function(rows,replicates=32L,
                                            seed_offset=3000L) {
  methods <- c("native_default","oracle_true_dispersion",
               "oracle_true_depth_and_dispersion")
  scenarios <- c("balanced_base","balanced_depth_dispersion",
                 "unequal_partial_batch")
  needed <- c("scenario","method","replicate","seed",
              "fdr_any_discovery","bh_discoveries","finite_p",
              "finite_q","wald_stat_squared_mean","dispersion_log_rmse",
              "depth_factor_log_rmse","native_bh_integrity_verified")
  if(!is.data.frame(rows)||!nrow(rows)||
     !all(needed %in% names(rows))||
     anyNA(rows[,needed])||
     !identical(sort(unique(as.character(rows$scenario))),sort(scenarios))||
     !identical(sort(unique(as.character(rows$method))),sort(methods))||
     length(replicates)!=1L||!is.finite(replicates)||
     replicates!=floor(replicates)||replicates<16L||
     length(seed_offset)!=1L||!is.finite(seed_offset)||
     seed_offset<1000L||
     any(!rows$fdr_any_discovery%in%c(0L,1L))||
     any(rows$bh_discoveries<0L)||
     any((rows$bh_discoveries>0L)!=(rows$fdr_any_discovery==1L))||
     any(!rows$native_bh_integrity_verified)||
     any(!is.finite(rows$dispersion_log_rmse))||
     any(!is.finite(rows$depth_factor_log_rmse))||
     any(!is.finite(rows$wald_stat_squared_mean))||
     any(rows$dispersion_log_rmse<0|rows$depth_factor_log_rmse<0)||
     any(rows$finite_p<300L|rows$finite_q<300L))
    stop("Invalid or incomplete native DESeq2 dispersion-Wald evidence.")
  cells <- split(rows,list(rows$scenario,rows$method),drop=TRUE)
  if(length(cells)!=9L||
     any(vapply(cells,nrow,integer(1))!=replicates))
    stop("Missing predeclared global-null scenario/method cells.")
  records <- list()
  for(s in scenarios) {
    local <- rows[rows$scenario==s,,drop=FALSE]
    blocks <- lapply(methods,function(m) {
      dat <- local[local$method==m,,drop=FALSE]
      dat <- dat[order(dat$replicate),,drop=FALSE]
      if(nrow(dat)!=replicates||
         anyDuplicated(dat$replicate)>0L||
         anyDuplicated(dat$seed)>0L||
         !identical(as.integer(dat$replicate),
                    seq.int(as.integer(seed_offset+1L),
                            as.integer(seed_offset+replicates))))
        stop("Replicates duplicated, missing or recycled from previous benchmarks.")
      dat
    })
    if(!identical(blocks[[1L]]$seed,blocks[[2L]]$seed)||
       !identical(blocks[[1L]]$seed,blocks[[3L]]$seed))
      stop("Three normalization/dispersion arms do not use EXACTLY the same cohort.")
    for(j in seq_along(methods)) {
      data <- blocks[[j]]
      n <- nrow(data)
      events <- sum(data$fdr_any_discovery)
      ci <- stats::binom.test(events,n,conf.level=.95)$conf.int
      records[[length(records)+1L]]<-data.frame(
        scenario=s,method=methods[[j]],
        independent_global_null_datasets=n,
        datasets_with_any_false_discovery=events,
        empirical_global_null_fdr=events/n,
        exact_binomial_ci95_low=ci[[1L]],
        exact_binomial_ci95_high=ci[[2L]],
        mean_number_of_bh_false_calls=mean(data$bh_discoveries),
        mean_native_null_wald_stat_squared=mean(data$wald_stat_squared_mean),
        mean_dispersion_log_rmse=mean(data$dispersion_log_rmse),
        mean_depth_factor_log_rmse=mean(data$depth_factor_log_rmse),
        both_arms_computed_on_same_source=TRUE,
        true_dispersion_available_in_human_samples=FALSE,
        universal_FDR_control_certified=FALSE,
        mechanism_causally_identified=FALSE,
        stringsAsFactors=FALSE)
    }
  }
  ans<-do.call(rbind,records)
  rownames(ans)<-NULL
  ans
}
