# Whole-study paired diagnostic of DESeq2's default Cook's distance screening
# and independently chosen BH gene-filter. Nothing here changes production.
rnaseq_cook_if_summary <- function(rows, repetitions=32L) {
  original <- c("balanced_base","balanced_depth_dispersion",
                "unequal_partial_batch")
  independent <- c("gp_balanced_heterogeneous","gp_depth_heterogeneous",
                   "gp_unbalanced_batch")
  scenarios <- c(original,independent)
  arms <- c("cooks_on_if_on","cooks_on_if_off",
            "cooks_off_if_on","cooks_off_if_off")
  required <- c("scenario","source_generator","replicate","seed","arm",
       "source_genes","modeled_genes","total_native_fit_count",
       "finite_p","finite_q","discovery_count","any_discovery",
       "p_removed_by_cooks","q_removed_by_filter",
       "raw_p_005","native_bh_replay_ok","same_model_pvalue_check")
  if(!is.data.frame(rows)||!nrow(rows)||!all(required %in% names(rows))||
     anyNA(rows[,required])||
     !setequal(unique(as.character(rows$scenario)),scenarios)||
     !setequal(unique(as.character(rows$arm)),arms)||
     any(!rows$source_generator%in%
             c("original_smooth_NB","independent_heterogeneous_GP"))||
     any(rows$source_genes!=400L)||
     any(rows$modeled_genes<300L|rows$modeled_genes>400L)||
     any(rows$total_native_fit_count!=1L)||
     any(rows$finite_p<0L|rows$finite_q<0L)||
     any(rows$finite_p>rows$modeled_genes|
         rows$finite_q>rows$finite_p)||
     any(rows$discovery_count<0L|rows$discovery_count>rows$finite_q)||
     any(!rows$any_discovery%in%c(0L,1L))||
     any((rows$discovery_count>0L)!=(rows$any_discovery==1L))||
     any(rows$p_removed_by_cooks<0L|
         rows$p_removed_by_cooks>rows$modeled_genes)||
     any(rows$q_removed_by_filter<0L|
         rows$q_removed_by_filter>rows$modeled_genes)||
     any(rows$raw_p_005<0L|rows$raw_p_005>rows$finite_p)||
     any(!rows$native_bh_replay_ok)||
     any(!rows$same_model_pvalue_check))
    stop("Incomplete or forged native Cook/filter evidence.")
  if(length(repetitions)!=1L||!is.finite(repetitions)||
     repetitions!=floor(repetitions)||repetitions<16L)
    stop("Invalid fixed experiment size.")
  cells<-split(rows,list(rows$scenario,rows$arm),drop=TRUE)
  if(length(cells)!=24L||
     any(vapply(cells,nrow,integer(1))!=repetitions))
    stop("All twenty-four source×scenario×screen cells are mandatory.")
  summaries<-list()
  for(sc in scenarios) {
    origin<-if(sc%in%original)"original_smooth_NB" else
      "independent_heterogeneous_GP"
    first<-if(sc%in%original)6001L else 5001L
    k<-match(sc,if(origin=="original_smooth_NB")original else independent)
    paired<-lapply(arms,function(arm) {
      x<-rows[rows$scenario==sc & rows$arm==arm,,drop=FALSE]
      x<-x[order(x$replicate),,drop=FALSE]
      if(nrow(x)!=repetitions||anyDuplicated(x$replicate)||
         !identical(as.integer(x$replicate),
            as.integer(seq.int(first,length.out=repetitions)))||
         !identical(as.character(x$source_generator),
            rep(origin,repetitions)))
         stop("Study generator labels or fresh original/Gamma–Poisson seeds invalid.")
      # Generator seed contracts intentionally not interchangeable:
      # original rnbinom seed 780000 + scenario*10000 + replicate;
      # independent Gamma–Poisson 910000 + scenario*10000 + replicate.
      seed_start<-if(origin=="original_smooth_NB")780000L else 910000L
      expected<-seed_start+as.integer(k*10000L)+as.integer(x$replicate)
      if(!identical(as.integer(x$seed),as.integer(expected)))
        stop("A scenario was relabeled or an original per-cohort seed was recycled.")
      x
    })
    for(x in paired) {
      if(!identical(x$seed,paired[[1L]]$seed)||
         !identical(x$modeled_genes,paired[[1L]]$modeled_genes))
        stop("Four screens did not run on identical native DESeq2 model fits.")
    }
    # Cook screening cannot magically restore p-values. Independent
    # filtering only affects q values (p values identical at same Cook flag).
    for(j in c(1L,3L)) {
      x<-paired[[j]]
      y<-paired[[j+1L]]
      if(!identical(x$finite_p,y$finite_p)||
         !identical(x$p_removed_by_cooks,y$p_removed_by_cooks)||
         any(x$finite_q>y$finite_q))
        stop("Independent filtering incorrectly changed raw Wald p-value eligibility.")
    }
    for(j in c(1L,2L)) {
      cook<-paired[[j]]
      no_cook<-paired[[j+2L]]
      if(any(cook$finite_p>no_cook$finite_p)||
         any(cook$p_removed_by_cooks<
             no_cook$p_removed_by_cooks))
        stop("Cook's filtering increased eligible p values or misreported exclusions.")
    }
    for(j in seq_along(arms)) {
      x<-paired[[j]]
      event<-sum(x$any_discovery)
      ci<-stats::binom.test(event,repetitions)$conf.int
      summaries[[length(summaries)+1L]]<-data.frame(
        scenario=sc,source_generator=origin,
        arm=arms[[j]],whole_null_studies=repetitions,
        studies_with_any_false_discovery=event,
        empirical_global_null_FDR=event/repetitions,
        individual_exact_ci95_low=ci[[1L]],
        individual_exact_ci95_high=ci[[2L]],
        mean_discoveries=mean(x$discovery_count),
        mean_native_eligible_p=mean(x$finite_p),
        mean_native_eligible_q=mean(x$finite_q),
        mean_raw_p_under_005=mean(x$raw_p_005),
        mean_genes_p_removed_by_cooks=mean(x$p_removed_by_cooks),
        mean_genes_q_removed_by_IF=mean(x$q_removed_by_filter),
        model_refitted_per_screen=FALSE,
        production_DESeq2_default_changed=FALSE,
        FDR_5pct_certified=FALSE,
        stringsAsFactors=FALSE)
    }
  }
  result<-do.call(rbind,summaries)
  rownames(result)<-NULL
  result
}
