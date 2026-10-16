# Pure base R metrics. Independent simulated STUDIES, not genes, are the
# resampling unit for the mean false-discovery proportion.
mixed_truth_feature_metrics <- function(result, planted_log2_fc, method, scenario,
                                        replicate, seed, alpha=0.05) {
  required <- c("feature","q","log2fc")
  if (!is.data.frame(result) || !all(required %in% names(result)) ||
      !nrow(result) || is.null(names(planted_log2_fc)) ||
      anyDuplicated(names(planted_log2_fc)) ||
      anyDuplicated(result$feature) ||
      anyNA(result$feature) || any(!nzchar(result$feature)) ||
      anyNA(planted_log2_fc) ||
      any(!is.finite(planted_log2_fc)) ||
      !length(planted_log2_fc) ||
      any(!result$feature %in% names(planted_log2_fc)) ||
      length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("Invalid mixed-truth feature table or prespecified ground truth.")
  if (any(!is.na(result$q) & (!is.finite(result$q) | result$q<0 | result$q>1)))
    stop("Adjusted p-values must lie within [0,1], or be missing.")
  if (any(!is.na(result$log2fc) & !is.finite(result$log2fc)))
    stop("Invalid log2 fold-change estimates.")
  truth <- unname(planted_log2_fc[result$feature])
  selected <- !is.na(result$q) & result$q<=alpha
  actual <- truth!=0
  false_discoveries <- sum(selected & !actual)
  true_discoveries <- sum(selected & actual)
  selections <- sum(selected)
  n_true <- sum(planted_log2_fc!=0)
  observed_signal <- selected & actual & is.finite(result$log2fc)
  wrong_sign <- selected & actual & is.finite(result$log2fc) &
                (sign(result$log2fc)!=sign(truth))
  # Wrong-direction alternative calls are *not* classical null discoveries;
  # report them separately rather than altering standard FDR.
  data.frame(
    scenario=as.character(scenario),method=as.character(method),
    replicate=as.integer(replicate),seed=as.integer(seed),
    total_planted_features=length(planted_log2_fc),
    tested_features=nrow(result),
    finite_q=sum(is.finite(result$q)),n_true=n_true,
    n_null=length(planted_log2_fc)-n_true,
    selected_features=selections,
    true_discoveries=true_discoveries,
    false_discoveries=false_discoveries,
    false_discovery_proportion=if(selections) false_discoveries/selections else 0,
    power=true_discoveries/n_true,
    wrong_direction_true_discoveries=sum(wrong_sign),
    wrong_direction_fraction=if(sum(observed_signal))
      sum(wrong_sign)/sum(observed_signal) else NA_real_,
    mean_absolute_log2fc_error=if(any(actual & is.finite(result$log2fc)))
      mean(abs(result$log2fc[actual & is.finite(result$log2fc)]-
               truth[actual & is.finite(result$log2fc)])) else NA_real_,
    signed_log2fc_bias=if(any(actual & is.finite(result$log2fc)))
      mean(result$log2fc[actual & is.finite(result$log2fc)]-
             truth[actual & is.finite(result$log2fc)]) else NA_real_,
    stringsAsFactors=FALSE)
}

# Nonasymptotic simultaneous Hoeffding bounds for the MEAN FDP across
# independent and identically distributed synthetic studies. Every study's
# FDP lies in [0,1] and is one unit of replication. Bonferroni union control
# covers arbitrarily dependent methods applied to the same simulated studies.
# These are bounds for a FIXED synthetic-generating mechanism only.
# They are intentionally conservative; zero observed FDP NEVER gives zero
# upper uncertainty with a small number of simulated studies.
rnaseq_fdr_hoeffding_bounds <- function(fdp, cells=1L,
                                       family_alpha=.05) {
  if (!is.numeric(fdp) || length(fdp)<2L ||
      anyNA(fdp) || any(!is.finite(fdp)) ||
      any(fdp<0 | fdp>1) ||
      length(cells)!=1L || !is.finite(cells) ||
      cells<1L || cells!=floor(cells) ||
      length(family_alpha)!=1L || !is.finite(family_alpha) ||
      family_alpha<=0 || family_alpha>=1)
    stop("FDR bounds require >=2 independent bounded study FDPs and valid simultaneous family.")
  epsilon <- sqrt(log(2*cells/family_alpha)/(2*length(fdp)))
  center <- mean(fdp)
  list(
    lower=max(0,center-epsilon),
    upper=min(1,center+epsilon),
    radius=epsilon,
    independent_studies=length(fdp),
    simultaneous_cells=as.integer(cells),
    confidence=1-family_alpha,
    method="two_sided_hoeffding_bonferroni",
    scope="iid_simulated_studies_only",
    # Never conflate calibration of a simulation with validation of a
    # particular uploaded biological cohort.
    biological_study_certified=FALSE)
}

mixed_truth_study_summary <- function(records, min_replicates=8L,
                                      bootstrap_repeats=1000L,
                                      seed=20261011L) {
  cols <- c("scenario","method","replicate","seed",
            "false_discovery_proportion","power",
            "selected_features","false_discoveries","tested_features",
            "total_planted_features","n_true","mean_absolute_log2fc_error",
            "signed_log2fc_bias")
  if (!is.data.frame(records) || !nrow(records) ||
      !all(cols %in% names(records)) ||
      anyNA(records[,cols]) ||
      any(!nzchar(records$scenario)) || any(!nzchar(records$method)) ||
      any(!is.finite(records$false_discovery_proportion)) ||
      any(records$false_discovery_proportion<0 |
          records$false_discovery_proportion>1) ||
      any(!is.finite(records$power)) || any(records$power<0 | records$power>1) ||
      any(records$n_true<1L) || any(records$tested_features<1L) ||
      any(records$tested_features>records$total_planted_features) ||
      any(records$false_discoveries>records$selected_features) ||
      any(records$selected_features>records$tested_features) ||
      any(records$replicate!=floor(records$replicate)) ||
      !is.finite(min_replicates) || min_replicates<2L ||
      !is.finite(bootstrap_repeats) || bootstrap_repeats<100L)
    stop("Invalid mixed-truth per-study metrics.")
  # Stable local RNG and restore caller's state: CI is reproducible without
  # contaminating the next simulated biological dataset.
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) original_seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(had_seed) assign(".Random.seed",original_seed,envir=.GlobalEnv)
          else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
            rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(seed)
  cells <- split(records,list(records$scenario,records$method),drop=TRUE)
  comparison_count <- length(cells)
  summaries <- lapply(cells,function(cell) {
    n <- nrow(cell)
    if(n<min_replicates || anyDuplicated(cell$replicate) ||
       anyDuplicated(cell$seed))
      stop("Insufficient or non-independent simulated datasets.")
    fdp <- cell$false_discovery_proportion
    bound <- rnaseq_fdr_hoeffding_bounds(fdp,cells=comparison_count)
    draws <- replicate(as.integer(bootstrap_repeats),
      mean(sample(fdp,size=n,replace=TRUE)))
    ci <- unname(stats::quantile(draws,probs=c(.025,.975),names=FALSE))
    data.frame(scenario=cell$scenario[[1L]],method=cell$method[[1L]],
      independent_datasets=n,estimated_fdr=mean(fdp),
      bootstrap_fdr_95_low=ci[[1L]],
      bootstrap_fdr_95_high=ci[[2L]],
      conservative_fdr_95_low=bound$lower,
      conservative_fdr_95_high=bound$upper,
      conservative_fdr_95_radius=bound$radius,
      conservative_fdr_family_cells=bound$simultaneous_cells,
      conservative_fdr_bound_method=bound$method,
      conservative_fdr_scope=bound$scope,
      scenario_calibration_upper_below_5pct=(bound$upper<=.05),
      mean_power=mean(cell$power),
      mean_false_discoveries=mean(cell$false_discoveries),
      mean_selected_features=mean(cell$selected_features),
      mean_tested_features=mean(cell$tested_features),
      mean_absolute_log2fc_error=mean(cell$mean_absolute_log2fc_error),
      mean_signed_log2fc_bias=mean(cell$signed_log2fc_bias),
      nominal_bh_q=0.05,
      scientific_certified=FALSE,
      stringsAsFactors=FALSE)
  })
  out <- do.call(rbind,summaries)
  rownames(out) <- NULL
  out[order(out$scenario,out$method),,drop=FALSE]
}
