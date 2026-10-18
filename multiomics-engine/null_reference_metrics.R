# Pure base-R audit helpers for the global-null reference benchmark.
# FDR = FWER under an all-null experiment because any selected feature is false.
# Intervals measure Monte Carlo uncertainty only. The simultaneous interval
# Bonferroni-corrects across all scenario × method cells from the same run;
# models fitted on shared data may be dependent, which union control permits.
null_reference_summary <- function(records, min_replicates=8L, confidence=0.95,
                                   simultaneous_family_alpha=0.05) {
  required <- c("scenario","replicate","method","seed","features_tested",
                "finite_q","discoveries_q05","false_discovery_fraction")
  if (!is.data.frame(records) || !nrow(records) ||
      !all(required %in% names(records)))
    stop("Null-calibration records must be a nonempty data frame with all required fields.")
  if (length(min_replicates)!=1L || !is.finite(min_replicates) ||
      min_replicates<2L || min_replicates!=as.integer(min_replicates))
    stop("Invalid minimum number of independent simulated datasets.")
  if (length(confidence)!=1L || !is.finite(confidence) ||
      confidence<=0 || confidence>=1)
    stop("Confidence level must lie strictly between zero and one.")
  if(length(simultaneous_family_alpha)!=1L ||
     !is.numeric(simultaneous_family_alpha) ||
     !is.finite(simultaneous_family_alpha) ||
     simultaneous_family_alpha<=0 || simultaneous_family_alpha>=1)
    stop("Simultaneous Monte Carlo family alpha must be in (0,1).")
  if (anyNA(records[, required]) ||
      any(!nzchar(as.character(records$scenario))) ||
      any(!nzchar(as.character(records$method))) ||
      any(!is.finite(records$replicate)) ||
      any(!is.finite(records$seed)) ||
      any(!is.finite(records$features_tested)) ||
      any(!is.finite(records$finite_q)) ||
      any(!is.finite(records$discoveries_q05)) ||
      any(!is.finite(records$false_discovery_fraction)) ||
      any(records$replicate<1 | records$replicate!=floor(records$replicate)) ||
      any(records$features_tested<1 | records$features_tested!=floor(records$features_tested)) ||
      any(records$finite_q<0 | records$finite_q>records$features_tested) ||
      any(records$finite_q!=floor(records$finite_q)) ||
      any(records$discoveries_q05<0 |
          records$discoveries_q05>records$finite_q |
          records$discoveries_q05!=floor(records$discoveries_q05)) ||
      any(!records$false_discovery_fraction %in% c(0,1)) ||
      any((records$discoveries_q05>0)!=as.logical(records$false_discovery_fraction)))
    stop("Invalid or inconsistent null-discovery records.")
  split_groups <- split(records, list(records$scenario,records$method), drop=TRUE)
  family_cells <- length(split_groups)
  if (family_cells < 1L) stop("No statistical comparison cells.")
  simultaneous_cell_confidence <- 1 - simultaneous_family_alpha/family_cells
  out <- lapply(split_groups, function(cell) {
    n <- nrow(cell)
    if (n < min_replicates || anyDuplicated(cell$replicate) ||
        anyDuplicated(cell$seed))
      stop("Each scenario-method cell needs unique independent datasets and enough replicates.")
    k <- sum(cell$false_discovery_fraction)
    interval <- stats::binom.test(k,n,conf.level=confidence)$conf.int
    simultaneous <- stats::binom.test(k,n,
      conf.level=simultaneous_cell_confidence)$conf.int
    data.frame(scenario=as.character(cell$scenario[[1L]]),
      method=as.character(cell$method[[1L]]),
      replicates=n, datasets_with_any_false_discovery=k,
      empirical_global_null_fdr=k/n,
      mc_exact_ci_low=unname(interval[[1L]]),
      mc_exact_ci_high=unname(interval[[2L]]),
      mc_confidence=confidence,
      simultaneous_ci_low=unname(simultaneous[[1L]]),
      simultaneous_ci_high=unname(simultaneous[[2L]]),
      simultaneous_family_cells=family_cells,
      simultaneous_family_confidence=1-simultaneous_family_alpha,
      simultaneous_interval_method="Clopper-Pearson Bonferroni",
      nominal_bh_q=0.05,
      simultaneous_inflation_detected=(simultaneous[[1L]]>0.05),
      simultaneous_inconclusive=(simultaneous[[1L]]<=0.05 &&
                                  simultaneous[[2L]]>=0.05),
      scientific_certified=FALSE,
      mean_bh_discoveries=mean(cell$discoveries_q05),
      mean_finite_q=mean(cell$finite_q),
      mean_features_tested=mean(cell$features_tested),
      stringsAsFactors=FALSE)
  })
  result <- do.call(rbind,out)
  rownames(result) <- NULL
  result[order(result$scenario,result$method),,drop=FALSE]
}

# A CI containing 5% cannot demonstrate calibration; even 0/12 gives a very
# wide upper bound. Only gross empirical inflation is a failure gate.
null_reference_gross_failures <- function(summary, threshold=0.35) {
  if (!is.data.frame(summary) ||
      !all(c("scenario","method","empirical_global_null_fdr") %in% names(summary)) ||
      !is.finite(threshold) || threshold<0 || threshold>1)
    stop("Invalid global-null failure gate.")
  summary[summary$empirical_global_null_fdr>threshold,,drop=FALSE]
}
