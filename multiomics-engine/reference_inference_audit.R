# Independent base-R statistical integrity audit for the *actual* native
# DESeq2 / edgeR-QL / limma-voom reference results in validation CI.
#
# Recompute BH directly from the unadjusted p-values that the R engine
# REPORTED as eligible for adjustment (finite native q values). This handles
# DESeq2 independent filtering explicitly and MUST NOT erase that filtering.
#
# Independent implementation of the Benjamini-Hochberg step-up formula.
# Deliberately does NOT call stats::p.adjust: that function may already be
# used by the reference packages and would not be an independent algorithm.
independent_bh_adjust <- function(p) {
  if (!is.numeric(p) || !length(p) ||
      anyNA(p) || any(!is.finite(p)) || any(p<0 | p>1))
    stop("BH requires a nonempty vector of finite unadjusted probabilities.")
  n <- length(p)
  ord <- order(p,seq_along(p),method="radix")
  sorted <- p[ord]
  ranked <- sorted * (as.numeric(n)/seq_len(n))
  adjusted <- pmin(1,rev(cummin(rev(ranked))))
  result <- numeric(n)
  result[ord] <- adjusted
  names(result) <- names(p)
  result
}


# Important: exact BH replay establishes numerical/provenance integrity, NOT
# calibrated null p-values or nominal FDR control in real studies.
audit_native_r_bh <- function(tab, method, scenario, replicate,
                              feature_col, p_col, q_col,
                              expected_features, tolerance=1e-8) {
  if (!is.data.frame(tab) || !nrow(tab) || !is.character(method) ||
      length(method)!=1L || !nzchar(method) ||
      !is.character(scenario) || length(scenario)!=1L || !nzchar(scenario) ||
      length(replicate)!=1L || !is.finite(replicate) ||
      replicate<1L || replicate!=floor(replicate) ||
      !all(c(feature_col,p_col,q_col) %in% names(tab)) ||
      anyDuplicated(c(feature_col,p_col,q_col)) ||
      !is.character(expected_features) || !length(expected_features) ||
      anyNA(expected_features) || any(!nzchar(expected_features)) ||
      anyDuplicated(expected_features) ||
      length(tolerance)!=1L || !is.finite(tolerance) ||
      tolerance<=0 || tolerance>1e-5)
    stop("Invalid native R BH integrity audit specification.")
  features <- as.character(tab[[feature_col]])
  if (anyNA(features) || any(!nzchar(features)) ||
      anyDuplicated(features) ||
      length(features)!=length(expected_features) ||
      !setequal(features,expected_features))
    stop("Native result features do not match the independent source manifest.")
  p <- tab[[p_col]]
  q <- tab[[q_col]]
  if (!is.numeric(p) || !is.numeric(q) || length(p)!=nrow(tab) ||
      length(q)!=nrow(tab) ||
      any(!is.na(p) & (!is.finite(p) | p<0 | p>1)) ||
      any(!is.na(q) & (!is.finite(q) | q<0 | q>1)))
    stop("Native p-values or BH q-values are not finite probabilities.")
  adjusted <- is.finite(q)
  if (!any(adjusted) || any(adjusted & !is.finite(p)))
    stop("Native finite q-values require at least one finite raw p-value.")
  replay <- independent_bh_adjust(p[adjusted])
  deviation <- max(abs(replay-q[adjusted]))
  if (!is.finite(deviation) || deviation>tolerance)
    stop(sprintf(
      "Native BH integrity mismatch %s/%s/rep%d: max deviation %.8g exceeds %.8g",
      scenario,method,replicate,deviation,tolerance
    ))
  data.frame(
    scenario=scenario,method=method,replicate=as.integer(replicate),
    source_features=length(expected_features),
    finite_raw_p=sum(is.finite(p)),
    finite_adjusted_q=sum(adjusted),
    excluded_raw_p_with_no_q=sum(is.finite(p) & !adjusted),
    raw_p_below_005=sum(is.finite(p) & p<=0.05),
    bh_discoveries_005=sum(adjusted & q<=0.05),
    replay_max_abs_q_difference=deviation,
    native_bh_integrity_verified=TRUE,
    fdr_calibration_certified=FALSE,
    biological_study_certified=FALSE,
    stringsAsFactors=FALSE
  )
}

# Independently replay a KNOWN true-uniform all-null benchmark, then a biased
# anti-conservative p-value negative control. This verifies that the BH
# diagnostic can detect failure; it does NOT validate DESeq2/edgeR/limma models.
run_independent_bh_oracle <- function(replicates=2048L, features=200L,
                                      seed=20261018L, alpha=0.05) {
  if (length(replicates)!=1L || !is.finite(replicates) ||
      replicates<512L || replicates!=floor(replicates) ||
      length(features)!=1L || !is.finite(features) ||
      features<20L || features!=floor(features) ||
      length(seed)!=1L || !is.finite(seed) || seed!=floor(seed) ||
      length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("Invalid prespecified BH reference benchmark.")
  had_seed <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  if(had_seed) prior <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
  on.exit({
    if(had_seed) assign(".Random.seed",prior,.GlobalEnv)
    else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  set.seed(as.integer(seed))
  draws <- vector("list",2L * as.integer(replicates))
  next_index <- 0L
  for(rep in seq_len(as.integer(replicates))) {
    # Distinct seeds/synthetic studies, same fixed generating mechanism.
    for(source in c("uniform_null","anticonservative_beta_half")) {
      p <- if(source=="uniform_null") stats::runif(features)
           else stats::rbeta(features,shape1=0.5,shape2=1)
      q <- independent_bh_adjust(p)
      # Orthogonal implementation comparison: analytic BH step-up formula vs
      # the R package utility. This is algorithmic, not model calibration.
      if (max(abs(q-stats::p.adjust(p,method="BH")))>1e-12)
        stop("Independent BH step-up formula disagrees with stats::p.adjust.")
      next_index <- next_index+1L
      draws[[next_index]] <- data.frame(
        mechanism=source,replicate=as.integer(rep),
        discoveries=sum(q<=alpha),any_false_discovery=as.integer(any(q<=alpha)),
        stringsAsFactors=FALSE)
    }
  }
  rows <- do.call(rbind,draws)
  rownames(rows) <- NULL
  summary <- do.call(rbind,lapply(split(rows,rows$mechanism),function(x) {
    n <- nrow(x)
    k <- sum(x$any_false_discovery)
    ci <- stats::binom.test(k,n,conf.level=.95)$conf.int
    data.frame(mechanism=x$mechanism[[1L]],studies=n,
      datasets_with_any_false_discovery=k,
      empirical_global_null_fdr=k/n,exact_mc_95_low=ci[[1L]],
      exact_mc_95_high=ci[[2L]],
      mean_bh_discoveries=mean(x$discoveries),
      prespecified_bh_alpha=alpha,
      native_method_validated=FALSE,biological_study_certified=FALSE)
  }))
  rownames(summary) <- NULL
  baseline <- summary[summary$mechanism=="uniform_null",]
  negative <- summary[summary$mechanism=="anticonservative_beta_half",]
  # Deliberately generous regression envelopes; NEVER a test proving exactly
  # 5% FDR in DESeq2 or any user experiment.
  if(nrow(baseline)!=1L || nrow(negative)!=1L ||
     abs(baseline$empirical_global_null_fdr-alpha)>0.03 ||
     negative$empirical_global_null_fdr<=0.20)
    stop("BH uniform-null oracle or anti-conservative detection failed.")
  list(per_study=rows,summary=summary,
    generating_mechanisms=c("independent Uniform(0,1) null p-values",
                            "independent Beta(0.5,1) anti-conservative negative control"),
    scope="BH algorithmic integrity only; not native-method model calibration")
}
