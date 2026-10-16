#!/usr/bin/env Rscript
# Test metrics and experimental truth without loading any Bioconductor package.
source("multiomics-engine/mixed_truth_reference_metrics.R")
truth <- setNames(c(0,0,1.25,-1.25,0,0),paste0("G",seq_len(6L)))
tab <- data.frame(feature=names(truth),q=c(.8,.02,.001,.01,NA,.9),
                  log2fc=c(.1,.2,1.1,-1.3,0,.4))
r <- mixed_truth_feature_metrics(tab,truth,"DESeq2","five_pct",1L,123L)
stopifnot(nrow(r)==1L,r$total_planted_features==6L,
  r$tested_features==6L,r$finite_q==5L,r$n_true==2L,
  r$n_null==4L,r$selected_features==3L,
  r$false_discoveries==1L,r$true_discoveries==2L,
  isTRUE(all.equal(r$false_discovery_proportion,1/3)),
  r$power==1,r$wrong_direction_true_discoveries==0L,
  !("scientific_certified" %in% names(r)))
wrong <- tab
wrong$log2fc[3L] <- -1.5
bad_sign <- mixed_truth_feature_metrics(wrong,truth,"DESeq2","five_pct",1L,123L)
stopifnot(bad_sign$wrong_direction_true_discoveries==1L)
# Pure null significant discoveries, with no reported alternatives in output.
missing <- tab[-4L,]
miss <- mixed_truth_feature_metrics(missing,truth,"DESeq2","five_pct",1L,123L)
stopifnot(miss$power==.5,miss$n_true==2L,
  miss$tested_features==5L)
reject <- function(t=tab,known=truth) inherits(
  try(mixed_truth_feature_metrics(t,known,"DESeq2","5%",1L,123L),silent=TRUE),
  "try-error")
bad <- tab;bad$feature[2L] <- bad$feature[1L]
stopifnot(reject(t=bad))
bad <- tab;bad$q[1L] <- -0.01
stopifnot(reject(t=bad))
bad <- tab;bad$q[1L] <- Inf
stopifnot(reject(t=bad))
bad <- tab;bad$log2fc[2L] <- Inf
stopifnot(reject(t=bad))
bad <- tab;bad$feature[2L] <- "EXTRA"
stopifnot(reject(t=bad))
bad_truth <- truth; names(bad_truth)[2L] <- names(bad_truth)[1L]
stopifnot(reject(known=bad_truth))

simulated_rows <- list()
for(scenario in c("one_pct","five_pct","twenty_pct")) {
  for(method in c("DESeq2","edgeR_QL","limma_voom")) {
    for(i in seq_len(20L)) {
      x <- r
      x$scenario <- scenario
      x$method <- method
      x$replicate <- i
      x$seed <- i+match(scenario,c("one_pct","five_pct","twenty_pct"))*1000L
      if(i>10L) {
        x$false_discoveries <- 0L
        x$false_discovery_proportion <- 0
      }
      simulated_rows[[length(simulated_rows)+1L]] <- x
    }
  }
}
simulated <- do.call(rbind,simulated_rows)
set.seed(9876L)
seed_before <- get(".Random.seed",envir=.GlobalEnv)
out <- mixed_truth_study_summary(simulated)
# A bootstrap percentile interval collapses to [0,0] for all-zero FDP.
# A finite-sample bound MUST remain strictly positive for n=8, 12 or 120.
for(n in c(8L,12L,120L)) {
  zero <- rep(0,n)
  conservative <- rnaseq_fdr_hoeffding_bounds(zero,cells=9L)
  stopifnot(conservative$lower==0,
    conservative$upper>0,
    conservative$radius>0,
    conservative$upper<=1,
    conservative$simultaneous_cells==9L,
    identical(conservative$scope,"iid_simulated_studies_only"),
    !conservative$biological_study_certified)
}
eight <- rnaseq_fdr_hoeffding_bounds(rep(0,8L),cells=9L)
hundred_twenty <- rnaseq_fdr_hoeffding_bounds(rep(0,120L),cells=9L)
stopifnot(hundred_twenty$upper<eight$upper)
unadjusted <- rnaseq_fdr_hoeffding_bounds(rep(0,12L),cells=1L)
adjusted <- rnaseq_fdr_hoeffding_bounds(rep(0,12L),cells=9L)
stopifnot(adjusted$upper>=unadjusted$upper,
  identical(adjusted$method,"two_sided_hoeffding_bonferroni"))
for(invalid in list(c(0,NA),c(-.01,.1),c(.1,1.1),c(Inf,0),0)) {
  stopifnot(inherits(try(rnaseq_fdr_hoeffding_bounds(invalid),
                         silent=TRUE),"try-error"))
}
stopifnot(inherits(try(rnaseq_fdr_hoeffding_bounds(c(0,0),cells=0L),
                       silent=TRUE),"try-error"))

stopifnot(identical(seed_before,get(".Random.seed",envir=.GlobalEnv)),
  nrow(out)==9L,all(out$independent_datasets==20L),
  all(abs(out$estimated_fdr-1/6)<1e-12),
  all(out$bootstrap_fdr_95_low<=out$estimated_fdr),
  all(out$bootstrap_fdr_95_high>=out$estimated_fdr),
  all(out$mean_power==1),all(!out$scientific_certified),
  all(out$conservative_fdr_95_high>=out$estimated_fdr),
  all(out$conservative_fdr_95_low<=out$estimated_fdr),
  all(out$conservative_fdr_family_cells==9L),
  all(out$conservative_fdr_95_radius>0),
  all(out$conservative_fdr_bound_method=="two_sided_hoeffding_bonferroni"),
  all(!out$scenario_calibration_upper_below_5pct))
bad <- simulated
bad$seed[2L] <- bad$seed[1L]
stopifnot(inherits(try(mixed_truth_study_summary(bad),silent=TRUE),"try-error"))
bad <- simulated
bad$replicate[2L] <- bad$replicate[1L]
stopifnot(inherits(try(mixed_truth_study_summary(bad),silent=TRUE),"try-error"))

# Check that a downstream summary cannot silently report a zero upper
# uncertainty simply because zero false discoveries occurred in 8 studies.
only_zero <- simulated[simulated$scenario=="one_pct" &
  simulated$method=="DESeq2" & simulated$replicate<=8L,,drop=FALSE]
only_zero$false_discovery_proportion <- 0
only_zero$false_discoveries <- 0L
zero_out <- mixed_truth_study_summary(only_zero,min_replicates=8L,
                                      bootstrap_repeats=200L)
stopifnot(zero_out$estimated_fdr==0,
  zero_out$bootstrap_fdr_95_high==0,
  zero_out$conservative_fdr_95_high>.05,
  !zero_out$scenario_calibration_upper_below_5pct)
cat("RNA-seq mixed truth FDR, power and conservative finite-sample bounds: PASS (no biological certification)\n")
