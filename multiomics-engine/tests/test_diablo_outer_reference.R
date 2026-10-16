#!/usr/bin/env Rscript
# Actual mixOmics DIABLO on independently held-out subjects.
# NOT a substitute for external clinical validation or nested selection among
# many candidate multiblock methods. Synthetic known effects/null controls.
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/validation_benchmarks.R")
if(!requireNamespace("mixOmics",quietly=TRUE))stop("mixOmics must be installed for this test.")
set.seed(20261009)
n <- 64L
ids <- sprintf("P%03d",seq_len(n))
y <- setNames(rep(c("control","treated"),each=32L),ids)
a <- matrix(stats::rnorm(n*14),nrow=n,ncol=14,
  dimnames=list(ids,paste0("gene",seq_len(14))))
b <- matrix(stats::rnorm(n*12),nrow=n,ncol=12,
  dimnames=list(ids,paste0("protein",seq_len(12))))
# Strong signal is deliberate, because this verifies the real package
# code path rather than making a claim about weak-effect power.
a[y=="treated",1:4] <- a[y=="treated",1:4]+3.0
b[y=="treated",1:4] <- b[y=="treated",1:4]+3.0
blocks <- list(rna=a,proteomics=b)
out <- tempfile("diablo_outer_")
result <- run_diablo_outer_holdout(
  blocks,y,out,seeds=c(20261009L,20261010L),tune=TRUE
)
stopifnot(
  identical(result$status,"ok"),
  length(result$summary$folds)==2L,
  all(is.finite(c(result$summary$meanBER,result$summary$maxBER))),
  result$summary$meanBER>=0 && result$summary$meanBER<=1,
  result$summary$meanBER<=0.40,
  all(c("sampleId","actual","predicted","correct","seed")%in%names(result$predictions)),
  length(unique(result$predictions$seed))==2L,
  file.exists(file.path(out,"diablo_outer_holdout.csv")),
  file.exists(file.path(out,"diablo_signature_stability.rds")),
  identical(result$summary$signatureStatus,
    "internal_repeatability_only_not_validated_biomarker"),
  length(result$summary$signatureStability)==2L,
  all(vapply(result$summary$signatureStability,
    function(x)is.finite(x$meanPairwiseJaccard) &&
      x$meanPairwiseJaccard>=0 && x$meanPairwiseJaccard<=1 &&
      x$nPairwiseComparisons==1L,logical(1)))
)
# Negative control with permuted labels, no label-informed preprocessing
set.seed(20261009)
shuffled <- sample(unname(y))
names(shuffled) <- names(y)
control <- run_diablo_outer_holdout(
  blocks,shuffled,tempfile("diablo_null_"),
  seeds=c(20261009L,20261010L),tune=FALSE
)
stopifnot(
  identical(control$status,"ok"),
  is.finite(control$summary$meanBER),
  control$summary$meanBER>=0 && control$summary$meanBER<=1,
  all(vapply(control$summary$signatureStability,
    function(x) is.finite(x$meanPairwiseJaccard),logical(1)))
)
# Eight complete subject-label permutations: fold-specific preprocessing,
# model selection and fitting are recomputed for EVERY null dataset. This is
# a reproducible *smoke* null distribution, not a 5% significance analysis:
# resolution cannot be lower than 1/9 (0.111...).
permutation_dir <- tempfile("diablo_permutations_")
negative <- run_diablo_permutation_benchmark(
  blocks,y,permutation_dir,n_permutations=8L,
  permutation_seed=20261009L,
  split_seeds=c(20261009L,20261010L),
  tune=FALSE
)
stopifnot(
  identical(negative$status,"exploratory_null_distribution_not_formal_calibration"),
  identical(negative$nPermutations,8L),
  abs(negative$minimumPossibleP-1/9)<1e-12,
  length(negative$nullBERValues)==8L,
  all(is.finite(negative$nullBERValues)),
  all(vapply(negative$featureStability,function(x)
    length(x$nullValues)==8L && all(is.finite(x$nullValues)) &&
    x$upperTailEmpiricalP>=negative$minimumPossibleP,logical(1))),
  file.exists(file.path(permutation_dir,"diablo_permutation_null.rds"))
)
unlink(permutation_dir,recursive=TRUE)

cat(sprintf(
  "Actual mixOmics DIABLO | known signal heldout BER %.3f | permuted-outcome heldout BER %.3f | per-train tuning/preprocessing and complete predictions PASS\n",
  result$summary$meanBER,control$summary$meanBER
))
unlink(out,recursive=TRUE)
