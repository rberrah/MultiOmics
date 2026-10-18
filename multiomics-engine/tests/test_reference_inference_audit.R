#!/usr/bin/env Rscript
# Package-free, adversarial tests of native feature identity, BH provenance
# and algorithmic positive/negative controls. Not an FDR certificate.
source("multiomics-engine/reference_inference_audit.R")

for(vector in list(
  c(0,0,.01,.01,1),
  c(.9,.9,.8,.1,.05,.02,.02),
  c(.05),
  c(.00001,.0001,.5,.9,1))) {
  stopifnot(isTRUE(all.equal(
    unname(independent_bh_adjust(vector)),
    unname(stats::p.adjust(vector,method="BH")),
    tolerance=1e-14)))
}
for(bad_input in list(c(.1,NA),c(.1,Inf),c(-.1,.2),c(.1,1.1),numeric())) {
  stopifnot(inherits(try(independent_bh_adjust(bad_input),silent=TRUE),
    "try-error"))
}

features <- sprintf("GENE%03d",seq_len(10L))
p <- c(.001,.002,.02,.08,.1,.3,.4,.5,.7,.95)
q <- stats::p.adjust(p,method="BH")
base <- data.frame(feature=features,p=p,q=q,stringsAsFactors=FALSE)
check <- function(x,expected_features=features) audit_native_r_bh(
  x,"DESeq2","test",1L,"feature","p","q",expected_features)
out <- check(base)
stopifnot(nrow(out)==1L,out$source_features==10L,
  out$finite_raw_p==10L,out$finite_adjusted_q==10L,
  out$excluded_raw_p_with_no_q==0L,
  out$bh_discoveries_005==sum(q<=.05),
  out$replay_max_abs_q_difference==0,
  out$native_bh_integrity_verified,
  !out$fdr_calibration_certified,!out$biological_study_certified)

# DESeq2 independent filtering: valid p remains measured, but q is missing.
filtered <- base
filtered$q[8:10] <- NA_real_
filtered$q[1:7] <- stats::p.adjust(filtered$p[1:7],method="BH")
filtered_result <- check(filtered)
stopifnot(filtered_result$finite_raw_p==10L,
  filtered_result$finite_adjusted_q==7L,
  filtered_result$excluded_raw_p_with_no_q==3L)
reordered <- filtered[10:1,,drop=FALSE]
stopifnot(check(reordered)$replay_max_abs_q_difference<1e-12)

reject <- function(bad,ids=features) inherits(try(
  check(bad,ids),silent=TRUE),"try-error")
bad <- base;bad$q[1L] <- .6;stopifnot(reject(bad))
bad <- base;bad$q[2L] <- NA_real_;stopifnot(reject(bad))
# That q has been removed requires recomputing BH on kept p:
# otherwise a screenshot-like q-table with dropped rows passes silently.
bad <- base;bad$feature[2L] <- "UNKNOWN";stopifnot(reject(bad))
bad <- base;bad$feature[2L] <- bad$feature[1L];stopifnot(reject(bad))
bad <- base;bad$p[1L] <- Inf;stopifnot(reject(bad))
bad <- base;bad$p[1L] <- -0.01;stopifnot(reject(bad))
bad <- base;bad$q[1L] <- Inf;stopifnot(reject(bad))
bad <- base;bad$p[1L] <- NA_real_;stopifnot(reject(bad))
bad <- base;bad$q[] <- NA_real_;stopifnot(reject(bad))
bad <- base;bad$p[1L] <- .9;stopifnot(reject(bad))
bad <- base;bad$q[1L] <- NA_real_;bad$p[1L] <- NA_real_;
stopifnot(reject(bad))
bad <- base;bad$q[1L] <- 0;stopifnot(reject(bad))
stopifnot(reject(base,features[-1L]))

# Preserve the caller's RNG state despite seeded synthetic experiments.
set.seed(12345L);before <- get(".Random.seed",envir=.GlobalEnv)
oracle <- run_independent_bh_oracle(replicates=512L,features=100L)
stopifnot(identical(before,get(".Random.seed",envir=.GlobalEnv)),
  nrow(oracle$per_study)==1024L,nrow(oracle$summary)==2L,
  length(unique(oracle$per_study$replicate))==512L,
  all(oracle$summary$studies==512L),
  all(oracle$summary$exact_mc_95_high>
      oracle$summary$exact_mc_95_low),
  all(!oracle$summary$native_method_validated),
  all(!oracle$summary$biological_study_certified))
stopifnot(inherits(try(run_independent_bh_oracle(
  replicates=4L),silent=TRUE),"try-error"))
cat("PASS: R-native BH inference integrity, independent filtering and seeded oracle negative controls; no scientific certification\n")
