#!/usr/bin/env Rscript
# No Bioconductor dependency: inputs, pairing and additive accounting.
source("multiomics-engine/rnaseq_normalization_lock.R")
samples <- sprintf("S%03d",seq_len(12L))
x <- matrix(100L,12L,10L,dimnames=list(samples,paste0("G",1:10)))
meta <- data.frame(batch=factor(rep(c("B1","B2"),6L)),
  condition=factor(rep(c("control","treated"),each=6L),
                   levels=c("control","treated")),row.names=samples)
ref <- stats::setNames(seq(100,120,length.out=12L),samples)
reject <- function(x,meta,ref,method="edgeR_QL")
  inherits(try(rnaseq_locked_scale_fit(method,x,meta,ref),
               silent=TRUE),"try-error")
bad <- x
bad[1,1] <- -.1
stopifnot(reject(bad,meta,ref))
bad <- x
bad[1,1] <- 2.7
stopifnot(reject(bad,meta,ref))
stopifnot(reject(x,meta,ref[-1L]))
misnamed <- ref
names(misnamed)[1L] <- "UNKNOWN_SUBJECT"
stopifnot(reject(x,meta,misnamed))
badmeta <- meta
rownames(badmeta) <- rev(rownames(meta))
stopifnot(reject(x,badmeta,ref))
meta_confounded <- meta
meta_confounded$batch <- factor(as.character(meta_confounded$condition))
stopifnot(reject(x,meta_confounded,ref))
stopifnot(reject(x,meta,ref,method="nonreference"))
dummy <- list(dge=list(samples=data.frame(
  lib.size=rep(120L,12L),norm.factors=rep(1,12L),
  row.names=samples)))
stopifnot(identical(
  rnaseq_reference_sample_scales("edgeR_QL",dummy),
  stats::setNames(rep(120,12L),samples)))
records <- list()
for(method in c("DESeq2","edgeR_QL","limma_voom"))
 for(rep in seq_len(10L))
  for(arm in c("reference_native","strict_frozen","strict_native")) {
    pp <- switch(arm,reference_native=.5,strict_frozen=.3,strict_native=.35)
    ff <- switch(arm,reference_native=.15,strict_frozen=.07,strict_native=.05)
    records[[length(records)+1L]] <- data.frame(scenario="case",
      method=method,replicate=rep,seed=1000L+rep,arm=arm,
      power=pp,false_discovery_proportion=ff,n_true=20L)
  }
rows <- do.call(rbind,records)
result<-rnaseq_scale_path_summary(rows)
stopifnot(nrow(result)==3L,
 all(abs(result$power_total_change+.15)<1e-10),
 all(abs(result$power_filter_and_refit_at_frozen_scales+.2)<1e-10),
 all(abs(result$power_renormalize_on_strict_filter-.05)<1e-10),
 all(abs(result$fdr_total_change+.1)<1e-10),
 all(abs(result$fdr_filter_and_refit_at_frozen_scales+.08)<1e-10),
 all(abs(result$fdr_renormalize_on_strict_filter+.02)<1e-10),
 all(!result$certified))
bad <- rows[-1L,,drop=FALSE]
stopifnot(inherits(try(rnaseq_scale_path_summary(bad),silent=TRUE),
  "try-error"))
bad <- rows
bad$seed[1L] <- 99999L
stopifnot(inherits(try(rnaseq_scale_path_summary(bad),silent=TRUE),
  "try-error"))
bad <- rows
bad$arm[1L] <- "unknown"
stopifnot(inherits(try(rnaseq_scale_path_summary(bad),silent=TRUE),
  "try-error"))

# Study-level bootstrap must retain pairing, reject invalid input, restore RNG,
# and NEVER turn uncertainty intervals into nominal FDR certification.
set.seed(1971L)
rng_before<-get(".Random.seed",envir=.GlobalEnv)
ci<-rnaseq_scale_path_paired_uncertainty(rows,minimum=10L,
                                          bootstrap_repeats=400L,seed=73L)
stopifnot(identical(rng_before,get(".Random.seed",envir=.GlobalEnv)),
          nrow(ci)==18L,all(ci$independent_studies==10L),
          all(ci$bootstrap_repeats==400L),
          all(ci$bootstrap_seed==73L),
          all(ci$percentile_95_lower<=ci$mean_paired_difference+1e-12),
          all(ci$percentile_95_upper>=ci$mean_paired_difference-1e-12),
          !any(ci$clinical_certification),
          !any(ci$nominal_fdr_certified),
          all(ci$interval_scope=="paired_iid_simulated_studies_descriptive_only"))
# Constant paired differences should have zero bootstrap width even under
# independent resampling of the *study* IDs.
stopifnot(all(abs(ci$percentile_95_lower-ci$mean_paired_difference)<1e-12),
          all(abs(ci$percentile_95_upper-ci$mean_paired_difference)<1e-12),
          identical(ci,rnaseq_scale_path_paired_uncertainty(
              rows,minimum=10L,bootstrap_repeats=400L,seed=73L)))
bad_ci<-rows[-1L,,drop=FALSE]
stopifnot(inherits(try(rnaseq_scale_path_paired_uncertainty(
  bad_ci,bootstrap_repeats=400L),silent=TRUE),"try-error"))
stopifnot(inherits(try(rnaseq_scale_path_paired_uncertainty(
  rows,bootstrap_repeats=0L),silent=TRUE),"try-error"))
stopifnot(inherits(try(rnaseq_scale_path_paired_uncertainty(
  rows,bootstrap_repeats=400L,seed=NA_integer_),silent=TRUE),"try-error"))
# Deliberate between-study heterogeneity creates a strictly positive width.
vary<-rows
vary$power[vary$arm=="strict_native"]<-
  pmax(0,pmin(1,.35+rep(seq(-.15,.15,length.out=10L),3L)))
ci_var<-rnaseq_scale_path_paired_uncertainty(vary,minimum=10L,
                                             bootstrap_repeats=1000L,seed=73L)
s<-ci_var[ci_var$method=="DESeq2"&ci_var$metric=="power"&
          ci_var$comparison=="total",]
stopifnot(nrow(s)==1L,s$percentile_95_upper>s$percentile_95_lower)

cat("RNA-seq frozen-scale validation, paired decomposition and fail-closed guards: PASS\n")
