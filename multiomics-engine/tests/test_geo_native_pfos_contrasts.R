#!/usr/bin/env Rscript
# Independent native R computation for GSE22871 fixed genotype-stratified contrasts.
# Works only from original individual measurements and registered design; not PK.
args <- commandArgs(trailingOnly=TRUE)
if (identical(args,"--self-test")) {
  stopifnot(identical(p.adjust(c(46/252,6/252),method="BH"),c(46/252,12/252)),
            identical(p.adjust(c(.2,.01,.01),method="BH"),c(.2,.015,.015)))
  cat("PASS: independent base-R BH unsorted/duplicate/null safeguards\n")
  quit(save="no")
}
if (length(args)!=1L) stop("Usage: Rscript test_geo_native_pfos_contrasts.R GEO_OUTPUT_DIR")
dir <- args[[1]]
infile <- file.path(dir,"independent-bioc-concordance","native_external_probe_observations",
                    "native_Cyp3a11_probe_measurements_not_model_scores.csv")
outdir <- file.path(dir,"external_native_PFOS_exploratory_contrasts")
reportfile <- file.path(outdir,"native_contrast_report.json")
if (!all(file.exists(c(infile,reportfile))) ||
    !requireNamespace("jsonlite",quietly=TRUE))
  stop("Missing qualified source data or jsonlite")
source <- utils::read.csv(infile,stringsAsFactors=FALSE,check.names=FALSE)
report <- jsonlite::fromJSON(reportfile)
stopifnot(nrow(source)==30L,length(unique(source$assay_id))==30L,
          length(unique(source$probe_id))==1L,
          all(source$cross_platform_model_score=="NOT_COMPUTED"),
          report$frozen_Nutrimouse_model_scores_generated==0L,
          report$n_primary_contrasts==2L,report$n_mice==30L)
results <- report$results
contrasts <- c("wild-type","PPARalpha-null")
stopifnot(nrow(results)==2L,identical(as.character(results$genotype),contrasts))
p <- numeric(2L)
effects <- numeric(2L)
for (i in seq_along(contrasts)) {
  g <- contrasts[[i]]
  high <- source$native_processed_log2_RMA_intensity[
    source$genotype==g & source$pfos_dose_mg_kg_day==10]
  zero <- source$native_processed_log2_RMA_intensity[
    source$genotype==g & source$pfos_dose_mg_kg_day==0]
  stopifnot(length(high)==5L,length(zero)==5L,all(is.finite(c(high,zero))))
  x <- c(high,zero)
  observed <- abs(mean(high)-mean(zero))
  assignments <- utils::combn(10L,5L)
  null <- apply(assignments,2L,function(sel)
    abs(mean(x[sel])-mean(x[-sel])))
  stopifnot(length(null)==252L)
  p[[i]] <- sum(null >= observed-1e-12*max(1,observed))/252
  effects[[i]] <- mean(high)-mean(zero)
}
q <- stats::p.adjust(p,method="BH")
stopifnot(max(abs(p-results$exact_two_sided_permutation_p)) < 1e-12,
          max(abs(effects-results$high_minus_zero_native_log2_RMA)) < 1e-12,
          max(abs(q-results$BH_q_two_predefined_contrasts)) < 1e-12)
lines <- sprintf("%s: delta=%.10f exact p=%.10f BH q=%.10f",
                 contrasts,effects,p,q)
writeLines(c("Independent R exact-permutation/BH replay PASS",lines,
             "No frozen Nutrimouse score or PK endpoint"),
           file.path(outdir,"independent_R_contrast_replay.txt"))
cat("PASS: independent native R exact p-values, log2 effects and correctly ordered BH q-values\n")
cat(paste(lines,collapse="\n"),"\n")
