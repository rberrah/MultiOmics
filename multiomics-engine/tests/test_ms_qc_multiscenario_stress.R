#!/usr/bin/env Rscript
# Independent synthetic known-truth MS QC stress fixtures. No package needed.
# Genuine four-file faahKO LC-MS data are validated in a different workflow
# and do NOT contain pooled-QC or blank injections.
source("multiomics-engine/ms_quantitative_assurance.R")
ids<-sprintf("INJ%03d",seq_len(24L))
features<-c("clean","contamination","drift","batch_shift","treated_missing","qc_absent")
meta<-data.frame(
 assay_id=ids,
 sample_type=c(rep("biological",12L),rep("pooled_qc",8L),rep("blank",4L)),
 condition=c(rep("control",6L),rep("treated",6L),rep(NA_character_,12L)),
 injection_order=seq_len(24L),
 batch=c(rep(c("B1","B2"),6L),rep(c("B1","B2"),each=4L),rep(c("B1","B2"),each=2L)),
 stringsAsFactors=FALSE
)
x<-matrix(100,24,6,dimnames=list(ids,features))
qc<-which(meta$sample_type=="pooled_qc")
blanks<-which(meta$sample_type=="blank")
treated<-which(meta$sample_type=="biological" & meta$condition=="treated")
x[blanks,]<-1
x[qc,]<-matrix(rep(c(100,99,101,100,98,102,101,99),length(features)),
                length(qc),length(features))
# Known anomalies are planted in *different features*, with exact
# quantitative expectations that do not depend on a random seed.
x[blanks,"contamination"]<-40
x[qc,"drift"]<-seq(70,112,length.out=8L)
x[qc,"batch_shift"]<-c(99,100,101,100,159,160,161,160)
x[treated,"treated_missing"]<-NA_real_
x[qc,"qc_absent"]<-NA_real_
summary<-ms_quantitative_assurance(x,meta)
f<-summary$features
row<-function(id)match(id,f$feature_id)
stopifnot(!summary$scientific_certified,
 !summary$feature_exclusions_automatically_applied,
 !summary$coverage$batch_condition_fully_confounded,
 summary$coverage$biological_assays==12L,
 summary$coverage$pooled_qc_assays==8L,
 summary$coverage$blank_assays==4L,
 summary$coverage$pooled_qc_batch_screen_evaluable,
 summary$status=="scientific_review_required",
 !any(unlist(f[row("clean"),
    c("flag_blank_signal","flag_qc_variation","flag_qc_drift_screen",
      "flag_qc_batch_shift_screen","flag_missing","flag_group_missingness")])))

stopifnot(
 isTRUE(f$flag_blank_signal[row("contamination")]),
 isTRUE(all.equal(f$blank_over_biological_ratio[row("contamination")],
                  0.4,tolerance=1e-12)),
 f$pooled_qc_drift_status[row("drift")]=="screen_evaluable",
 isTRUE(f$flag_qc_drift_screen[row("drift")]),
 !f$flag_qc_variation[row("drift")],
 isTRUE(all.equal(f$pooled_qc_order_spearman[row("drift")],1,tolerance=1e-12)),
 f$pooled_qc_batch_status[row("batch_shift")]=="screen_evaluable",
 isTRUE(f$flag_qc_batch_shift_screen[row("batch_shift")]),
 isTRUE(all.equal(f$pooled_qc_batch_median_ratio[row("batch_shift")],
                  1.6,tolerance=1e-12)),
 isTRUE(f$flag_group_missingness[row("treated_missing")]),
 isTRUE(f$flag_missing[row("treated_missing")]),
 isTRUE(all.equal(f$missingness_group_difference[row("treated_missing")],
                   1,tolerance=1e-12)),
 isTRUE(f$flag_qc_not_evaluable[row("qc_absent")]),
 isTRUE(f$flag_qc_drift_not_evaluable[row("qc_absent")]),
 isTRUE(f$flag_qc_batch_not_evaluable[row("qc_absent")]),
 !f$flag_qc_batch_shift_screen[row("qc_absent")]
)
# Independent scenario 2: insufficient QC/blank coverage must be prominent.
no_controls<-ms_quantitative_assurance(
 x[1:12,,drop=FALSE],meta[1:12,,drop=FALSE])
stopifnot(!no_controls$coverage$pooled_qc_precision_evaluable,
 !no_controls$coverage$blank_background_evaluable,
 !no_controls$coverage$pooled_qc_drift_screen_evaluable,
 !no_controls$coverage$pooled_qc_batch_screen_evaluable,
 !no_controls$scientific_certified,
 "pooled_qc_not_sufficient_for_rsd"%in%no_controls$alerts,
 "blanks_not_sufficient_for_background"%in%no_controls$alerts)
# Independent scenario 3: completely confounded biological groups/batches.
confounded<-meta
confounded$batch[1:12]<-confounded$condition[1:12]
bad<-ms_quantitative_assurance(x,confounded)
stopifnot(isTRUE(bad$coverage$batch_condition_fully_confounded),
 "batch_perfectly_confounded_with_biological_group"%in%bad$alerts,
 !bad$scientific_certified)
# Scenario 4: a pooled QC label does not make a missing feature evaluable.
partial<-x
partial[qc[1:6],"clean"]<-NA_real_
flag<-ms_quantitative_assurance(partial,meta)$features
stopifnot(flag$pooled_qc_precision_status[row("clean")] ==
 "too_few_finite_pooled_qc_values",
 isTRUE(flag$flag_qc_not_evaluable[row("clean")]),
 isTRUE(flag$flag_qc_batch_not_evaluable[row("clean")]))
# Reject meaningless thresholds and invalid injection data.
bad_threshold<-try(ms_quantitative_assurance(
 x,meta,qc_batch_ratio_screen_threshold=0.9),silent=TRUE)
stopifnot(inherits(bad_threshold,"try-error"))
bad_order<-meta
bad_order$injection_order[24]<-bad_order$injection_order[23]
stopifnot(inherits(try(ms_quantitative_assurance(x,bad_order),silent=TRUE),"try-error"))
out<-Sys.getenv("MULTIOMICS_MS_QC_OUTPUT","")
if(nzchar(out)){
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 utils::write.csv(f,file.path(out,"synthetic_known_truth_feature_qc.csv"),row.names=FALSE,na="")
 utils::write.csv(meta,file.path(out,"synthetic_known_truth_manifest.csv"),row.names=FALSE,na="")
 writeLines(c(
  "These are SYNTHETIC known-truth intensity matrices, NOT human patient data.",
  "Three isolated experimental stress scenarios plus a partially unobserved QC feature.",
  "QC batch/monotone-run-order correlations are prespecified DESCRIPTIVE screens.",
  "No p-values, FDR or clinical/compound validation are claimed.",
  paste("source_commit:",Sys.getenv("GITHUB_SHA","local")),
  "A separate run validates four REAL faahKO CDF files without QC or blanks."
 ),file.path(out,"INTERPRETATION.txt"))
}
cat("PASS: synthetic MS QC blank contamination, monotone drift, replicated QC batches, no controls, confounding, feature observability, thresholds\n")
