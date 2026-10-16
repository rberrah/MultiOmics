#!/usr/bin/env Rscript
source("multiomics-engine/ms_quantitative_assurance.R")
ids<-paste0("MS",seq_len(17L))
meta<-data.frame(assay_id=ids,
 sample_type=c(rep("biological",10),rep("pooled_qc",5),rep("blank",2)),
 condition=c(rep(c("control","treated"),each=5L),rep("",7)),
 batch=c(rep(c("B1","B2"),5),rep("B1",7)),
 injection_order=seq_len(17L),stringsAsFactors=FALSE)
x<-matrix(100,17L,4L,dimnames=list(ids,c("clean","background","unstable_qc","missing_in_treated")))
x[1:10,"background"]<-50
x[16:17,"background"]<-35
x[11:15,"background"]<-50
x[11:15,"unstable_qc"]<-c(30,150,60,190,10)
x[6:10,"missing_in_treated"]<-NA_real_
a<-ms_quantitative_assurance(x,meta)
# Independently hand-check numerical estimands against reference R formulas.
reference<-x
reference[11:15,"background"]<-c(20,25,30,35,40)
numeric_check<-ms_quantitative_assurance(reference,meta)
stopifnot(isTRUE(all.equal(numeric_check$features$pooled_qc_rsd[2],
  stats::sd(c(20,25,30,35,40))/mean(c(20,25,30,35,40)),
  tolerance=1e-12)),
 isTRUE(all.equal(numeric_check$features$blank_over_biological_ratio[2],
  35/50,tolerance=1e-12)),
 isTRUE(all.equal(numeric_check$features$missingness_group_difference[4],
  1.0,tolerance=1e-12)),
 isTRUE(all.equal(numeric_check$features$pooled_qc_order_spearman[2],
  1.0,tolerance=1e-12)))

# A documented QC injection does not mean each molecular feature has
# enough finite measurements to estimate precision/background.
stopifnot(all(c("pooled_qc_precision_status",
  "blank_background_status","pooled_qc_drift_status",
  "finite_biological_values","finite_pooled_qc_values",
  "finite_blank_values") %in% names(a$features)),
  a$features$pooled_qc_precision_status[1]=="rsd_evaluable",
  a$features$blank_background_status[1]=="ratio_evaluable")

# An analyte abundant in blanks but absent in the biological sample has
# an undefined blank/biological ratio; it must be flagged, not silently NA=OK.
blank_only<-x
blank_only[1:10,"clean"]<-0
blank_only[16:17,"clean"]<-50
blank_audit<-ms_quantitative_assurance(blank_only,meta)
stopifnot(is.na(blank_audit$features$blank_over_biological_ratio[1]),
 blank_audit$features$blank_background_status[1]=="positive_blank_with_zero_biological_median",
 isTRUE(blank_audit$features$flag_blank_signal[1]))

# Control injections recorded but the tested feature itself unobservable.
qc_missing<-x
qc_missing[11:15,"clean"]<-NA_real_
qc_audit<-ms_quantitative_assurance(qc_missing,meta)
stopifnot(is.na(qc_audit$features$pooled_qc_rsd[1]),
 qc_audit$features$finite_pooled_qc_values[1]==0,
 qc_audit$features$pooled_qc_precision_status[1]=="too_few_finite_pooled_qc_values",
 qc_audit$features$pooled_qc_precision_status[2]=="rsd_evaluable",
 length(qc_audit$features$pooled_qc_precision_status)==ncol(x),
 isTRUE(qc_audit$features$flag_qc_not_evaluable[1]),
 !qc_audit$features$flag_qc_not_evaluable[2],
 any(grepl("flag_qc_not_evaluable",qc_audit$alerts)))

blank_missing<-x
blank_missing[16:17,"clean"]<-NA_real_
blank_missing_audit<-ms_quantitative_assurance(blank_missing,meta)
stopifnot(blank_missing_audit$features$blank_background_status[1]=="too_few_finite_blank_values",
 isTRUE(blank_missing_audit$features$flag_blank_not_evaluable[1]),
 any(grepl("flag_blank_not_evaluable",blank_missing_audit$alerts)))

zero_qc<-x
zero_qc[11:15,"clean"]<-0
zero_qc_audit<-ms_quantitative_assurance(zero_qc,meta)
stopifnot(zero_qc_audit$features$pooled_qc_precision_status[1]=="pooled_qc_mean_zero_or_not_estimable",
 isTRUE(zero_qc_audit$features$flag_qc_not_evaluable[1]))

both_zero<-x
both_zero[1:10,"clean"]<-0
both_zero[16:17,"clean"]<-0
zero_blank_audit<-ms_quantitative_assurance(both_zero,meta)
stopifnot(zero_blank_audit$features$blank_background_status[1]=="both_blank_and_biological_medians_zero",
 isTRUE(zero_blank_audit$features$flag_blank_not_evaluable[1]))

stopifnot(a$status=="scientific_review_required",
 a$coverage$biological_assays==10L,
 a$coverage$pooled_qc_assays==5L,
 a$coverage$blank_assays==2L,
 isFALSE(a$coverage$batch_condition_fully_confounded),
 !a$scientific_certified,!a$mnar_diagnosis_possible_from_this_table,
 !a$pathway_or_compound_identity_certified,
 !a$feature_exclusions_automatically_applied,
 a$features$flag_blank_signal[2],
 a$features$flag_qc_variation[3],
 a$features$flag_group_missingness[4],
 a$features$flag_missing[4],
 all(!a$features$flag_missing[1:3]))
reject<-function(x,meta)
 inherits(try(ms_quantitative_assurance(x,meta),silent=TRUE),"try-error")
bad<-meta;bad$assay_id[2]<-bad$assay_id[1]
stopifnot(reject(x,bad))
bad<-meta;bad$sample_type[2]<-"unknown"
stopifnot(reject(x,bad))
bad<-meta;bad$injection_order[2]<-bad$injection_order[1]
stopifnot(reject(x,bad))
bad<-x;bad[1,1]<--2
stopifnot(reject(bad,meta))
bad<-x;rownames(bad)[2]<-rownames(bad)[1]
stopifnot(reject(bad,meta))
bad<-meta;bad$batch[1]<-""
stopifnot(reject(x,bad))
# Real laboratory manifests legitimately leave biological condition empty
# on pooled-QC and solvent blank injections. Never discard or reject those
# controls just because they have no treatment group.
optional_groups<-meta
optional_groups$condition[11:17]<-NA_character_
stopifnot(ms_quantitative_assurance(x,optional_groups)$coverage$blank_assays==2L)
bad_group<-meta
bad_group$condition[1L]<-NA_character_
stopifnot(reject(x,bad_group))
without_controls<-ms_quantitative_assurance(x[1:10,,drop=FALSE],meta[1:10,,drop=FALSE])
stopifnot(without_controls$status=="scientific_review_required",
 "pooled_qc_not_sufficient_for_rsd"%in%without_controls$alerts,
 "blanks_not_sufficient_for_background"%in%without_controls$alerts,
 !without_controls$coverage$pooled_qc_precision_evaluable)
confounded<-meta
confounded$batch[1:10]<-confounded$condition[1:10]
flag<-ms_quantitative_assurance(x,confounded)
stopifnot(isTRUE(flag$coverage$batch_condition_fully_confounded),
 "batch_perfectly_confounded_with_biological_group"%in%flag$alerts)
# Exercise the real CSV command-line entry point, not only the helper.
temp<-tempfile("ms_audit_cli_")
dir.create(temp)
peak_path<-file.path(temp,"peak.csv")
manifest_path<-file.path(temp,"sample.csv")
out_path<-file.path(temp,"audit")
peak_csv<-data.frame(feature_id=colnames(x),t(x),check.names=FALSE,
                     stringsAsFactors=FALSE)
utils::write.csv(peak_csv,peak_path,row.names=FALSE,na="")
csv_manifest<-meta
csv_manifest$condition[11:17]<-NA_character_
utils::write.csv(csv_manifest,manifest_path,row.names=FALSE,na="")
exec<-suppressWarnings(system2(Sys.which("Rscript"),
   c("multiomics-engine/run_ms_quantitative_assurance.R",
     shQuote(peak_path),shQuote(manifest_path),shQuote(out_path)),
   stdout=TRUE,stderr=TRUE))
stopifnot(is.null(attr(exec,"status"))||attr(exec,"status")==0L,
  file.exists(file.path(out_path,"ms_quantitative_qc_features.csv")),
  file.exists(file.path(out_path,"ms_quantitative_qc_manifest.R")),
  file.exists(file.path(out_path,"ms_quantitative_qc_alerts.csv")))
report_csv<-utils::read.csv(file.path(out_path,"ms_quantitative_qc_features.csv"))
stopifnot(nrow(report_csv)==4L,any(report_csv$flag_group_missingness))
unlink(temp,recursive=TRUE)

# Actual public processed-peak tables may not contain documented QC injection
# order. Never treat arbitrary row position as measured temporal drift.
unverified<-meta
unverified$injection_order_verified<-rep(FALSE,nrow(unverified))
no_order<-ms_quantitative_assurance(x,unverified)
stopifnot(!no_order$coverage$pooled_qc_injection_order_verified,
 !no_order$coverage$pooled_qc_drift_screen_evaluable,
 all(is.na(no_order$features$pooled_qc_order_spearman)),
 all(no_order$features$pooled_qc_drift_status=="injection_order_not_verified"),
 all(!no_order$features$flag_qc_drift_screen),
 "pooled_qc_injection_order_not_verified"%in%no_order$alerts,
 !no_order$scientific_certified)
known_order<-meta
known_order$injection_order_verified<-rep(TRUE,nrow(known_order))
known<-ms_quantitative_assurance(x,known_order)
stopifnot(known$coverage$pooled_qc_injection_order_verified,
 known$coverage$pooled_qc_drift_screen_evaluable)
mixed<-unverified
mixed$injection_order_verified[1]<-TRUE
stopifnot(reject(x,mixed))
na_flag<-unverified
na_flag$injection_order_verified[1]<-NA
stopifnot(reject(x,na_flag))
text_flag<-unverified
text_flag$injection_order_verified<-rep("false",nrow(text_flag))
stopifnot(reject(x,text_flag))

cat("MS QUANTITATIVE QA: pooled QC, blanks, MNAR caution, group missingness, batch and fail-closed IDs PASS\n")
