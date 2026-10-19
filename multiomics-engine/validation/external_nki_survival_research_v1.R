#!/usr/bin/env Rscript
# Real public external execution qualification for MultiOmics native Cox path.
# Dataset: breastCancerNKI. Endpoint is prespecified DMFS (t.dmfs/e.dmfs).
# Feature family is selected by identifier/measurement validity only, before
# any survival fit. This is wrapper/model-contract qualification, NOT biomarker
# validation, FDR certification, causal inference or clinical validation.

suppressPackageStartupMessages({
  library(Biobase)
  library(breastCancerNKI)
  library(survival)
  library(jsonlite)
  library(digest)
})
source("multiomics-engine/advanced_methods.R")

out <- Sys.getenv("MULTIOMICS_EXTERNAL_SURVIVAL_OUT", "tmp/external-survival-nki")
if (dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE)))
  stop("Refusing to overwrite prior external NKI evidence.")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

data(nki)
x_all <- Biobase::exprs(nki)
pd <- Biobase::pData(nki)
if (!all(c("t.dmfs","e.dmfs") %in% names(pd)))
  stop("NKI DMFS fields are absent; do not infer an alternative endpoint.")

time <- suppressWarnings(as.numeric(pd$t.dmfs))
event <- suppressWarnings(as.numeric(pd$e.dmfs))
eligible <- is.finite(time) & time > 0 & is.finite(event) & event %in% c(0,1)
if (sum(eligible) < 250L || sum(event[eligible]) < 40L)
  stop("Unexpected NKI DMFS completeness/event count; source contract changed.")

sample_ids <- sampleNames(nki)[eligible]
expr <- x_all[, eligible, drop=FALSE]
if (!identical(colnames(expr), sample_ids))
  stop("NKI expression/clinical subject order mismatch.")
if (anyDuplicated(sample_ids) || any(!nzchar(sample_ids)))
  stop("NKI subject identifiers must be unique and nonempty.")

# Prespecified family: first 64 lexicographically sorted source feature IDs
# that are completely finite and nonconstant in the eligible DMFS population.
# No endpoint value, coefficient or p-value participates in feature selection.
ids <- rownames(expr)
measurement_ok <- vapply(seq_len(nrow(expr)), function(i) {
  z <- expr[i,]
  all(is.finite(z)) && stats::sd(z) > 1e-10
}, logical(1))
candidate <- sort(ids[measurement_ok])
if (length(candidate) < 64L) stop("Too few complete nonconstant NKI features.")
features <- candidate[seq_len(64L)]
matrix <- t(expr[features,,drop=FALSE])
storage.mode(matrix) <- "double"
rownames(matrix) <- sample_ids
colnames(matrix) <- features

meta <- data.frame(
  sample_id=sample_ids,
  subject_id=sample_ids,
  survival_time=time[eligible],
  survival_event=event[eligible],
  stringsAsFactors=FALSE
)
rownames(meta) <- sample_ids

wrapper <- run_cox_survival_reference(
  matrix=matrix, metadata=meta,
  output_dir=file.path(out,"wrapper"),
  covariates=character(), ties="efron", min_events=12L
)
tab <- wrapper$results
if (nrow(tab) != length(features) || !setequal(tab$feature,features))
  stop("Native Cox wrapper did not preserve the complete prespecified feature family.")

# Independent direct survival::coxph replay feature by feature.
direct <- lapply(features, function(feature) {
  dat <- meta
  dat$assay_feature <- matrix[,feature]
  fit <- survival::coxph(
    survival::Surv(survival_time,survival_event) ~ assay_feature,
    data=dat,ties="efron",x=TRUE,
    control=survival::coxph.control(iter.max=60,timefix=FALSE)
  )
  beta <- unname(stats::coef(fit)["assay_feature"])
  se <- sqrt(stats::vcov(fit)["assay_feature","assay_feature"])
  zph <- survival::cox.zph(fit,transform="km")
  data.frame(
    feature=feature,
    effect=beta,
    standardError=se,
    pValue=2*stats::pnorm(-abs(beta/se)),
    hazardRatio=exp(beta),
    ciLow=beta-1.959963984540054*se,
    ciHigh=beta+1.959963984540054*se,
    proportionalHazardsP=unname(zph$table["assay_feature","p"]),
    globalPhP=unname(zph$table["GLOBAL","p"]),
    stringsAsFactors=FALSE
  )
})
direct <- do.call(rbind,direct)
direct$qValue <- stats::p.adjust(direct$pValue,method="BH")
utils::write.csv(direct,file.path(out,"independent_direct_coxph.csv"),row.names=FALSE)

w <- tab[match(features,tab$feature),,drop=FALSE]
d <- direct[match(features,direct$feature),,drop=FALSE]
if (anyNA(w$feature) || anyNA(d$feature)) stop("External Cox feature matching failed.")
fields <- c("effect","standardError","pValue","qValue","hazardRatio",
            "ciLow","ciHigh","proportionalHazardsP","globalPhP")
max_abs <- setNames(numeric(length(fields)),fields)
for (name in fields) {
  delta <- abs(as.numeric(w[[name]])-as.numeric(d[[name]]))
  if (any(!is.finite(delta))) stop(paste("Non-finite external Cox comparison:",name))
  max_abs[[name]] <- max(delta)
  if (max_abs[[name]] > 1e-9)
    stop(sprintf("NKI Cox wrapper/direct mismatch %s: %.12g",name,max_abs[[name]]))
}

ph_warning <- sum(w$status=="ph_assumption_warning")
if (any(!w$status %in% c("ph_test_not_rejected","ph_assumption_warning")))
  stop("Unexpected failed Cox fit in the frozen NKI feature family.")

feature_manifest <- data.frame(
  rank=seq_along(features), feature=features,
  selection_rule="lexicographic_first_64_complete_nonconstant_before_outcome_fit",
  stringsAsFactors=FALSE
)
utils::write.csv(feature_manifest,file.path(out,"prespecified_feature_family.csv"),row.names=FALSE)
utils::write.csv(meta,file.path(out,"dmfs_metadata_used.csv"),row.names=FALSE)

session <- capture.output(sessionInfo())
writeLines(session,file.path(out,"sessionInfo.txt"))
report <- list(
  schema="multiomics.external-survival-nki-research-v1.v1",
  source=list(
    package="breastCancerNKI",
    packageVersion=as.character(utils::packageVersion("breastCancerNKI")),
    doi="10.18129/B9.bioc.breastCancerNKI",
    originalStudies=c("van't Veer et al. 2002","van de Vijver et al. 2002"),
    endpoint="distant metastasis-free survival (t.dmfs/e.dmfs)"
  ),
  prespecifiedSelection="first 64 lexicographically sorted feature IDs with finite nonconstant measurements among subjects with complete DMFS; no outcome statistics used",
  subjects=nrow(meta),
  events=sum(meta$survival_event),
  attemptedFeatures=length(features),
  estimableFeatures=sum(is.finite(w$pValue)),
  phWarnings=ph_warning,
  maxAbsoluteWrapperDirectDifference=as.list(max_abs),
  wrapperMethod=wrapper$summary$method,
  resultInterpretation="External real-data numerical execution qualification only. PH warnings are preserved; non-significance does not prove PH. No biomarker, FDR, causal or clinical validity claim.",
  softwareQualified=TRUE,
  statisticalCalibrationCertified=FALSE,
  biologicalValidation=FALSE,
  clinicalCertification=FALSE
)
jsonlite::write_json(report,file.path(out,"external_survival_nki_evidence.json"),
                     pretty=TRUE,auto_unbox=TRUE,null="null")
cat(sprintf(
  "PASS NKI external Cox qualification: %d subjects, %d events, %d prespecified features, %d PH warnings, max numeric delta %.3g; NO biomarker/FDR/clinical certification.\n",
  nrow(meta),sum(meta$survival_event),length(features),ph_warning,
  max(unlist(max_abs))
))
