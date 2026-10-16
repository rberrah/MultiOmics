#!/usr/bin/env Rscript
# Professional optional MS QC from exported peak-area matrices and a declared
# sample sheet. Does NOT run imputation, select "good" features or change UI.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)
  stop("Usage: Rscript multiomics-engine/run_ms_quantitative_assurance.R <peak_areas.csv> <sample_manifest.csv> <output_dir>")
source("multiomics-engine/ms_quantitative_assurance.R")
data<-utils::read.csv(args[[1L]],check.names=FALSE,stringsAsFactors=FALSE,
 na.strings=c("","NA","NaN"))
metadata<-utils::read.csv(args[[2L]],check.names=FALSE,
 stringsAsFactors=FALSE,na.strings=c("","NA","NaN"))
if(!nrow(data)||ncol(data)<3L||
   names(data)[1L]!="feature_id"||
   anyNA(data$feature_id)||anyDuplicated(data$feature_id)||
   any(!nzchar(as.character(data$feature_id))))
  stop("Peak-area matrix needs a unique first feature_id column.")
expected<-as.character(metadata$assay_id)
if(is.null(expected)||anyNA(expected)||anyDuplicated(expected)||
  length(expected)<4L||length(names(data)[-1L])!=length(expected)||
  anyDuplicated(names(data)[-1L])||!setequal(names(data)[-1L],expected))
  stop("Peak-area assay columns must match exactly the MS sample manifest.")
ids<-as.character(data$feature_id)
raw_values<-as.matrix(data[,-1L,drop=FALSE])
feature_by_assay<-matrix(suppressWarnings(as.numeric(raw_values)),
  nrow=nrow(raw_values),ncol=ncol(raw_values),dimnames=dimnames(raw_values))
if(any(!is.na(raw_values)&is.na(feature_by_assay)))
  stop("MS peak-area CSV contains nonnumeric non-missing measurements.")
rownames(feature_by_assay)<-ids
x<-t(feature_by_assay)
x<-x[expected,,drop=FALSE]
results<-ms_quantitative_assurance(x,metadata)
dir.create(args[[3L]],recursive=TRUE,showWarnings=FALSE)
utils::write.csv(results$features,
 file.path(args[[3L]],"ms_quantitative_qc_features.csv"),
 row.names=FALSE,na="")
utils::write.csv(data.frame(code=results$alerts),
 file.path(args[[3L]],"ms_quantitative_qc_alerts.csv"),row.names=FALSE)
manifest<-list(status=results$status,coverage=results$coverage,
 parameters=results$parameters,scientific_certified=FALSE,
 feature_exclusions_automatically_applied=FALSE,
 limitation=results$limitation,r_version=R.version.string,
 source_files=basename(args[1:2]))
dput(manifest,file=file.path(args[[3L]],"ms_quantitative_qc_manifest.R"))
if(requireNamespace("jsonlite",quietly=TRUE))
 jsonlite::write_json(manifest,
  file.path(args[[3L]],"ms_quantitative_qc_manifest.json"),
  auto_unbox=TRUE,pretty=TRUE,null="null")
cat("QUANTITATIVE MS QC:",results$status,
    "| results archived, zero auto-exclusions; scientific review still required\n")
