#!/usr/bin/env Rscript
# True pooled-QC external mzML data, published via Zenodo + Galaxy GTN.
# 6 biological human urine runs + 3 pooled QC runs, no blank measurements.
source("multiomics-engine/raw_ms_processing.R")
source("multiomics-engine/ms_quantitative_assurance.R")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: real_sacurine_zenodo_qc_xcms.R <source dir> <output dir>")
names<-c("HU_neg_048","HU_neg_090","HU_neg_123","HU_neg_157",
         "HU_neg_173","HU_neg_192","QC1_002","QC1_008","QC1_014")
role<-c(rep("biological",6L),rep("pooled_qc",3L))
files<-file.path(args[[1L]],paste0(names,".mzML"))
if(any(!file.exists(files)))stop("One or more original Zenith/Galaxy mzML files missing")
out<-args[[2L]]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Refuse reuse of stale native raw results")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
metadata<-data.frame(sample_id=names,assay_id=names,sample_type=role,
 condition=c(rep("urine",6L),rep(NA_character_,3L)),
 injection_order=seq_along(names),injection_order_verified=rep(FALSE,9L),
 stringsAsFactors=FALSE)
params<-list(peak_detection=list(ppm=25,peakwidth=c(10,60),snthresh=8,
                                     prefilter=c(3,100)),
 retention_time=list(enabled=TRUE,binSize=.6),
 correspondence=list(bw=8,minFraction=.5,minSamples=1L,binSize=.025,ppm=15))
res<-run_xcms_raw_pipeline(files,metadata,out,parameters=params,chunk_size=2L)
stopifnot(res$summary$samples==9L,res$summary$features>0L,
 res$summary$chromatographic_peaks>0L,
 ncol(res$matrix)==9L,nrow(res$matrix)==res$summary$features,
 !anyDuplicated(rownames(res$matrix)))
native<-t(res$matrix)
rownames(native)<-as.character(colnames(res$matrix))
qa<-ms_quantitative_assurance(native,metadata)
stopifnot(qa$coverage$biological_assays==6L,
 qa$coverage$pooled_qc_assays==3L,
 qa$coverage$blank_assays==0L,
 !qa$coverage$pooled_qc_drift_screen_evaluable,
 !qa$coverage$pooled_qc_injection_order_verified,
 !qa$coverage$blank_background_evaluable,
 !qa$scientific_certified,
 !qa$feature_exclusions_automatically_applied,
 "pooled_qc_injection_order_not_verified"%in%qa$alerts || 
   "pooled_qc_not_sufficient_for_drift"%in%qa$alerts)
rds<-readRDS(file.path(out,"metabolomics_feature_definitions.rds"))
csv<-utils::read.csv(file.path(out,"metabolomics_peak_area.csv"),
                     check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(csv)==nrow(rds),setequal(csv$feature_id,rownames(res$matrix)))
actual<-as.matrix(csv[,-1L,drop=FALSE])
storage.mode(actual)<-"double"
stopifnot(isTRUE(all.equal(unname(actual),unname(res$matrix[csv$feature_id,,drop=FALSE]),
                           check.attributes=FALSE,tolerance=1e-8)))
f<-qa$features
if(!requireNamespace("jsonlite",quietly=TRUE))stop("Missing reference jsonlite")
utils::write.csv(f,file.path(out,"real_sacurine_qc_quality_features.csv"),
                 row.names=FALSE,na="")
jsonlite::write_json(list(
  schema="multiomics.real-sacurine-pooled-qc-xcms-evidence.v1",
  record="https://zenodo.org/records/3244991",
  source_context="Galaxy Training Network real nine-file LC-MS tutorial",
  real_biological_runs=6,real_pooled_qc_runs=3,real_blank_runs=0,
  native_chromatographic_peaks=res$summary$chromatographic_peaks,
  native_grouped_features=res$summary$features,
  qc_rsd_evaluable_features=sum(is.finite(f$pooled_qc_rsd)),
  real_qc_rsd_gt30pct=sum(is.finite(f$pooled_qc_rsd)&f$pooled_qc_rsd>.3),
  injection_order_verified=FALSE,drift_evaluable=FALSE,
  blank_contamination_evaluable=FALSE,
  biological_group_statistic_performed=FALSE,
  scientific_certified=FALSE,
  source_sha=Sys.getenv("GITHUB_SHA","local"),
  limitation="Original six urine plus three pooled QCs, no blanks or verified injection order. QC precision is descriptive, not an independent clinical/FDR validation."
),file.path(out,"real_sacurine_qc_evidence.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS real original QC mzML: ",res$summary$chromatographic_peaks,
    " chromatographic peaks; ",res$summary$features," grouped features; ",
    sum(is.finite(f$pooled_qc_rsd))," evaluated real QC RSDs\n",sep="")
