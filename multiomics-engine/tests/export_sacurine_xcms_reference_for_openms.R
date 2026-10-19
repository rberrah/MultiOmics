#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: export_sacurine_xcms_reference_for_openms.R <xcms-dir> <out-dir>")
root<-args[[1L]];out<-args[[2L]]
required<-c("metabolomics_peak_area.csv","metabolomics_feature_definitions.rds",
            "real_sacurine_12raw_qc_blank_feature_screens.csv",
            "real_sacurine_12raw_external_evidence.json")
if(any(!file.exists(file.path(root,required))))stop("Incomplete xcms baseline evidence")
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Refuse to overwrite comparator reference")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
area<-utils::read.csv(file.path(root,"metabolomics_peak_area.csv"),check.names=FALSE,stringsAsFactors=FALSE)
defs<-readRDS(file.path(root,"metabolomics_feature_definitions.rds"))
if(!all(c("feature_id","mzmed","rtmed")%in%names(defs))||
   anyDuplicated(defs$feature_id)||anyDuplicated(area$feature_id)||
   !setequal(as.character(defs$feature_id),as.character(area$feature_id)))
  stop("xcms feature identity/coordinate contract failed")
coords<-defs[match(area$feature_id,defs$feature_id),c("feature_id","mzmed","rtmed"),drop=FALSE]
if(anyNA(coords)||any(!is.finite(coords$mzmed))||any(coords$mzmed<=0)||
   any(!is.finite(coords$rtmed))||any(coords$rtmed<0))stop("Invalid xcms coordinates")
utils::write.csv(coords,file.path(out,"xcms_feature_coordinates.csv"),row.names=FALSE,na="")
file.copy(file.path(root,"real_sacurine_12raw_qc_blank_feature_screens.csv"),
          file.path(out,"xcms_qc_blank_screens.csv"))
file.copy(file.path(root,"real_sacurine_12raw_external_evidence.json"),
          file.path(out,"xcms_external_evidence.json"))
cat("PASS: complete xcms feature coordinates and QC screens exported for independent OpenMS comparison\n")
