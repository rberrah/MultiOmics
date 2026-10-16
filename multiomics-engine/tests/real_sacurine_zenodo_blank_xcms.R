#!/usr/bin/env Rscript
# Sacurine independent real open-format external chromatogram test.
# 2 independent biological urine runs + 2 *actual* solvent blank runs.
# No pooled-QC exists; never fabricate temporal drift or batch correction.
source("multiomics-engine/raw_ms_processing.R")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: Rscript real_sacurine_zenodo_blank_xcms.R <raw source dir> <evidence dir>")
root<-args[[1L]]
out<-args[[2L]]
types<-c("blank","blank","biological","biological")
names<-c("Blanc12","Blanc17","HU_neg_099","HU_neg_178")
files<-file.path(root,paste0(names,".mzXML"))
stopifnot(all(file.exists(files)),length(unique(files))==4L)
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
 stop("Refuse to overwrite existing evidence")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sample_data<-data.frame(
 sample_id=names,
 assay_id=names,
 sample_type=types,
 condition=c(NA_character_,NA_character_,"urine","urine"),
 injection_order=seq_along(names),  # arbitrary local key; NOT real acquisition time
 injection_order_verified=rep(FALSE,4L),
 stringsAsFactors=FALSE
)
parameters<-list(
 peak_detection=list(ppm=25,peakwidth=c(10,60),snthresh=8,prefilter=c(3,100)),
 retention_time=list(enabled=TRUE,binSize=.6),
 correspondence=list(bw=8,minFraction=.5,minSamples=1L,binSize=.025,ppm=15)
)
res<-run_xcms_raw_pipeline(files=files,sample_data=sample_data,
                           output_dir=out,parameters=parameters,chunk_size=1L)
stopifnot(identical(res$summary$status,"ok"),
 res$summary$samples==4L,
 res$summary$chromatographic_peaks>0L,
 res$summary$features>0L,
 res$summary$blank_samples_excluded_from_feature_definition==2L,
 ncol(res$matrix)==4L,nrow(res$matrix)==res$summary$features,
 !anyDuplicated(rownames(res$matrix)),
 identical(as.character(res$sample_data$sample_type),types))
intensities<-res$matrix
stopifnot(all(is.finite(intensities)|is.na(intensities)),
 all(intensities[is.finite(intensities)]>=0))
csv<-utils::read.csv(file.path(out,"metabolomics_peak_area.csv"),
                     check.names=FALSE,stringsAsFactors=FALSE)
defs<-utils::read.csv(file.path(out,"metabolomics_feature_definitions.csv"),
                     check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(csv)==nrow(intensities),nrow(defs)==nrow(intensities),
 setequal(csv$feature_id,rownames(intensities)),!anyDuplicated(csv$feature_id))
saved<-as.matrix(csv[,-1L,drop=FALSE])
storage.mode(saved)<-"double"
stopifnot(isTRUE(all.equal(unname(saved),unname(intensities[csv$feature_id,,drop=FALSE]),
                           tolerance=1e-8,check.attributes=FALSE)))
# This is an *independent genuine blank calculation*, not a clinical
# contamination false-positive estimate. No unverified batch or RT correction.
bio_ix<-which(types=="biological")
blank_ix<-which(types=="blank")
blank_ratios<-rep(NA_real_,nrow(intensities))
not_evaluable<-character(nrow(intensities))
for(i in seq_len(nrow(intensities))){
 b<-intensities[i,blank_ix]
 v<-intensities[i,bio_ix]
 if(sum(is.finite(b))<2L||sum(is.finite(v))<2L){
   not_evaluable[[i]]<-"insufficient_finite_blank_or_biological"
 }else if(stats::median(v)<=0){
   not_evaluable[[i]]<-"zero_biological_median"
 }else{
   blank_ratios[[i]]<-stats::median(b)/stats::median(v)
   not_evaluable[[i]]<-"observed_descriptive_blank_to_urine_ratio"
 }
}
stopifnot(all(blank_ratios[is.finite(blank_ratios)]>=0))
metrics<-data.frame(feature_id=rownames(intensities),
                    observed_blank_to_biological_median_ratio=blank_ratios,
                    status=not_evaluable,
                    stringsAsFactors=FALSE)
utils::write.csv(metrics,file.path(out,"true_blank_external_feature_ratios.csv"),
                 row.names=FALSE,na="")
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite required for real source reporting")
jsonlite::write_json(list(
 schema="multiomics.sacurine-toy-real-blank-evidence.v1",
 source="Zenodo 1346742 published LC-MS mzXML",
 source_sha=Sys.getenv("GITHUB_SHA","local"),
 real_chromatograms=4L,
 real_biological_urine=2L,
 real_solvent_blanks=2L,
 pooled_qc_injections=0L,
 chromatographic_peaks=res$summary$chromatographic_peaks,
 grouped_features=res$summary$features,
 evaluable_true_blank_ratios=sum(is.finite(blank_ratios)),
 analytical_blank_flag_threshold=.2,
 descriptive_features_blank_ratio_gte_20pct=sum(is.finite(blank_ratios)&blank_ratios>=.2),
 blank_feature_definitions_excluded_count=2L,
 output_alignment_checked=TRUE,
 laboratory_replicate_study_sufficient=FALSE,
 qc_precision_evaluable=FALSE,
 injection_order_verified=FALSE,
 injection_drift_evaluable=FALSE,
 batch_correction_performed=FALSE,
 synthetic_measurements=FALSE,
 biological_statistical_test_performed=FALSE,
 limitation="Only 2 urine and 2 blank real runs; raw true blank diagnostic, not proof of specificity, clinical significance, FDR or population performance."
),file.path(out,"sacurine_real_blank_evidence.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: external Sacurine raw mzXML chromatograms",res$summary$chromatographic_peaks,
    " peaks,",res$summary$features," grouped features,",
    sum(is.finite(blank_ratios))," real blank-ratio denominators\n")
