#!/usr/bin/env Rscript
# External WHOLE 12-file Sacurine Zenodo v3 LC-MS source replay.
# Original raw mzML authenticated by the separate Zenodo v3 MD5 verifier.
# 6 real human urine samples, 3 pooled QC injections, 3 solvent blanks.
# Native R xcms feature tables and independent MS numerical QC checked.
#
# Unsupervised FEATURE extraction only. DO NOT test the reported age/sex/BMI
# associations here or make clinical, metabolite-ID or PK claims.
source("multiomics-engine/raw_ms_processing.R")
source("multiomics-engine/ms_quantitative_assurance.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)
  stop("Usage: Rscript real_sacurine_zenodo_12raw_xcms.R <verified raw directory> <empty evidence directory>")
root <- args[[1L]]
out <- args[[2L]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Refuse to overwrite existing real instrumental evidence.")
source_tsv <- file.path(root,"sampleMetadata_12samp_completed.tsv")
source_manifest <- file.path(root,"sacurine_v3_12_original_source_manifest.json")
if(!file.exists(source_tsv) || !file.exists(source_manifest))
  stop("Original data and source checksum/provenance report are required.")
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("Exact original source manifest parsing requires jsonlite.")
auth <- jsonlite::read_json(source_manifest,simplifyVector=TRUE)
if(!identical(auth$raw_source_verified,TRUE) ||
   !identical(auth$metadata_verified,TRUE) ||
   auth$raw_count!=12L || !identical(auth$source_injection_orders_verified,TRUE))
  stop("All 12 original published files and source chronology must be authenticated first.")
actual_tsv_md5 <- unname(tools::md5sum(source_tsv))
if(!identical(actual_tsv_md5,"f88044818627dee797c75d5a310f5c8e"))
  stop("Source sample metadata does not match the immutable publisher MD5.")
md <- utils::read.delim(source_tsv,check.names=FALSE,
                        stringsAsFactors=FALSE,na.strings=c("NA",""))
required <- c("sample_name","class","sampleType","injectionOrder","batch")
if(nrow(md)!=12L || any(!required %in% names(md)) ||
   anyDuplicated(names(md)) || anyDuplicated(md$sample_name) ||
   anyNA(md[,required]) || any(!is.finite(md$injectionOrder)) ||
   anyDuplicated(md$injectionOrder))
  stop("Incomplete/invalid authentic sample sheet.")
group_type <- c(sample="biological",QC="pooled_qc",blk="blank")
role <- unname(group_type[md$class])
if(anyNA(role) ||
   sum(role=="biological")!=6L ||
   sum(role=="pooled_qc")!=3L ||
   sum(role=="blank")!=3L ||
   any(md$sampleType != c(sample="sample",QC="pool",blk="blank")[md$class]) ||
   !all(md$batch=="ne1") ||
   !identical(range(md$injectionOrder),c(27L,185L)))
  stop("Source-defined sample roles, numbers, batch or true injection positions changed.")
files <- file.path(root,paste0(md$sample_name,".mzML"))
if(any(!file.exists(files)) || anyDuplicated(files) ||
   !setequal(basename(files),as.character(auth$files$name)))
  stop("The authenticated source file list and biological run manifest differ.")
# Only retain non-identifying sample and injection metadata in the output.
# Biological strata are deliberately an undifferentiated 'urine' label;
# do not infer sex, treatment groups or donor-level PK pairings.
manifest <- data.frame(
  sample_id=md$sample_name,assay_id=md$sample_name,
  sample_type=role,
  condition=ifelse(role=="biological","urine",NA_character_),
  injection_order=as.integer(md$injectionOrder),
  injection_order_verified=rep(TRUE,12L),
  batch=as.character(md$batch),stringsAsFactors=FALSE)
parameters <- list(
  peak_detection=list(ppm=25,peakwidth=c(10,60),snthresh=8,
                      prefilter=c(3,100)),
  retention_time=list(enabled=TRUE,binSize=.6),
  correspondence=list(bw=8,minFraction=.5,minSamples=1L,
                      binSize=.025,ppm=15))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
res <- run_xcms_raw_pipeline(
  files=files,sample_data=manifest,output_dir=out,
  parameters=parameters,chunk_size=2L)
stopifnot(identical(res$summary$status,"ok"),
  res$summary$samples==12L,
  res$summary$chromatographic_peaks>0L,
  res$summary$features>0L,
  res$summary$blank_samples_excluded_from_feature_definition==3L,
  ncol(res$matrix)==12L,
  nrow(res$matrix)==res$summary$features,
  !anyDuplicated(rownames(res$matrix)),
  identical(as.character(res$sample_data$assay_id),md$sample_name),
  identical(as.character(res$sample_data$raw_file),basename(files)),
  identical(as.character(res$sample_data$sample_type),role))
native <- t(res$matrix)
rownames(native) <- manifest$assay_id
stopifnot(!any(is.infinite(native)),all(native>=0,na.rm=TRUE),
  nrow(native)==12L,ncol(native)==res$summary$features,
  !anyDuplicated(colnames(native)),all(nzchar(colnames(native))))
csv <- utils::read.csv(file.path(out,"metabolomics_peak_area.csv"),
                       check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(csv)==ncol(native),
  setequal(csv$feature_id,colnames(native)),
  !anyDuplicated(csv$feature_id))
saved <- as.matrix(csv[,-1L,drop=FALSE])
storage.mode(saved) <- "double"
stopifnot(isTRUE(all.equal(unname(saved),
     unname(t(native)[csv$feature_id,,drop=FALSE]),
     tolerance=1e-8,check.attributes=FALSE)))
def <- readRDS(file.path(out,"metabolomics_feature_definitions.rds"))
stopifnot(nrow(def)==ncol(native))

# Independently recompute each instrumental screen from the original spectra
# output, rather than just comparing two filenames produced by our own code.
qa <- ms_quantitative_assurance(native,manifest)
coverage <- qa$coverage
stopifnot(coverage$biological_assays==6L,
          coverage$pooled_qc_assays==3L,
          coverage$blank_assays==3L,
          isTRUE(coverage$pooled_qc_precision_evaluable),
          isTRUE(coverage$blank_background_evaluable),
          isTRUE(coverage$pooled_qc_injection_order_verified),
          !coverage$pooled_qc_drift_screen_evaluable,
          !coverage$pooled_qc_batch_screen_evaluable,
          !coverage$missingness_group_imbalance_evaluable,
          !qa$scientific_certified,
          !qa$feature_exclusions_automatically_applied)
bio <- which(role=="biological")
qc <- which(role=="pooled_qc")
blank <- which(role=="blank")
independent_rsd <- vapply(seq_len(ncol(native)),function(j) {
  z <- native[qc,j]
  if(sum(is.finite(z))!=3L || mean(z,na.rm=TRUE)<=0) return(NA_real_)
  stats::sd(z)/mean(z)
},numeric(1))
independent_blank <- vapply(seq_len(ncol(native)),function(j) {
  b <- native[blank,j]
  x <- native[bio,j]
  if(sum(is.finite(b))<2L || sum(is.finite(x))<2L ||
     stats::median(x,na.rm=TRUE)<=0) return(NA_real_)
  stats::median(b,na.rm=TRUE)/stats::median(x,na.rm=TRUE)
},numeric(1))
eq_finite_and_missing <- function(a,b,tolerance=1e-10) {
  if(length(a)!=length(b)||any(is.na(a)!=is.na(b)))return(FALSE)
  ix <- which(!is.na(a))
  !length(ix) || isTRUE(all.equal(unname(a[ix]),unname(b[ix]),
                                 tolerance=tolerance))
}
f <- qa$features
stopifnot(identical(as.character(f$feature_id),colnames(native)),
  eq_finite_and_missing(independent_rsd,f$pooled_qc_rsd),
  eq_finite_and_missing(independent_blank,f$blank_over_biological_ratio),
  all(!is.finite(f$pooled_qc_order_spearman)),
  all(f$pooled_qc_drift_status=="too_few_pooled_qc_injections"))
utils::write.csv(f,file.path(out,"real_sacurine_12raw_qc_blank_feature_screens.csv"),
                 row.names=FALSE,na="")
jsonlite::write_json(list(
  schema="multiomics.sacurine-v3-twelve-real-raw-xcms-qa.v1",
  record="https://zenodo.org/records/10130758",
  doi="10.5281/zenodo.10130758",
  source_sha=Sys.getenv("QUALIFIED_SOURCE_SHA",
                        Sys.getenv("GITHUB_SHA","local_unverified")),
  source_published_md5_verified=TRUE,
  raw_spectra_md5_verified=TRUE,
  real_biological=6L,real_pooled_qc=3L,real_blanks=3L,
  all_original_injection_orders_verified=TRUE,
  original_acquisition_positions=as.integer(md$injectionOrder),
  native_chromatographic_peaks=res$summary$chromatographic_peaks,
  native_grouped_features=res$summary$features,
  numeric_export_replay_equal=TRUE,
  qc_rsd_features_evaluable=sum(is.finite(independent_rsd)),
  qc_rsd_gt30pct=sum(is.finite(independent_rsd)&independent_rsd>.30),
  blank_ratio_features_evaluable=sum(is.finite(independent_blank)),
  blank_ratio_ge20pct=sum(is.finite(independent_blank)&independent_blank>=.20),
  original_injection_order_is_sparse=TRUE,
  qc_injection_count_too_low_for_drift_screen=TRUE,
  biological_group_test_performed=FALSE,
  feature_identification_validated=FALSE,
  omics_pk_subject_linkage_validated=FALSE,
  clinical_certification=FALSE,
  scientific_certified=FALSE,
  limitation=paste("Authentic raw-data compatibility and descriptive feature-level",
    "QC/blank metrics only. Three pooled QC injections at sparse genuine",
    "injection positions cannot validate drift correction; no analytical",
    "identification, biological hypothesis test, FDR, metabolomic",
    "biomarker or linked PK outcome has been certified.")
),file.path(out,"real_sacurine_12raw_external_evidence.json"),
  pretty=TRUE,auto_unbox=TRUE,null="null")
cat("PASS: 12 genuine LC-MS raw files xcms processed:",
    res$summary$chromatographic_peaks,"chromatographic peaks;",
    res$summary$features,"features;",
    sum(is.finite(independent_rsd)),"RSDs;",
    sum(is.finite(independent_blank)),"blank/urine ratios;",
    "NO biological/FDR/clinical certification.\n")
