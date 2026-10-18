#!/usr/bin/env Rscript
# One PREDECLARED perturbation, same authentic 12 Zenodo v3 raw injections.
# This is parameter-sensitivity evidence, NOT a separate instrument pipeline.
# Baseline: ppm25 / snthresh8; perturbation: ppm20 / snthresh10.
# All other processing parameters, samples, and source are identical.
source("multiomics-engine/raw_ms_processing.R")
source("multiomics-engine/ms_quantitative_assurance.R")
source("multiomics-engine/ms_feature_sensitivity.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)
  stop("Usage: real_sacurine_xcms_sensitivity.R <verified raw dir> <complete baseline fit dir> <empty sensitivity dir>")
root <- args[[1L]]
baseline <- args[[2L]]
out <- args[[3L]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Refuse overwriting prior parameter-sensitivity evidence.")
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("Native R evidence requires jsonlite.")
required <- c("raw_ms_sample_sheet.csv","metabolomics_peak_area.csv",
              "metabolomics_feature_definitions.rds",
              "real_sacurine_12raw_external_evidence.json")
if(any(!file.exists(file.path(baseline,required))))
  stop("Complete genuine 12-raw baseline xcms experiment not found.")
source_manifest <- jsonlite::read_json(
  file.path(root,"sacurine_v3_12_original_source_manifest.json"),
  simplifyVector=TRUE)
if(!identical(source_manifest$raw_source_verified,TRUE) ||
   !identical(source_manifest$metadata_verified,TRUE) ||
   source_manifest$raw_count!=12L)
  stop("Parameter sensitivity requires an authenticated twelve-file source.")
md5 <- unname(tools::md5sum(file.path(root,"sampleMetadata_12samp_completed.tsv")))
if(!identical(md5,"f88044818627dee797c75d5a310f5c8e"))
  stop("Original Zenodo v3 sample metadata changed.")
sheet <- utils::read.csv(file.path(baseline,"raw_ms_sample_sheet.csv"),
                         stringsAsFactors=FALSE,check.names=FALSE)
orig <- utils::read.delim(file.path(root,"sampleMetadata_12samp_completed.tsv"),
                          stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(sheet)!=12L || nrow(orig)!=12L ||
   !identical(as.character(sheet$assay_id),as.character(orig$sample_name)) ||
   !identical(as.integer(sheet$injection_order),as.integer(orig$injectionOrder)) ||
   anyDuplicated(sheet$assay_id) ||
   !setequal(basename(file.path(root,paste0(sheet$sample_id,".mzML"))),
             as.character(source_manifest$files$name)) ||
   !identical(as.character(sheet$sample_type),
       unname(c(sample="biological",QC="pooled_qc",blk="blank")[orig$class])) ||
   any(!sheet$injection_order_verified))
  stop("Different subjects, original source order or study design: abort paired comparison.")
baseline_qc <- jsonlite::read_json(file.path(baseline,
                           "real_sacurine_12raw_external_evidence.json"),
                           simplifyVector=TRUE)
if(baseline_qc$native_grouped_features<1L ||
   baseline_qc$real_biological!=6L ||
   baseline_qc$real_pooled_qc!=3L ||
   baseline_qc$real_blanks!=3L ||
   !identical(baseline_qc$source_sha,
              Sys.getenv("QUALIFIED_SOURCE_SHA",Sys.getenv("GITHUB_SHA","local_unverified"))) ||
   !identical(baseline_qc$raw_spectra_md5_verified,TRUE))
  stop("Baseline R study integrity or code SHA differs from source run.")
paths <- file.path(root,paste0(sheet$sample_id,".mzML"))
if(any(!file.exists(paths)))stop("Paired study lost original chromatograms.")

baseline_area <- utils::read.csv(
  file.path(baseline,"metabolomics_peak_area.csv"),
  check.names=FALSE,stringsAsFactors=FALSE)
if(ncol(baseline_area)!=13L || anyDuplicated(baseline_area$feature_id))
  stop("Unusable baseline xcms feature table.")
baseline_defs <- readRDS(file.path(baseline,"metabolomics_feature_definitions.rds"))
reference <- ms_feature_coordinates(baseline_defs,baseline_area$feature_id)
# Convert feature-by-assay to assay-by-feature; same source *row* order.
ref <- t(as.matrix(baseline_area[,-1L,drop=FALSE]))
storage.mode(ref) <- "double"
rownames(ref) <- sheet$assay_id
colnames(ref) <- baseline_area$feature_id

# ONE prespecified, combined sensitivity variant. No tuning after seeing QA.
params <- list(
  peak_detection=list(ppm=20,peakwidth=c(10,60),
                      snthresh=10,prefilter=c(3,100)),
  retention_time=list(enabled=TRUE,binSize=.6),
  correspondence=list(bw=8,minFraction=.5,minSamples=1L,
                      binSize=.025,ppm=15))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
variant_dir <- file.path(out,"variant-ppm20-snthresh10")
variant <- run_xcms_raw_pipeline(
  files=paths,sample_data=sheet,output_dir=variant_dir,
  parameters=params,chunk_size=2L)
stopifnot(variant$summary$status=="ok",
          variant$summary$samples==12L,
          variant$summary$blank_samples_excluded_from_feature_definition==3L,
          variant$summary$features>0L,
          ncol(variant$matrix)==12L)
# Source/assay order must be identical; row IDs are NOT presumed stable.
stopifnot(identical(as.character(variant$sample_data$assay_id),
                    as.character(sheet$assay_id)),
          identical(as.character(variant$sample_data$sample_type),
                    as.character(sheet$sample_type)),
          identical(as.integer(variant$sample_data$injection_order),
                    as.integer(sheet$injection_order)))
alternate <- ms_feature_coordinates(
  readRDS(file.path(variant_dir,"metabolomics_feature_definitions.rds")),
  rownames(variant$matrix))
match <- ms_reciprocal_unique_matches(reference,alternate,
        ppm_tolerance=20,rt_seconds_tolerance=30)
stopifnot(nrow(match)==nrow(reference),
  !anyDuplicated(match$reference_feature),
  all(match$match_status %in% c("reciprocal_unique",
        "no_coordinate_candidate","ambiguous_multiple_candidates",
        "not_reciprocally_unique")))
alt <- t(variant$matrix)
rownames(alt) <- as.character(sheet$assay_id)
qa_ref <- ms_quantitative_assurance(ref,sheet)
qa_alt <- ms_quantitative_assurance(alt,sheet)
qc <- ms_paired_feature_sensitivity(match,qa_ref$features,qa_alt$features)
# No feature excluded; signed coordinate deltas and missing matches archived.
utils::write.csv(match,file.path(out,
   "xcms_prespecified_parameter_reciprocal_feature_matches.csv"),
   row.names=FALSE,na="")
summary <- list(
  schema="multiomics.sacurine-v3-xcms-peak-detection-sensitivity.v1",
  source="https://zenodo.org/records/10130758",
  source_sha=baseline_qc$source_sha,
  source_metadata_md5=md5,
  source_raw_checksums_verified=TRUE,
  twelve_identical_original_acquisitions=TRUE,
  original_injection_order_preserved=TRUE,
  baseline=list(ppm=25,snthresh=8,peakwidth_seconds=c(10,60),
                original_grouped_features=nrow(reference)),
  variant=list(ppm=20,snthresh=10,peakwidth_seconds=c(10,60),
               original_grouped_features=nrow(alternate)),
  xcms_alignment_and_grouping_parameters_unchanged=TRUE,
  matching_policy=list(
    one_to_one="reciprocal UNIQUE only, by native median feature coordinates",
    m_z_tolerance_ppm=20,retention_time_tolerance_seconds=30,
    unmatched_and_ambiguous_features_are_not_counted_as_stable=TRUE),
  coordinate_match_status_counts=as.list(table(factor(match$match_status,
    levels=c("reciprocal_unique","no_coordinate_candidate",
             "ambiguous_multiple_candidates","not_reciprocally_unique")))),
  matched_feature_count=qc$pairs,
  fraction_reference_features_reciprocally_matched=qc$pairs/nrow(reference),
  fraction_variant_features_reciprocally_matched=qc$pairs/nrow(alternate),
  matched_both_qc_evaluable=qc$both_qc_evaluable,
  matched_both_blank_evaluable=qc$both_blank_evaluable,
  paired_qc_rsd_30pct_screen_agreement=qc$qc_screen_agreement,
  paired_blank_20pct_screen_agreement=qc$blank_screen_agreement,
  independent_analytical_pipeline_validated=FALSE,
  parameter_causal_attribution_possible=FALSE,
  parameter_optimisation_performed=FALSE,
  feature_molecular_identity_confirmed=FALSE,
  qc_screen_is_scientific_certification=FALSE,
  biological_inference_performed=FALSE,
  omics_PK_pairing_available=FALSE,
  review_status="scientific_sensitivity_review_required",
  limitations=paste("This prespecified combined peak-picking perturbation",
    "does not identify the effect of ppm versus S/N separately.",
    "One-to-one matched m/z and RT features are not chemical identities;",
    "ambiguous/unmatched peaks remain visible. No inference of",
    "biomarker performance, controlled FDR or clinical PK associations.")
)
jsonlite::write_json(summary,file.path(out,
  "xcms_prespecified_sensitivity_summary.json"),
  pretty=TRUE,auto_unbox=TRUE,null="null")
cat("PASS genuine source sensitivity, baseline",nrow(reference),"variant",
    nrow(alternate),"unique reciprocal coordinate matches",qc$pairs,
    "INFERENCE NOT CERTIFIED\n")
