#!/usr/bin/env Rscript
source(file.path("multiomics-engine", "raw_ms_processing.R"))

required <- c("xcms", "MsExperiment", "Spectra", "mzR", "faahKO")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop(paste("Missing integration-test packages:", paste(missing, collapse = ", ")), call. = FALSE)

cdf_root <- system.file("cdf", package = "faahKO")
if (!nzchar(cdf_root)) stop("faahKO CDF directory was not found.", call. = FALSE)
files <- list.files(cdf_root, pattern = "\\.CDF$", recursive = TRUE, full.names = TRUE)
if (length(files) < 4L) stop("Expected at least four faahKO CDF files.", call. = FALSE)

# Keep the integration test representative but bounded: two KO + two WT files.
ko <- files[grepl("/KO/|\\\\KO\\\\", files, ignore.case = TRUE)][1:2]
wt <- files[grepl("/WT/|\\\\WT\\\\", files, ignore.case = TRUE)][1:2]
selected <- c(ko, wt)
if (any(is.na(selected)) || length(selected) != 4L) stop("Could not select two KO and two WT faahKO files.", call. = FALSE)

sample_data <- data.frame(
  sample_id = tools::file_path_sans_ext(basename(selected)),
  assay_id = tools::file_path_sans_ext(basename(selected)),
  condition = c("KO", "KO", "WT", "WT"),
  sample_type = "biological",
  injection_order = seq_along(selected),
  stringsAsFactors = FALSE
)

# faahKO CDF files are a restricted retention-time/mass subset. Use permissive
# but still explicit parameters to keep this as an API/pipeline integration
# test rather than an attempt to re-optimise the published biological study.
parameters <- list(
  peak_detection = list(
    ppm = 25,
    peakwidth = c(10, 60),
    snthresh = 8,
    prefilter = c(3, 100)
  ),
  retention_time = list(
    enabled = TRUE,
    binSize = 0.6
  ),
  correspondence = list(
    bw = 8,
    minFraction = 0.5,
    minSamples = 1L,
    binSize = 0.025,
    ppm = 15
  )
)

output_from_ci <- Sys.getenv("MULTIOMICS_RAW_MS_OUTPUT","")
out <- if(nzchar(output_from_ci)) output_from_ci else tempfile("pmx_raw_ms_")
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Raw MS evidence directory is not empty; refusing to reuse stale artifacts")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
res <- run_xcms_raw_pipeline(
  files = selected,
  sample_data = sample_data,
  output_dir = out,
  parameters = parameters,
  chunk_size = 1L
)

stopifnot(identical(res$summary$status, "ok"))
stopifnot(res$summary$samples == 4L)
stopifnot(res$summary$chromatographic_peaks > 0L)
stopifnot(res$summary$features > 0L)
stopifnot(ncol(res$matrix) == 4L)
stopifnot(nrow(res$matrix) == res$summary$features)
stopifnot(file.exists(file.path(out, "metabolomics_peak_area.csv")))
stopifnot(file.exists(file.path(out, "metabolomics_feature_definitions.csv")))
stopifnot(file.exists(file.path(out, "raw_ms_sample_sheet.csv")))
stopifnot(file.exists(file.path(out, "xcms_processed_experiment.rds")))

# Quantify the actual saved output rather than only checking file existence.
saved <- utils::read.csv(file.path(out,"metabolomics_peak_area.csv"),
                         check.names=FALSE,stringsAsFactors=FALSE)
features <- utils::read.csv(file.path(out,"metabolomics_feature_definitions.csv"),
                            check.names=FALSE,stringsAsFactors=FALSE)
manifest <- utils::read.csv(file.path(out,"raw_ms_sample_sheet.csv"),
                            check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(saved)==res$summary$features,
          nrow(saved)==nrow(features),
          nrow(manifest)==4L,
          "feature_id"%in%names(saved),
          "feature_id"%in%names(features),
          !anyNA(saved$feature_id),
          !anyDuplicated(saved$feature_id),
          setequal(saved$feature_id,features$feature_id),
          setequal(names(saved)[-1L],colnames(res$matrix)),
          identical(as.character(manifest$condition),sample_data$condition),
          sum(manifest$condition=="KO")==2L,
          sum(manifest$condition=="WT")==2L)
actual<-as.matrix(saved[,-1L,drop=FALSE])
storage.mode(actual)<-"double"
stopifnot(all(is.finite(actual)|is.na(actual)),
          all(actual[is.finite(actual)]>=0),
          isTRUE(all.equal(unname(actual),unname(res$matrix[saved$feature_id,,drop=FALSE]),
                           check.attributes=FALSE,tolerance=1e-8)),
          isTRUE(all.equal(mean(!is.finite(res$matrix)),
                           res$summary$matrix_missing_fraction,tolerance=1e-12)))
if(nzchar(output_from_ci)){
  package_names<-c("xcms","MsExperiment","Spectra","mzR","faahKO","BiocManager")
  installed_versions<-vapply(package_names,function(pkg)
    as.character(utils::packageVersion(pkg)),character(1))
  writeLines(c(
    paste("github_sha:",Sys.getenv("GITHUB_SHA","not_recorded")),
    paste("reference_r:",getRversion()),
    paste("bioconductor:",BiocManager::version()),
    paste("package_versions:",paste(names(installed_versions),installed_versions,
                                    sep="=",collapse="; ")),
    paste("raw_files:",paste(basename(selected),collapse="; ")),
    "2 KO + 2 WT public faahKO subset; open netCDF CDF files",
    "Prespecified integration/peak-grouping parameters are NOT instrument-independent.",
    "4 real files form an API/reproducibility test; no FDR, group-effect or validated-biomarker inference.",
    capture.output(sessionInfo())),
    file.path(out,"real_faahKO_provenance.txt"))
}


cat(sprintf(
  "multiomics raw MS xcms integration: PASS — %d chromatographic peaks, %d grouped features\n",
  res$summary$chromatographic_peaks,
  res$summary$features
))
