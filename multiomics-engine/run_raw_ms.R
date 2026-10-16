#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: Rscript multiomics-engine/run_raw_ms.R",
      "<raw_ms_manifest.csv> <output_dir> <parameters.json|default>",
      "\nManifest columns: raw_file plus optional sample_id, assay_id, condition, sample_type, injection_order."
    ),
    call. = FALSE
  )
}

source(file.path("multiomics-engine", "raw_ms_processing.R"))
manifest_path <- args[[1]]
output_dir <- args[[2]]
parameters_path <- args[[3]]

if (!file.exists(manifest_path)) stop("Raw-MS manifest CSV was not found.", call. = FALSE)
manifest <- utils::read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"raw_file" %in% names(manifest)) stop("Manifest must contain a raw_file column.", call. = FALSE)

base <- dirname(normalizePath(manifest_path, mustWork = TRUE))
files <- as.character(manifest$raw_file)
relative <- !grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", files)
files[relative] <- file.path(base, files[relative])

parameters <- list()
if (!identical(tolower(parameters_path), "default")) {
  if (!file.exists(parameters_path)) stop("Parameter JSON was not found.", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite to read parameter JSON.", call. = FALSE)
  parameters <- jsonlite::read_json(parameters_path, simplifyVector = TRUE)
}

result <- run_xcms_raw_pipeline(
  files = files,
  sample_data = manifest,
  output_dir = output_dir,
  parameters = parameters
)

cat(sprintf(
  "Raw MS preprocessing complete: %d files -> %d features\nMatrix: %s\n",
  result$summary$samples,
  result$summary$features,
  file.path(output_dir, result$summary$outputs$matrix)
))
