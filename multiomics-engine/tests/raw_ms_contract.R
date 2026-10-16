#!/usr/bin/env Rscript
source(file.path("multiomics-engine", "raw_ms_processing.R"))

stopifnot(identical(raw_ms_format_status("sample.mzML")$status, "supported_open"))
stopifnot(identical(raw_ms_format_status("sample.mzXML")$status, "supported_open"))
stopifnot(identical(raw_ms_format_status("sample.CDF")$status, "supported_open"))
stopifnot(identical(raw_ms_format_status("thermo.raw")$status, "vendor_conversion_required"))
stopifnot(identical(raw_ms_format_status("waters.RAW")$status, "vendor_conversion_required"))
stopifnot(identical(raw_ms_format_status("agilent.d")$status, "vendor_conversion_required"))
stopifnot(identical(raw_ms_format_status("notes.txt")$status, "unsupported"))

params <- raw_ms_default_parameters()
stopifnot(isTRUE(raw_ms_validate_parameters(params)))
invalid <- params
invalid$peak_detection$ppm <- 0
stopifnot(inherits(try(raw_ms_validate_parameters(invalid), silent = TRUE), "try-error"))

sample_data <- data.frame(
  condition = c("A", "A", "B", "", ""),
  sample_type = c("biological", "biological", "biological", "pooled_qc", "blank"),
  stringsAsFactors = FALSE
)
groups <- raw_ms_sample_groups(sample_data)
stopifnot(identical(groups[1:4], c("A", "A", "B", "pooled_qc")))
stopifnot(is.na(groups[[5]]))

plan <- proteowizard_msconvert_plan("thermo.raw", "converted", centroid = TRUE, zlib = TRUE)
stopifnot(identical(plan$executable, "msconvert"))
stopifnot("--mzML" %in% plan$args)
stopifnot("--zlib" %in% plan$args)
stopifnot("peakPicking true 1-" %in% plan$args)

# Vendor files must be rejected by the xcms input contract before any package
# or filesystem-dependent processing is attempted.
err <- try(
  raw_ms_validate_inputs(
    c("thermo.raw"),
    data.frame(sample_id = "S1", stringsAsFactors = FALSE),
    require_existing = FALSE
  ),
  silent = TRUE
)
stopifnot(inherits(err, "try-error"))

cat("multiomics raw MS contract: PASS\n")
