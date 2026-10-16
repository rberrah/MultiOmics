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

# xcms returns genuine peak index list-columns. Verify that CSV serialization
# preserves exact structured membership without collapsing to "1,2" strings.
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("jsonlite required for authentic featureDefinitions CSV test")
fixture<-data.frame(feature_id=c("FT1","FT2","FT3"),
                    mzmed=c(150.2,175.5,200.7),
                    stringsAsFactors=FALSE)
fixture$peakidx<-I(list(as.integer(c(1,3,9)),integer(0),as.integer(7)))
flat<-raw_ms_feature_definitions_csv(fixture)
stopifnot(is.character(flat$peakidx),nrow(flat)==3L,
          identical(flat$feature_id,fixture$feature_id),
          isTRUE(all.equal(flat$mzmed,fixture$mzmed)),
          identical(as.integer(jsonlite::fromJSON(flat$peakidx[[1L]])),
                    as.integer(c(1,3,9))),
          length(jsonlite::fromJSON(flat$peakidx[[2L]]))==0L,
          identical(as.integer(jsonlite::fromJSON(flat$peakidx[[3L]])),
                    as.integer(7)))
csv<-tempfile(fileext=".csv")
utils::write.csv(flat,csv,row.names=FALSE,na="")
reread<-utils::read.csv(csv,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(identical(reread$peakidx,flat$peakidx),
          identical(reread$feature_id,fixture$feature_id))
unlink(csv)
cat("multiomics raw MS contract: PASS\n")
