# Standards-based LC-MS/GC-MS preprocessing for the multi-omics backend.
#
# Scientific boundary:
# - open exchange formats (mzML/mzXML/netCDF) are processed directly;
# - proprietary vendor formats are not parsed by this project. They should be
#   converted locally with ProteoWizard/msconvert (or the vendor's validated
#   exporter) before entering this pipeline;
# - default xcms parameters are reproducible starting values, not instrument-
#   independent optimal settings. Confirmatory analyses should prespecify or
#   justify instrument/method-specific parameters.

raw_ms_require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for raw MS processing.", pkg), call. = FALSE)
  }
}

raw_ms_supported_extensions <- function() c(".mzml", ".mzxml", ".cdf", ".netcdf")
raw_ms_vendor_extensions <- function() c(".raw", ".wiff", ".wiff2", ".d", ".baf", ".tdf", ".tdf_bin", ".yep", ".fid")

raw_ms_extension <- function(path) {
  value <- tolower(basename(as.character(path)))
  # Directory-like vendor formats such as *.d need an explicit suffix check.
  hits <- c(raw_ms_vendor_extensions(), raw_ms_supported_extensions())
  match <- hits[endsWith(value, hits)]
  if (length(match)) match[[which.max(nchar(match))]] else tolower(tools::file_ext(value)) |> paste0(".")
}

raw_ms_format_status <- function(path) {
  ext <- raw_ms_extension(path)
  if (ext %in% raw_ms_supported_extensions()) {
    return(list(status = "supported_open", extension = ext, action = "process_directly"))
  }
  if (ext %in% raw_ms_vendor_extensions()) {
    return(list(
      status = "vendor_conversion_required",
      extension = ext,
      action = "convert_to_mzML",
      note = "Convert locally with ProteoWizard/msconvert or a validated vendor exporter before xcms preprocessing."
    ))
  }
  list(status = "unsupported", extension = ext, action = "review_format")
}

raw_ms_default_parameters <- function() {
  list(
    peak_detection = list(
      method = "centWave",
      ppm = 15,
      peakwidth = c(5, 30),
      snthresh = 10,
      prefilter = c(3, 100),
      mzCenterFun = "wMean",
      integrate = 1L,
      mzdiff = -0.001,
      noise = 0
    ),
    retention_time = list(
      method = "obiwarp",
      enabled = TRUE,
      binSize = 0.6,
      response = 1L,
      distFun = "cor_opt",
      gapInit = numeric(),
      gapExtend = numeric(),
      factorDiag = 2,
      factorGap = 1,
      localAlignment = FALSE
    ),
    correspondence = list(
      method = "peak_density",
      bw = 5,
      minFraction = 0.5,
      minSamples = 1L,
      binSize = 0.025,
      ppm = 10
    ),
    gap_filling = list(
      enabled = TRUE,
      method = "ChromPeakAreaParam"
    ),
    quantification = list(
      value = "into",
      multiple_peak_method = "medret"
    )
  )
}

raw_ms_merge_parameters <- function(defaults, overrides = list()) {
  if (!length(overrides)) return(defaults)
  out <- defaults
  for (name in names(overrides)) {
    if (is.list(out[[name]]) && is.list(overrides[[name]])) {
      out[[name]] <- utils::modifyList(out[[name]], overrides[[name]])
    } else {
      out[[name]] <- overrides[[name]]
    }
  }
  out
}

raw_ms_validate_parameters <- function(params) {
  pd <- params$peak_detection
  corr <- params$correspondence
  if (!is.numeric(pd$ppm) || length(pd$ppm) != 1L || !is.finite(pd$ppm) || pd$ppm <= 0) {
    stop("centWave ppm must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(pd$peakwidth) || length(pd$peakwidth) != 2L || any(!is.finite(pd$peakwidth)) ||
      pd$peakwidth[[1]] <= 0 || pd$peakwidth[[2]] <= pd$peakwidth[[1]]) {
    stop("centWave peakwidth must be c(min,max) with 0 < min < max.", call. = FALSE)
  }
  if (!is.numeric(pd$snthresh) || length(pd$snthresh) != 1L || pd$snthresh < 0) {
    stop("centWave snthresh must be non-negative.", call. = FALSE)
  }
  if (!is.numeric(corr$minFraction) || length(corr$minFraction) != 1L ||
      corr$minFraction <= 0 || corr$minFraction > 1) {
    stop("Peak-density minFraction must be in (0,1].", call. = FALSE)
  }
  invisible(TRUE)
}

raw_ms_canonical_sample_type <- function(value) {
  key <- tolower(trimws(as.character(value %||% "")))
  if (key %in% c("blank", "blanc", "solvent_blank", "solvent blank")) return("blank")
  if (key %in% c("qc", "pooled_qc", "pooled qc", "quality_control", "quality control", "pool")) return("pooled_qc")
  "biological"
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x))) y else x

raw_ms_sample_groups <- function(sample_data, condition_column = "condition", sample_type_column = "sample_type") {
  sample_data <- as.data.frame(sample_data, stringsAsFactors = FALSE)
  n <- nrow(sample_data)
  if (!n) return(character())
  sample_types <- if (sample_type_column %in% names(sample_data)) {
    vapply(sample_data[[sample_type_column]], raw_ms_canonical_sample_type, character(1))
  } else rep("biological", n)
  condition <- if (condition_column %in% names(sample_data)) trimws(as.character(sample_data[[condition_column]])) else rep("", n)
  group <- ifelse(nzchar(condition), condition, "all_biological")
  group[sample_types == "pooled_qc"] <- "pooled_qc"
  # Blanks remain available for downstream blank-ratio QC, but do not define
  # chromatographic features by themselves during correspondence analysis.
  group[sample_types == "blank"] <- NA_character_
  group
}

raw_ms_validate_inputs <- function(files, sample_data, require_existing = TRUE) {
  files <- as.character(files)
  sample_data <- as.data.frame(sample_data, stringsAsFactors = FALSE)
  if (!length(files)) stop("No raw MS files were supplied.", call. = FALSE)
  if (nrow(sample_data) != length(files)) {
    stop("sample_data must contain exactly one row per raw MS file, in the same order.", call. = FALSE)
  }
  statuses <- lapply(files, raw_ms_format_status)
  vendor <- which(vapply(statuses, function(x) identical(x$status, "vendor_conversion_required"), logical(1)))
  unsupported <- which(vapply(statuses, function(x) identical(x$status, "unsupported"), logical(1)))
  if (length(vendor)) {
    stop(
      paste0(
        "Proprietary vendor MS files require conversion to mzML before preprocessing: ",
        paste(basename(files[vendor]), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(unsupported)) {
    stop(paste0("Unsupported raw MS format: ", paste(basename(files[unsupported]), collapse = ", ")), call. = FALSE)
  }
  if (isTRUE(require_existing) && any(!file.exists(files))) {
    stop(paste0("Raw MS file(s) not found: ", paste(basename(files[!file.exists(files)]), collapse = ", ")), call. = FALSE)
  }
  invisible(statuses)
}

proteowizard_msconvert_plan <- function(input, output_dir, centroid = TRUE, zlib = TRUE) {
  status <- raw_ms_format_status(input)
  if (!status$status %in% c("vendor_conversion_required", "supported_open")) {
    stop("Input format is not recognised as a supported MS source.", call. = FALSE)
  }
  args <- c(as.character(input), "--mzML", "--outdir", as.character(output_dir))
  if (isTRUE(zlib)) args <- c(args, "--zlib")
  if (isTRUE(centroid)) args <- c(args, "--filter", "peakPicking true 1-")
  list(
    executable = "msconvert",
    args = args,
    input_status = status$status,
    note = "Inspect centroiding requirements for the acquisition method before using peakPicking. Vendor readers may require Windows/vendor libraries."
  )
}

convert_vendor_to_mzml <- function(input, output_dir, centroid = TRUE, zlib = TRUE, execute = FALSE) {
  plan <- proteowizard_msconvert_plan(input, output_dir, centroid = centroid, zlib = zlib)
  if (!isTRUE(execute)) return(plan)
  exe <- Sys.which(plan$executable)
  if (!nzchar(exe)) stop("ProteoWizard msconvert was not found on PATH.", call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(exe, args = plan$args, stdout = TRUE, stderr = TRUE)
  attr(status, "plan") <- plan
  status
}

raw_ms_package_manifest <- function() {
  pkgs <- c("xcms", "MsExperiment", "Spectra", "mzR", "BiocParallel")
  versions <- lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
    as.character(utils::packageVersion(pkg))
  })
  names(versions) <- pkgs
  list(r = R.version.string, platform = R.version$platform, packages = versions)
}

run_xcms_raw_pipeline <- function(
  files,
  sample_data,
  output_dir,
  parameters = list(),
  condition_column = "condition",
  sample_type_column = "sample_type",
  chunk_size = 2L
) {
  raw_ms_require_namespace("xcms")
  raw_ms_require_namespace("MsExperiment")
  raw_ms_require_namespace("Spectra")
  raw_ms_require_namespace("mzR")
  raw_ms_validate_inputs(files, sample_data, require_existing = TRUE)

  params <- raw_ms_merge_parameters(raw_ms_default_parameters(), parameters)
  raw_ms_validate_parameters(params)
  sample_data <- as.data.frame(sample_data, stringsAsFactors = FALSE)
  sample_data$raw_file <- basename(files)
  groups <- raw_ms_sample_groups(sample_data, condition_column, sample_type_column)
  sample_data$xcms_feature_group <- groups

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mse <- MsExperiment::readMsExperiment(
    spectraFiles = files,
    sampleData = sample_data
  )

  pd <- params$peak_detection
  peak_param <- xcms::CentWaveParam(
    ppm = pd$ppm,
    peakwidth = as.numeric(pd$peakwidth),
    snthresh = pd$snthresh,
    prefilter = as.numeric(pd$prefilter),
    mzCenterFun = pd$mzCenterFun,
    integrate = as.integer(pd$integrate),
    mzdiff = pd$mzdiff,
    noise = pd$noise
  )
  xdata <- xcms::findChromPeaks(mse, param = peak_param)

  rt <- params$retention_time
  if (isTRUE(rt$enabled) && length(files) >= 2L) {
    align_param <- xcms::ObiwarpParam(
      binSize = rt$binSize,
      response = as.integer(rt$response),
      distFun = rt$distFun,
      gapInit = rt$gapInit,
      gapExtend = rt$gapExtend,
      factorDiag = rt$factorDiag,
      factorGap = rt$factorGap,
      localAlignment = isTRUE(rt$localAlignment)
    )
    xdata <- xcms::adjustRtime(xdata, param = align_param, chunkSize = as.integer(max(1L, chunk_size)))
  }

  corr <- params$correspondence
  correspondence_param <- xcms::PeakDensityParam(
    sampleGroups = groups,
    bw = corr$bw,
    minFraction = corr$minFraction,
    minSamples = as.integer(corr$minSamples),
    binSize = corr$binSize,
    ppm = corr$ppm
  )
  xdata <- xcms::groupChromPeaks(xdata, param = correspondence_param)

  if (isTRUE(params$gap_filling$enabled)) {
    xdata <- xcms::fillChromPeaks(xdata, param = xcms::ChromPeakAreaParam())
  }

  values <- xcms::featureValues(
    xdata,
    value = params$quantification$value,
    method = params$quantification$multiple_peak_method,
    filled = TRUE,
    missing = NA_real_
  )
  feature_defs <- as.data.frame(xcms::featureDefinitions(xdata), stringsAsFactors = FALSE)
  feature_ids <- rownames(values)
  if (is.null(feature_ids) || any(!nzchar(feature_ids))) feature_ids <- paste0("FT", seq_len(nrow(values)))
  rownames(values) <- feature_ids
  if (nrow(feature_defs) == length(feature_ids)) {
    feature_defs$feature_id <- rownames(feature_defs)
    if (any(!nzchar(feature_defs$feature_id))) feature_defs$feature_id <- feature_ids
  } else {
    feature_defs <- data.frame(feature_id = feature_ids, stringsAsFactors = FALSE)
  }

  intensity_table <- data.frame(feature_id = feature_ids, values, check.names = FALSE, row.names = NULL)
  utils::write.csv(intensity_table, file.path(output_dir, "metabolomics_peak_area.csv"), row.names = FALSE, na = "")
  utils::write.csv(feature_defs, file.path(output_dir, "metabolomics_feature_definitions.csv"), row.names = FALSE, na = "")
  utils::write.csv(sample_data, file.path(output_dir, "raw_ms_sample_sheet.csv"), row.names = FALSE, na = "")

  peaks <- xcms::chromPeaks(xdata)
  peak_counts <- integer(length(files))
  if (nrow(peaks) && "sample" %in% colnames(peaks)) {
    tab <- table(as.integer(peaks[, "sample"]))
    peak_counts[as.integer(names(tab))] <- as.integer(tab)
  }
  missing_fraction <- if (length(values)) mean(!is.finite(values)) else NA_real_
  filled_count <- NA_integer_
  peak_data <- try(xcms::chromPeakData(xdata), silent = TRUE)
  if (!inherits(peak_data, "try-error") && "is_filled" %in% colnames(peak_data)) {
    filled_count <- sum(as.logical(peak_data$is_filled), na.rm = TRUE)
  }

  summary <- list(
    status = "ok",
    method = "MsExperiment + xcms",
    files = lapply(seq_along(files), function(i) list(
      sample_index = i,
      file = basename(files[[i]]),
      format = raw_ms_extension(files[[i]]),
      chromatographic_peaks = unname(peak_counts[[i]])
    )),
    samples = length(files),
    features = nrow(values),
    chromatographic_peaks = nrow(peaks),
    filled_chromatographic_peaks = filled_count,
    matrix_missing_fraction = missing_fraction,
    blank_samples_excluded_from_feature_definition = sum(is.na(groups)),
    parameters = params,
    environment = raw_ms_package_manifest(),
    outputs = list(
      matrix = "metabolomics_peak_area.csv",
      feature_definitions = "metabolomics_feature_definitions.csv",
      sample_sheet = "raw_ms_sample_sheet.csv"
    ),
    interpretation = paste(
      "Peak picking, retention-time alignment, correspondence and gap filling create an untargeted feature table.",
      "They do not identify compounds. Parameter adequacy, pooled-QC behaviour, blank contamination and feature annotation must still be reviewed."
    )
  )

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(summary, file.path(output_dir, "raw_ms_processing_manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  }
  saveRDS(xdata, file.path(output_dir, "xcms_processed_experiment.rds"))
  invisible(list(experiment = xdata, matrix = values, features = feature_defs, sample_data = sample_data, summary = summary))
}
