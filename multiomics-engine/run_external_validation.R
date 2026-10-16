#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: Rscript multiomics-engine/run_external_validation.R",
      "<predictions.csv> <config.json> <output_dir>"
    ),
    call. = FALSE
  )
}

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite first.", call. = FALSE)
source(file.path("multiomics-engine", "external_validation.R"))

predictions_path <- args[[1]]
config_path <- args[[2]]
output_dir <- args[[3]]
if (!file.exists(predictions_path)) stop("Predictions CSV was not found.", call. = FALSE)
if (!file.exists(config_path)) stop("Validation config JSON was not found.", call. = FALSE)

data <- utils::read.csv(predictions_path, check.names = FALSE, stringsAsFactors = FALSE)
config <- jsonlite::read_json(config_path, simplifyVector = TRUE)

allowed <- c(
  "outcome_type", "prediction_column", "outcome_column", "prediction_kind",
  "time_column", "event_column", "probability_columns", "independent_cohort",
  "cohort_label", "bootstrap_repetitions", "seed"
)
unknown <- setdiff(names(config), allowed)
if (length(unknown)) stop(paste("Unknown config field(s):", paste(unknown, collapse = ", ")), call. = FALSE)
if (is.null(config$outcome_type)) stop("config.outcome_type is required.", call. = FALSE)

arguments <- c(list(data = data), config)
result <- do.call(validate_external_predictions, arguments)
result$input <- list(
  predictions_file = basename(predictions_path),
  rows = nrow(data),
  columns = names(data)
)
result$environment <- list(
  r = R.version.string,
  platform = R.version$platform,
  jsonlite = as.character(utils::packageVersion("jsonlite"))
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  result,
  file.path(output_dir, "external_validation_report.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  na = "null"
)

summary_rows <- data.frame(
  field = c("status", "evaluation_status", "cohort", "outcome_type", "primary_metric"),
  value = c(
    result$status,
    result$evaluation_status,
    result$cohort %||% "",
    result$outcome_type,
    result$primary_metric
  ),
  stringsAsFactors = FALSE
)
metric_rows <- if (is.list(result$metrics)) {
  numeric_metrics <- result$metrics[vapply(result$metrics, function(x) length(x) == 1L && (is.numeric(x) || is.integer(x)), logical(1))]
  if (length(numeric_metrics)) data.frame(field = paste0("metric.", names(numeric_metrics)), value = unlist(numeric_metrics), stringsAsFactors = FALSE) else data.frame()
} else data.frame()
utils::write.csv(rbind(summary_rows, metric_rows), file.path(output_dir, "external_validation_summary.csv"), row.names = FALSE)

cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"), "\n")
