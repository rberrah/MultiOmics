#!/usr/bin/env Rscript
# Independent native-R replay of real GEO processed Cyp3a11 probe observations.
# Does not fit, project or validate any frozen Nutrimouse/PK model.
args <- commandArgs(trailingOnly=TRUE)
rank_native <- function(x, value) {
  (sum(x < value) + (sum(x == value)-1)/2)/(length(x)-1)
}
if (identical(args, "--self-test")) {
  stopifnot(abs(rank_native(c(1,2,2,3,4),2)-.375)<1e-14)
  stopifnot(identical(rank_native(c(4,7,9,10,15),9),
                      rank_native(c(112,121,127,130,145),127)))
  cat("PASS: independently implemented R midrank replay test\n")
  quit(save="no")
}
if (length(args)!=2L) stop("Usage: Rscript test_geo_native_probe_replay.R GEO_DIR BIOC_DIR")
source_dir <- args[[1L]]
audit_dir <- args[[2L]]
data_file <- file.path(source_dir, "GSE22871_series_matrix.txt.gz")
observations <- file.path(audit_dir, "native_external_probe_observations",
                          "native_Cyp3a11_probe_measurements_not_model_scores.csv")
report_file <- file.path(audit_dir, "native_external_probe_observations",
                         "external_native_measurement_audit.json")
if (!file.exists(data_file) || !file.exists(observations) ||
    !file.exists(report_file)) stop("Independent replay missing original external evidence")
if (!requireNamespace("jsonlite",quietly=TRUE)) stop("Missing JSON evidence parser")
report <- jsonlite::fromJSON(report_file,simplifyVector=TRUE)
stopifnot(identical(report$status,"EXTERNAL_NATIVE_MEASUREMENT_ONLY_PROJECTION_BLOCKED"),
          identical(report$external_biological_signature_scores_exported,0L),
          identical(report$n_independent_mice,30L),
          identical(report$n_original_measured_probes,45101L))
values <- utils::read.csv(observations,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(values)==30L*report$n_anchor_candidates,
          all(values$cross_platform_model_score=="NOT_COMPUTED"),
          length(unique(values$assay_id))==30L)
connection <- gzfile(data_file,open="rt")
on.exit(close(connection),add=TRUE)
lines <- readLines(connection,warn=FALSE,encoding="UTF-8")
start <- which(lines=="!series_matrix_table_begin")
end <- which(lines=="!series_matrix_table_end")
if (length(start)!=1L || length(end)!=1L || end <= start+1L)
  stop("Original GEO table boundaries changed")
table_text <- paste(lines[seq.int(start+1L,end-1L)],collapse="\n")
matrix <- utils::read.delim(text=table_text,check.names=FALSE,
                             stringsAsFactors=FALSE,quote="\"",sep="\t")
if (nrow(matrix)!=45101L || ncol(matrix)!=31L ||
    !identical(colnames(matrix)[[1L]],"ID_REF") ||
    !identical(colnames(matrix)[-1L],paste0("GSM",565106:565135)) ||
    anyDuplicated(matrix[[1L]])) stop("Genuine 45101 × 30 source changed")
rownames(matrix) <- matrix[[1L]]
matrix <- as.matrix(matrix[,-1L,drop=FALSE])
storage.mode(matrix) <- "numeric"
stopifnot(all(is.finite(matrix)))
max_delta <- 0
for(i in seq_len(nrow(values))) {
  rec <- values[i,,drop=FALSE]
  probe <- rec$probe_id[[1L]]
  subject <- rec$assay_id[[1L]]
  if (!(probe %in% rownames(matrix)) || !(subject %in% colnames(matrix)))
    stop("Candidate probe or GSM missing from independently read real source")
  column <- matrix[,subject]
  intensity <- unname(matrix[probe,subject])
  percentile <- rank_native(column,intensity)
  difference <- max(abs(intensity-rec$native_processed_log2_RMA_intensity[[1L]]),
                    abs(percentile-rec$within_array_probe_intensity_percentile[[1L]]))
  if (!is.finite(difference) || difference > 1e-12)
    stop(sprintf("Real external native observation/rank changed for %s %s: %.9g",
                 probe,subject,difference))
  max_delta <- max(max_delta,difference)
}
writeLines(sprintf("R independent external native replay PASS; rows=%d; max_absolute_difference=%.17g; projected_model_scores=0",
                   nrow(values), max_delta),
           file.path(audit_dir,"native_external_probe_observations","independent_R_replay.txt"))
cat("PASS: independent R reread all 45101 original probes, 30 independent mice, native quantitative intensities and percentiles agree with Python; NO model score\n")
