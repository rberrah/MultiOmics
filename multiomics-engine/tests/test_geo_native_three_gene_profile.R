#!/usr/bin/env Rscript
# Independent base-R replay of GSE22871 original Affymetrix intensities and
# per-array 3-gene descriptive midranks. This is NOT a transported score,
# external biological replication, or PK/PD evidence.
args <- commandArgs(trailingOnly=TRUE)
midrank <- function(x, v) {
  (sum(x < v) + (sum(x == v) - 1)/2)/(length(x)-1)
}
if (identical(args, "--self-test")) {
  stopifnot(abs(midrank(c(1,2,2,3,4),2)-.375)<1e-15,
            abs(midrank(c(2,4,7),4)-midrank(c(200,400,700),400))<1e-15)
  cat("PASS: base R midrank ties and monotone within-array scaling\n")
  quit(save="no")
}
if (length(args)!=1L) stop("Usage: Rscript test_geo_native_three_gene_profile.R GEO_OUTPUT_DIR")
root <- args[[1L]]
source_path <- file.path(root,"GSE22871_series_matrix.txt.gz")
out <- file.path(root,"external_native_three_gene_descriptive_profile")
csv_raw <- file.path(out,"individual_native_three_gene_probe_measurements.csv")
csv_index <- file.path(out,"individual_native_three_gene_descriptive_index.csv")
audit_path <- file.path(out,"native_three_gene_audit.json")
if (!all(file.exists(c(source_path,csv_raw,csv_index,audit_path))) ||
    !requireNamespace("jsonlite",quietly=TRUE)) stop("Pinned original GEO data/audit absent")
audit <- jsonlite::fromJSON(audit_path)
stopifnot(audit$original_technical_panel_genes==4L,
          audit$unambiguously_mapped_native_probe_genes==3L,
          audit$n_independent_subjects==30L,
          audit$frozen_model_cross_platform_projections==0L,
          audit$pk_covariates_exported==0L,
          !audit$selection_uses_external_dose_genotype_outcomes,
          !audit$native_index_preprocessing_fits_across_subjects,
          !audit$nutrimouse_frozen_four_gene_model_applied)
source_sha <- sub(" .*","",system2("sha256sum",shQuote(source_path),stdout=TRUE))
stopifnot(identical(source_sha,audit$original_gse22871_sha256))
observed <- read.csv(csv_raw,stringsAsFactors=FALSE,check.names=FALSE)
indices <- read.csv(csv_index,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(observed)==90L,nrow(indices)==30L,
          length(unique(observed$subject_id))==30L,
          length(unique(indices$assay_id))==30L,
          identical(unique(observed$assay_id),indices$assay_id),
          identical(unique(observed$original_frozen_feature_label),
                    c("CYP3A11","FAS","CYP4A14")))
# Parse *original* table independently, not Python's serialized calculations.
raw_lines <- readLines(gzfile(source_path),warn=FALSE)
start <- which(raw_lines == "!series_matrix_table_begin")
end <- which(raw_lines == "!series_matrix_table_end")
stopifnot(length(start)==1L,length(end)==1L,start<end,
          end-start-2L==45101L)
original <- read.delim(text=paste(raw_lines[(start+1L):(end-1L)],collapse="\n"),
                       header=TRUE,sep="\t",check.names=FALSE,
                       quote="\"",comment.char="",stringsAsFactors=FALSE)
stopifnot(nrow(original)==45101L,ncol(original)==31L,
          names(original)[[1L]]=="ID_REF",
          identical(names(original)[-1L],indices$assay_id),
          !anyDuplicated(original$ID_REF))
matrix <- as.matrix(original[,-1L,drop=FALSE])
storage.mode(matrix) <- "double"
stopifnot(all(is.finite(matrix)))
gene_names <- c("CYP3A11","FAS","CYP4A14")
for (j in seq_len(30L)) {
  chosen <- observed[observed$assay_id==indices$assay_id[[j]],,drop=FALSE]
  stopifnot(nrow(chosen)==3L,identical(chosen$original_frozen_feature_label,gene_names),
            identical(unique(chosen$subject_id),indices$subject_id[[j]]))
  scores <- numeric(3L)
  for (k in seq_len(3L)) {
    p <- match(chosen$target_probe_id[[k]],original$ID_REF)
    stopifnot(length(p)==1L,!is.na(p))
    value <- matrix[p,j]
    rank <- midrank(matrix[,j],value)
    stopifnot(abs(value-chosen$original_processed_log2_RMA[[k]])<1e-12,
              abs(rank-chosen$within_array_midrank_percentile[[k]])<1e-12)
    scores[[k]] <- rank
  }
  stopifnot(abs(mean(scores)-indices$target_native_three_probe_mean_midrank[[j]])<1e-12,
            indices$classification[[j]]=="descriptive_not_transferred_model_score")
}
writeLines(c("PASS: original 45101 x 30 GEO matrix independent numeric replay",
             "30 genuine biological mice; 90 original target-probe values and midranks",
             "Three-gene index is descriptive, nonfitted and NOT frozen Nutrimouse transfer",
             "0 clinical/PK or externally validated signatures"),
           file.path(out,"independent_R_native_three_gene_replay.txt"))
cat("PASS: 90 native target-probe intensities and percentiles + 30 three-probe indexes reproduced in base R\n")
