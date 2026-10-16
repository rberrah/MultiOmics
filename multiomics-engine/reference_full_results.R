# Complete R reference-method result tables, for use from the consented LOCAL R
# backend only. No external web services or global caches are involved.
# Scientific contract: show a top-100 preview, but do NOT disguise it as
# an exhaustive result. Every fitted feature and every fitted result column
# must be preserved, including non-significant genes and missing q values.
# Optional, operator-controlled retention for an oversized, sensitive result.
# No HTTP download endpoint is exposed. The researcher must explicitly set
# PMX_MULTIOMICS_RESULTS_DIR on their own machine and retrieve the CSV locally.
# Random subdirectories have mode 0700; filenames never contain user-supplied IDs.
archive_large_reference_table <- function(csv_path, archive_dir, max_local_bytes=1024^3) {
  if(!is.character(archive_dir)||length(archive_dir)!=1L||
     is.na(archive_dir)||!nzchar(trimws(archive_dir)))
    return(list(saved=FALSE,reason="not_configured"))
  if(!dir.exists(archive_dir)||file.access(archive_dir,2L)!=0L)
    return(list(saved=FALSE,reason="directory_missing_or_unwritable"))
  bytes <- unname(file.info(csv_path)$size)
  if(!is.finite(bytes)||bytes>max_local_bytes)
    return(list(saved=FALSE,reason="local_archive_size_limit"))
  if(!requireNamespace("digest",quietly=TRUE))
    return(list(saved=FALSE,reason="sha256_digest_package_unavailable"))
  private <- tempfile("pmx_reference_",tmpdir=normalizePath(archive_dir,mustWork=TRUE))
  if(!dir.create(private,mode="0700",showWarnings=FALSE))
    return(list(saved=FALSE,reason="could_not_create_private_directory"))
  saved <- file.path(private,"fitted_reference.csv")
  if(!file.copy(csv_path,saved,overwrite=FALSE,copy.mode=FALSE)){
    unlink(private,recursive=TRUE,force=TRUE)
    return(list(saved=FALSE,reason="file_copy_failed"))
  }
  Sys.chmod(private,mode="0700")
  Sys.chmod(saved,mode="0600")
  sha <- digest::digest(file=saved,algo="sha256",serialize=FALSE)
  list(saved=TRUE,id=basename(private),bytes=bytes,sha256=sha)
}

reference_complete_table <- function(x, original_features=NA_integer_,
                                     max_wire_bytes=64L*1024L*1024L,
                                     local_archive_dir=Sys.getenv("PMX_MULTIOMICS_RESULTS_DIR",""),
                                     max_local_bytes=1024^3) {
  if(is.null(x))stop("Complete reference export requires fitted results.")
  tab <- as.data.frame(x, check.names=FALSE, stringsAsFactors=FALSE)
  if(!"feature" %in% names(tab))stop("Reference result lacks feature IDs.")
  ids <- as.character(tab$feature)
  if(anyNA(ids)||any(!nzchar(ids))||anyDuplicated(ids))
    stop("Complete reference result contains missing or duplicated feature IDs.")
  if(nrow(tab)!=length(ids))stop("Inconsistent reference feature count.")
  if(!is.numeric(max_wire_bytes) || length(max_wire_bytes)!=1L ||
     !is.finite(max_wire_bytes) || max_wire_bytes < 256L)
    stop("Invalid maximum wire export size.")
  value_types <- vapply(tab,function(y) {
    if(is.integer(y))"integer" else if(is.numeric(y))"numeric"
    else if(is.logical(y))"logical" else "string"
  },character(1))
  all_rows <- nrow(tab)
  # This estimate is only a preflight; the final exact byte count is
  # verified below. This guard prevents accidentally serializing a multi-GB
  # study into one HTTP JSON response.
  if(as.double(object.size(tab))>max_wire_bytes*3 && !nzchar(local_archive_dir)){
    return(list(referenceTableStatus="too_large",referenceTableComplete=FALSE,
      referenceTableRows=all_rows,referenceTableExportedRows=0L,
      referenceTableColumns=as.list(names(tab)),
      referenceTableTypes=as.list(value_types),
      referenceTableCsv=NULL,referenceTableMessage=
      "R reference result exceeds the HTTP budget and was NOT retained by this API. To preserve it, explicitly configure PMX_MULTIOMICS_RESULTS_DIR on the local R server and rerun. Top rows are only a preview."))
  }
  # Stream to a temporary LOCAL CSV first, so a large payload can be rejected
  # before an extra full-size string is allocated in the HTTP JSON response.
  export_path <- tempfile("pmx_reference_full_",fileext=".csv")
  on.exit(unlink(export_path),add=TRUE)
  utils::write.csv(tab,file=export_path,row.names=FALSE,na="NA",quote=TRUE)
  bytes <- unname(file.info(export_path)$size)
  if(!is.finite(bytes))stop("Unable to determine full R CSV size.")
  if(bytes>max_wire_bytes){
    local <- archive_large_reference_table(export_path,local_archive_dir,max_local_bytes)
    locally_saved <- isTRUE(local$saved)
    return(list(
      referenceTableStatus=if(locally_saved)"too_large_archived_locally" else "too_large",
      referenceTableComplete=FALSE,
      referenceTableRows=all_rows,referenceTableExportedRows=0L,
      referenceTableColumns=as.list(names(tab)),
      referenceTableTypes=as.list(value_types),
      referenceTableCsv=NULL,referenceTableBytes=bytes,
      referenceTableLocalArchiveId=if(locally_saved)local$id else NULL,
      referenceTableLocalArchiveSha256=if(locally_saved)local$sha256 else NULL,
      referenceTableLocalArchiveBytes=if(locally_saved)local$bytes else NULL,
      referenceTableMessage=if(locally_saved)
        paste("Too large for JSON: full CSV saved ONLY in the explicitly configured local",
              "PMX_MULTIOMICS_RESULTS_DIR directory under a private subfolder;",
              "not included in browser ZIP. Verify its SHA-256 and retrieve locally.")
      else paste("Too large for JSON, NOT retained by this API (reason:",
                 local$reason,"). Configure PMX_MULTIOMICS_RESULTS_DIR in a",
                 "private local R environment and rerun to preserve the full CSV.")
    ))
  }
  csv <- readChar(export_path,nchars=as.integer(bytes),useBytes=TRUE)
  if(!identical(nchar(csv,type="bytes"),as.integer(bytes)))
    stop("Full R CSV readback mismatch.")
  if(!is.na(original_features)){
    original_features<-as.integer(original_features)
    if(original_features<all_rows)
      stop("Fitted reference result unexpectedly has more feature rows than input.")
  }
  list(referenceTableStatus="complete",referenceTableComplete=TRUE,
    referenceTableRows=all_rows,referenceTableExportedRows=all_rows,
    referenceTableOriginalFeatures=if(is.na(original_features))NULL else original_features,
    referenceTableExcludedBeforeFit=if(is.na(original_features))NULL else
      original_features-all_rows,
    referenceTableColumns=as.list(names(tab)),referenceTableTypes=as.list(value_types),
    referenceTableCsv=csv,referenceTableBytes=bytes,
    referenceTableCsvFormat="RFC4180-quoted UTF-8, NA token 'NA', unmodified original feature IDs",
    referenceTableMessage="This is the complete fitted reference table, not a biological or statistical certification. Raw IDs are unchanged; do not open untrusted CSV in spreadsheet software without formula-safety controls."
  )
}
