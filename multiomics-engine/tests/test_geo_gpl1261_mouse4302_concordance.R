#!/usr/bin/env Rscript
# Genuine independent Bioconductor mouse4302.db concordance on the ORIGINAL
# GEO GPL1261 probe IDs for real GSE22871. This is an annotation audit ONLY:
# changing gene symbols, one-to-many mappings and platform scales are
# explicitly NOT resolved to invent validated external scores.
args <- commandArgs(trailingOnly=TRUE)

clean_set <- function(x) {
  x <- as.character(x)
  sort(unique(x[!is.na(x) & nzchar(trimws(x)) & !(toupper(x) %in% c("NA","---","NULL"))]))
}
signature <- function(x) paste(clean_set(x),collapse="|")
status_for <- function(a,b) {
  a <- clean_set(a); b <- clean_set(b)
  if (!length(a) && !length(b)) return("both_unmapped")
  if (!length(a)) return("bioc_only")
  if (!length(b)) return("geo_only")
  if (identical(a,b)) return("exact_symbol_set_agreement")
  if (length(intersect(a,b))) return("partial_symbol_overlap")
  "disjoint_symbol_sets"
}

self_test <- function() {
  stopifnot(
    status_for("Cyp3a11","Cyp3a11")=="exact_symbol_set_agreement",
    status_for(c("A","B"),c("B","C"))=="partial_symbol_overlap",
    status_for("A","C")=="disjoint_symbol_sets",
    status_for("A",NA_character_)=="geo_only",
    status_for(NA_character_,"B")=="bioc_only",
    status_for(NA_character_,NA_character_)=="both_unmapped",
    signature(c("B","A","A",NA_character_))=="A|B"
  )
  cat("PASS: independent Bioconductor comparison preserves both-side, partial, disjoint, missing and ambiguous probe mappings\n")
}
if (length(args)==1L && identical(args[[1L]],"--self-test")) {
  self_test()
  quit(save="no",status=0L)
}
if (length(args)!=2L) stop("Usage: Rscript .../test_geo_gpl1261_mouse4302_concordance.R GEO_OUTPUT_DIR AUDIT_OUTPUT_DIR")
in_dir <- args[[1L]]; out_dir <- args[[2L]]
if (!requireNamespace("AnnotationDbi",quietly=TRUE) ||
    !requireNamespace("mouse4302.db",quietly=TRUE) ||
    !requireNamespace("jsonlite",quietly=TRUE) ||
    !requireNamespace("BiocManager",quietly=TRUE))
  stop("True Bioconductor mouse4302.db, AnnotationDbi, jsonlite and BiocManager are mandatory; no surrogate.")
if (!identical(as.character(getRversion()),"4.5.2") ||
    !identical(as.character(BiocManager::version()),"3.22"))
  stop("R 4.5.2 and Bioconductor 3.22 must remain pinned")
if (!identical(as.character(utils::packageVersion("mouse4302.db")),"3.13.0"))
  stop("The Bioconductor annotation package version is not frozen at mouse4302.db 3.13.0")
db <- get("mouse4302.db",envir=asNamespace("mouse4302.db"))
if (!all(c("PROBEID","SYMBOL","ENTREZID") %in% AnnotationDbi::keytypes(db)) ||
    !all(c("SYMBOL","ENTREZID") %in% AnnotationDbi::columns(db)))
  stop("Published Bioconductor chip annotation lacks expected PROBEID/SYMBOL/ENTREZID contract")
source_file <- file.path(in_dir,"GPL1261_probe_to_gene_all.csv")
metadata_file <- file.path(in_dir,"GPL1261_annotation_mapping_audit.json")
cohort_file <- file.path(in_dir,"external_geo_intake_audit.json")
if (!all(file.exists(c(source_file,metadata_file,cohort_file))))
  stop("Original source-derived GEO mapping and real cohort audit are absent")
geo_audit <- jsonlite::read_json(metadata_file,simplifyVector=TRUE)
cohort <- jsonlite::read_json(cohort_file,simplifyVector=TRUE)
stopifnot(identical(geo_audit$official_annotation_gzip_sha256,
                    "394c103cf4ca99500d06369e4fb2b9b9bffc2c10fd3cbeb0b0d64c51c3fc95f8"),
          identical(geo_audit$external_cohort_source_sha256,
                    "630706c5d3778068a9e08af01ce315b4886de17abefbfaae32708e7ed3d2d861"),
          identical(cohort$n_biological_subjects,30L) ||
            isTRUE(all.equal(cohort$n_biological_subjects,30)))
geo <- utils::read.csv(source_file,stringsAsFactors=FALSE,check.names=FALSE,
                       colClasses="character",na.strings=c("NA"))
expected <- c("platform_probe_id","gene_symbol","probe_mapping_status")
if (!identical(colnames(geo),expected) ||
    anyNA(geo$platform_probe_id) || any(!nzchar(geo$platform_probe_id)))
  stop("GEO original platform crosswalk schema changed")
probes <- sort(unique(geo$platform_probe_id))
if (length(probes)!=45101L || anyDuplicated(geo[,c("platform_probe_id","gene_symbol")]))
  stop("Expected precisely 45,101 source probes and no duplicated symbol/probe pairs")
# Genuine chip database lookup. Keep its independent NCBI Entrez IDs, but do
# not silently overwrite the original GEO mapping with newer annotation.
chunk <- split(probes,ceiling(seq_along(probes)/700L))
results <- lapply(chunk,function(ids) {
  ans <- AnnotationDbi::select(db,keys=ids,columns=c("SYMBOL","ENTREZID"),
                               keytype="PROBEID")
  if (!all(c("PROBEID","SYMBOL","ENTREZID") %in% colnames(ans)))
    stop("Bioconductor returned a nonconforming annotation table")
  ans[,c("PROBEID","SYMBOL","ENTREZID"),drop=FALSE]
})
bioc <- do.call(rbind,results)
if (!all(probes %in% bioc$PROBEID))
  stop("Bioconductor returned no row for some measured platform probes")
bioc$PROBEID <- as.character(bioc$PROBEID)
bioc$SYMBOL <- as.character(bioc$SYMBOL)
bioc$ENTREZID <- as.character(bioc$ENTREZID)
geo_by_probe <- split(geo$gene_symbol,geo$platform_probe_id)
bio_symbols <- split(bioc$SYMBOL,bioc$PROBEID)
bio_entrez <- split(bioc$ENTREZID,bioc$PROBEID)
rows <- lapply(probes,function(probe) {
  geo_names <- clean_set(geo_by_probe[[probe]])
  bio_names <- clean_set(bio_symbols[[probe]])
  ids <- clean_set(bio_entrez[[probe]])
  data.frame(
    probe_id=probe,
    geo_symbols=signature(geo_names),
    bioc_symbols=signature(bio_names),
    bioc_entrez_ids=signature(ids),
    geo_symbol_count=length(geo_names),
    bioc_symbol_count=length(bio_names),
    bioc_entrez_count=length(ids),
    concordance_status=status_for(geo_names,bio_names),
    stringsAsFactors=FALSE
  )
})
audit <- do.call(rbind,rows)
stopifnot(nrow(audit)==45101L,!anyDuplicated(audit$probe_id))
states <- c("exact_symbol_set_agreement","partial_symbol_overlap",
            "disjoint_symbol_sets","geo_only","bioc_only","both_unmapped")
if (!all(audit$concordance_status %in% states))
  stop("Unknown independent annotation disagreement class")
counts <- table(factor(audit$concordance_status,levels=states))
# The correctness gate requires real independent comparison, not forced
# biological agreement: annotation revisions might legitimately disagree.
stopifnot(sum(counts)==45101L,
          any(audit$bioc_symbol_count>0L),
          any(audit$geo_symbol_count>0L),
          any(audit$bioc_entrez_count>0L))
anchor <- audit[tolower(audit$geo_symbols)=="cyp3a11" |
                  grepl("(^|\\|)Cyp3a11(\\||$)",audit$bioc_symbols,ignore.case=TRUE),
                ,drop=FALSE]
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(audit,file.path(out_dir,"GPL1261_GEO_vs_mouse4302_db_all_probe_concordance.csv"),
                 row.names=FALSE,na="NA")
utils::write.csv(audit[!audit$concordance_status %in% "exact_symbol_set_agreement",],
                 file.path(out_dir,"GPL1261_all_annotation_discrepancies.csv"),
                 row.names=FALSE,na="NA")
utils::write.csv(anchor,file.path(out_dir,"Cyp3a11_candidate_probe_concordance.csv"),
                 row.names=FALSE,na="NA")
pkg <- c("mouse4302.db","AnnotationDbi","org.Mm.eg.db")
installed <- pkg[vapply(pkg,requireNamespace,logical(1),quietly=TRUE)]
versions <- as.list(vapply(installed,function(p)
                           as.character(utils::packageVersion(p)),character(1)))
report <- list(
  study="GSE22871",platform="GPL1261",
  source_geo_annotation_sha256=geo_audit$official_annotation_gzip_sha256,
  source_external_cohort_sha256=geo_audit$external_cohort_source_sha256,
  reference="Bioconductor mouse4302.db PROBEID/SYMBOL/ENTREZID",
  pinned_r_version=as.character(getRversion()),
  pinned_bioconductor_version=as.character(BiocManager::version()),
  actual_package_versions=versions,
  probes_compared=nrow(audit),
  statuses=as.list(as.integer(counts)),
  status_names=names(counts),
  cyp3a11_candidate_probes=nrow(anchor),
  cyp3a11_candidate_rows=anchor,
  comparison_scope="GEO symbol-string equality vs independent updated Bioconductor SYMBOL; probe-wise without guessing aliases or gene ID identity",
  gene_identifier_caveat="Entrez IDs are stored for reference from Bioconductor, but they were NOT compared against a pinned GEO Gene ID crosswalk.",
  ambiguity_policy="Preserve all original and independent mappings; no label-based representative probe or consensus imputation",
  bioconductor_vs_geo_discrepancy_expected=TRUE,
  cross_platform_projection="BLOCKED_UNTIL_ANNOTATION_CONSENSUS_AND_SCALE_TRANSFER",
  score_rows_generated=0L,
  biological_or_pk_validation=FALSE
)
jsonlite::write_json(report,file.path(out_dir,"GPL1261_mouse4302_db_concordance_audit.json"),
                     pretty=TRUE,auto_unbox=TRUE,na="null",digits=12L)
cat(sprintf(
  "REAL mouse4302.db vs GEO PROBE concordance PASS: 45101 compared; exact=%d partial=%d disjoint=%d GEOonly=%d BiocOnly=%d bothUnmapped=%d; CYP3A11 candidates=%d; 0 projected scores\n",
  counts[[1]],counts[[2]],counts[[3]],counts[[4]],counts[[5]],counts[[6]],nrow(anchor)
))
