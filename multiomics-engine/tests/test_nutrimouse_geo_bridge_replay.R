#!/usr/bin/env Rscript
# Independent implementation of exact source-Nutrimouse/GSE22871 gene identities.
# All probes included; no external group labels, expression data, PK or scores.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)
  stop("Usage: Rscript test_nutrimouse_geo_bridge_replay.R SOURCE_ID_DIR GEO_ENTREZ_DIR PY_BRIDGE_DIR")
src <- args[[1L]]
ref <- args[[2L]]
out <- args[[3L]]
a <- utils::read.csv(file.path(src,"frozen_nutrimouse_genes_source_entrez.csv"),
                     stringsAsFactors=FALSE,check.names=FALSE)
p <- utils::read.csv(file.path(ref,"GPL1261_all_probe_entrez_concordance.csv"),
                     stringsAsFactors=FALSE,check.names=FALSE,na.strings=c("NA"))
py <- utils::read.csv(file.path(out,"locked_panel_to_external_probe_crosswalk.csv"),
                      stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(a)==4L,nrow(p)==45101L,nrow(py)==4L,
          all(!duplicated(p$probe_id)),identical(tolower(a$gene_symbol[[1L]]),"cyp3a11"),
          identical(as.character(a$gene_symbol),as.character(py$source_gene_symbol)))
symbols <- function(value) {
  if(is.na(value) || !nzchar(value))return(character())
  tolower(unlist(strsplit(value,"[,;|[:space:]/]+"),use.names=FALSE))
}
ids <- function(value) {
  value <- as.character(value)
  if(is.na(value) || !nzchar(value))return(character())
  values <- unlist(strsplit(value,"\\|",fixed=FALSE))
  if(any(!grepl("^[1-9][0-9]*$",values)))
    stop("Independent replay refuses nonnumeric source/target Entrez")
  sort(unique(values),method="radix")
}
lines <- c("GENUINE four-feature original Nutrimouse to GEO two-reference R audit")
for(i in seq_len(4L)) {
  source_label <- a$gene_symbol[[i]]
  name <- a$source_canonical_symbol[[i]]
  sid <- ids(a$source_entrez_ids[[i]])
  if(toupper(source_label)=="FAS") {
    stopifnot(identical(name,"Fasn"),identical(sid,"14104"),
              identical(a$original_primary_refseq[[i]],"NM_007988"),
              grepl("DOI_10.1002/hep.21510",
                    a$original_primary_refseq_citation[[i]],fixed=TRUE))
  }
  chosen <- character()
  if(length(sid)==1L && a$source_annotation_status[[i]]=="unique") {
    for(j in seq_len(nrow(p))) {
      geo <- symbols(p$geo_symbols_existing_audit[[j]])
      bioc <- symbols(p$bioc_symbols_existing_audit[[j]])
      if(!(tolower(name) %in% geo && tolower(name) %in% bioc))next
      if(p$entrez_concordance[[j]]!="exact_entrez_set" ||
         tolower(as.character(p$multi_entrez_ambiguous[[j]]))!="false")next
      if(!identical(ids(p$geo_entrez_ids[[j]]),sid) ||
         !identical(ids(p$bioc_entrez_ids[[j]]),sid))next
      chosen <- c(chosen,p$probe_id[[j]])
    }
  }
  chosen <- sort(unique(chosen),method="radix")
  expected <- py$matching_geo_probe_ids[[i]]
  expected <- if(is.na(expected) || !nzchar(expected))character() else
    sort(unlist(strsplit(expected,"\\|")),method="radix")
  if(!identical(chosen,expected) || length(chosen)!=py$strict_probe_count[[i]] ||
     !identical(py$resolved_target_gene_symbol[[i]],name))
    stop(sprintf("Independent four-gene real crosswalk mismatch for %s (%s)",source_label,name))
  lines <- c(lines,sprintf("%s -> %s: %d strictly acceptable target probes [%s]",
                           source_label,name,length(chosen),paste(chosen,collapse="|")))
}
lines <- c(lines,"Zero cross-platform standardised scores; no calibration, PK or group-label optimisation")
writeLines(lines,file.path(out,"independent_R_source_target_gene_bridge_replay.txt"))
cat("PASS: independent R re-evaluated original 45,101 target probes and 4 frozen Nutrimouse genes; full Python identity parity; no PK/external score\n")
cat(paste(lines,collapse="\n"),"\n")
