#!/usr/bin/env Rscript
# Map the REAL frozen Nutrimouse 4-gene panel to an independently versioned
# mouse Entrez gene reference, without consulting GSE22871 phenotypes.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L)
  stop("Usage: Rscript test_nutrimouse_source_gene_identity.R FROZEN_NUTRIMOUSE_DIR OUT_DIR")
root <- args[[1L]]
out <- args[[2L]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if (!requireNamespace("jsonlite",quietly=TRUE) ||
    !requireNamespace("AnnotationDbi",quietly=TRUE) ||
    !requireNamespace("org.Mm.eg.db",quietly=TRUE))
  stop("The exact R and mouse gene reference packages are required")
model_file <- file.path(root,"frozen_panel_model.json")
definition_file <- file.path(root,"panel_feature_definition.json")
report_file <- file.path(root,"frozen_score_audit.json")
if (!file.exists(model_file) || !file.exists(definition_file))
  stop("Training-frozen public Nutrimouse model evidence absent")
model <- jsonlite::fromJSON(model_file)
definition <- jsonlite::fromJSON(definition_file)
genes <- as.character(unlist(model$feature_ids))
if (length(genes)!=4L || anyNA(genes) || anyDuplicated(tolower(genes)) ||
    tolower(genes[[1L]])!="cyp3a11" ||
    !identical(genes,as.character(unlist(definition$selected_gene_ids))) ||
    !identical(model$pinned_source_commit,
              "ef3e760526623d9e91e945e4b50be48d764efc40") ||
    length(model$training_subject_ids)!=30L ||
    length(model$heldout_subject_ids)!=10L ||
    length(intersect(model$training_subject_ids,model$heldout_subject_ids))!=0L ||
    !identical(model$phenotype_labels_used_for_selection_or_fitting,FALSE) ||
    !identical(model$heldout_values_used_for_selection_or_fitting,FALSE))
  stop("Model gene feature/split identity is not the original frozen 30/10 training-only artifact")
if (!all(c("SYMBOL","ENTREZID") %in% AnnotationDbi::columns(org.Mm.eg.db::org.Mm.eg.db)))
  stop("Mouse source annotation reference lacks SYMBOL/ENTREZID")
ref <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
               keys=genes,keytype="SYMBOL",columns="ENTREZID"))
if (!all(c("SYMBOL","ENTREZID") %in% names(ref)))
  stop("Independent mouse reference did not supply symbol/Entrez IDs")
id_by_gene <- setNames(lapply(genes,function(symbol) {
  ids <- as.character(ref$ENTREZID[!is.na(ref$SYMBOL) & ref$SYMBOL==symbol])
  ids <- ids[!is.na(ids) & grepl("^[0-9]+$",ids)]
  sort(unique(ids),method="radix")
}), genes)
rows <- data.frame(
  model_feature_index=seq_along(genes),
  gene_symbol=genes,
  source_entrez_ids=vapply(id_by_gene,paste,collapse="|",FUN.VALUE=character(1)),
  source_entrez_count=vapply(id_by_gene,length,FUN.VALUE=integer(1)),
  source_annotation_status=vapply(id_by_gene,function(ids) {
    if(length(ids)==0L) "unmapped"
    else if(length(ids)==1L) "unique"
    else "multi_gene_ambiguous"
  },character(1)),
  source_annotation="org.Mm.eg.db SYMBOL->ENTREZID",
  source_annotation_version=as.character(utils::packageVersion("org.Mm.eg.db")),
  stringsAsFactors=FALSE)
utils::write.csv(rows,file.path(out,"frozen_nutrimouse_genes_source_entrez.csv"),
                 row.names=FALSE,na="NA")
source_sha <- function(f) {
  x <- suppressWarnings(system2("sha256sum",shQuote(normalizePath(f)),
                         stdout=TRUE,stderr=TRUE))
  if(length(x)!=1L || !grepl("^[a-f0-9]{64}  ",x))
    stop("Cannot independently audit SHA-256 of original source/model")
  substr(x,1L,64L)
}
aud <- list(source_dataset="real_public_mixOmics_nutrimouse",
  pinned_source_commit=model$pinned_source_commit,
  source_rda_sha256=model$source_rda_sha256,
  model_sha256=source_sha(model_file),
  panel_definition_sha256=source_sha(definition_file),
  training_n=30L,heldout_n=10L,
  selected_feature_names=genes,
  features_selected_training_only=TRUE,
  biological_feature_selection=FALSE,
  source_annotation_version=as.character(utils::packageVersion("org.Mm.eg.db")),
  source_mapping_is_pheno_blind=TRUE,
  any_incomplete_source_mapping=any(rows$source_entrez_count!=1L),
  status="SOURCE_GENE_IDS_ONLY_NO_CROSS_PLATFORM_EXPRESSION_CALIBRATION",
  external_signature_scores=0L)
jsonlite::write_json(aud,file.path(out,"frozen_source_gene_identity_audit.json"),
                     auto_unbox=TRUE,pretty=TRUE)
cat("REAL frozen Nutrimouse 30/10 model genes + independent org.Mm.eg.db: ",
    paste(paste(genes,rows$source_entrez_ids,sep=":"),collapse=", "),
    "; mapped uniquely=",sum(rows$source_entrez_count==1L),
    "; no external signatures or PK claims\n",sep="")
