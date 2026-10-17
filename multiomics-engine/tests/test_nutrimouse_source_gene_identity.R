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
# The REAL frozen macroarray stores source gene symbols in UPPERCASE
# (e.g. CYP3A11), while official mouse SYMBOL uses Cyp3a11. Resolve by
# strictly case-insensitive EXACT official SYMBOL equality, never ALIAS or
# substring, and reject multiple keys instead of choosing the first.
official <- AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db,keytype="SYMBOL")
official <- unique(as.character(official[!is.na(official)]))
keys_by_source <- setNames(lapply(genes,function(symbol) {
  official[tolower(official)==tolower(symbol)]
}),genes)
official_unique <- unique(unlist(keys_by_source[vapply(keys_by_source,length,integer(1))==1L]))
ref <- if(length(official_unique)) {
  suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                 keys=official_unique,keytype="SYMBOL",columns="ENTREZID"))
} else {
  data.frame(SYMBOL=character(),ENTREZID=character())
}
if (!all(c("SYMBOL","ENTREZID") %in% names(ref)))
  stop("Independent mouse reference did not supply symbol/Entrez IDs")
# The source may use an acronym that is ALSO an alias of a different
# official mouse gene. Example: FAS can refer to Fas/CD95 (14102), or
# fatty-acid synthase Fasn (14104). Do not silently promote the exact
# official-SYMBOL match to biological identity without the original
# Nutrimouse macroarray probe sequence or primary probe metadata.
if(!("ALIAS" %in% AnnotationDbi::keytypes(org.Mm.eg.db::org.Mm.eg.db)))
  stop("Cannot audit alternate source names; refuse identity inference")
official_aliases <- AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db,keytype="ALIAS")
official_aliases <- unique(as.character(official_aliases[!is.na(official_aliases)]))
alias_keys_by_source <- setNames(lapply(genes,function(symbol) {
  official_aliases[tolower(official_aliases)==tolower(symbol)]
}),genes)
alias_lookup <- unique(unlist(alias_keys_by_source,use.names=FALSE))
alias_ref <- if(length(alias_lookup)) {
  suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                keys=alias_lookup,keytype="ALIAS",columns=c("SYMBOL","ENTREZID")))
} else data.frame(ALIAS=character(),SYMBOL=character(),ENTREZID=character())
if(!all(c("ALIAS","ENTREZID") %in% names(alias_ref)))
  stop("Cannot trace original source alias collisions")
numeric_ids <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & grepl("^[0-9]+$",values)]
  sort(unique(values),method="radix")
}
exact_ids <- setNames(lapply(genes,function(symbol) {
  exact_keys <- keys_by_source[[symbol]]
  if(length(exact_keys)!=1L)return(character())
  numeric_ids(ref$ENTREZID[!is.na(ref$SYMBOL) & ref$SYMBOL==exact_keys[[1L]]])
}),genes)
alias_ids <- setNames(lapply(genes,function(symbol) {
  matches <- alias_keys_by_source[[symbol]]
  if(!length(matches))return(character())
  numeric_ids(alias_ref$ENTREZID[!is.na(alias_ref$ALIAS) & alias_ref$ALIAS %in% matches])
}),genes)
id_by_gene <- setNames(lapply(genes,function(symbol) {
  sort(unique(c(exact_ids[[symbol]],alias_ids[[symbol]])),method="radix")
}),genes)
# Primary original experimental source (Martin et al., Hepatology 2007,
# doi:10.1002/hep.21510, Table 4, journal page 772 / PDF page 6) explicitly
# annotates the original INRARRAY "FAS" row with mouse RefSeq NM_007988.
# This transcript belongs to fatty-acid synthase Fasn (14104), NOT to the
# cell-death receptor Fas (14102). Cross-verify BOTH with the pinned mouse
# gene database; refuse any disagreement rather than guessing from biology.
primary_refseq <- "NM_007988"
primary_entrez <- "14104"
primary_symbol <- "Fasn"
if(sum(toupper(genes)=="FAS")!=1L ||
   !("REFSEQ" %in% AnnotationDbi::columns(org.Mm.eg.db::org.Mm.eg.db)) ||
   !("ENTREZID" %in% AnnotationDbi::keytypes(org.Mm.eg.db::org.Mm.eg.db)))
  stop("Original FAS primary RefSeq cannot be independently verified")
primary_annotation <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db, keys=c("14102","14104"),
  keytype="ENTREZID", columns=c("REFSEQ","SYMBOL")))
transcript <- as.character(primary_annotation$REFSEQ)
transcript <- sub("\\.[0-9]+$", "", transcript)
verified_entrez <- unique(as.character(primary_annotation$ENTREZID[
  !is.na(transcript) & transcript==primary_refseq]))
verified_symbol <- unique(as.character(primary_annotation$SYMBOL[
  !is.na(transcript) & transcript==primary_refseq]))
if(!identical(verified_entrez, primary_entrez) ||
   !identical(verified_symbol,primary_symbol) ||
   !(primary_entrez %in% alias_ids[["FAS"]]) ||
   !(primary_entrez %in% id_by_gene[["FAS"]]) ||
   !( "14102" %in% exact_ids[["FAS"]]) )
  stop("Primary original FAS NM_007988 is not uniquely assigned to Fasn/14104")
id_by_gene[["FAS"]] <- primary_entrez
canonical_by_gene <- setNames(vapply(genes,function(symbol) {
  if(toupper(symbol)=="FAS")return(primary_symbol)
  matches <- keys_by_source[[symbol]]
  if(length(matches)==1L)matches[[1L]] else symbol
},character(1)),genes)
source_refseq <- setNames(ifelse(toupper(genes)=="FAS", primary_refseq, ""),genes)

statuses <- vapply(genes,function(symbol) {
  matches <- keys_by_source[[symbol]]
  official_ids <- exact_ids[[symbol]]
  alternatives <- alias_ids[[symbol]]
  # RefSeq in the ORIGINAL publication resolves the "FAS" label.
  # This makes gene ID technically traceable, but NOT its exact nylon probe
  # sequence/specificity, nor the Affymetrix measurement-scale equivalence.
  if(toupper(symbol)=="FAS")return("unique")
  if(length(matches)==1L && length(official_ids)==1L &&
     !length(setdiff(alternatives,official_ids)))return("unique")
  if(length(matches)==1L && length(setdiff(alternatives,official_ids))>0L)
    return("official_symbol_vs_alias_conflict")
  if(length(matches)>1L)return("multiple_casefold_official_symbols")
  if(length(matches)==0L && length(alternatives)>0L)
    return("alias_only_requires_primary_probe_identity")
  if(length(id_by_gene[[symbol]])>1L)return("multi_gene_ambiguous")
  "unmapped"
},character(1))
rows <- data.frame(
  model_feature_index=seq_along(genes),
  gene_symbol=genes,
  source_canonical_symbol=unname(canonical_by_gene),
  original_primary_refseq=unname(source_refseq),
  original_primary_refseq_citation=ifelse(toupper(genes)=="FAS",
    "Martin_2007_Hepatology_45_767_777_Table4_p772_DOI_10.1002/hep.21510",
    "not_specified_in_primary_table_for_this_frozen_feature"),
  source_entrez_ids=vapply(id_by_gene,paste,collapse="|",FUN.VALUE=character(1)),
  source_entrez_count=vapply(id_by_gene,length,FUN.VALUE=integer(1)),
  source_annotation_status=unname(statuses),
  exact_official_entrez=vapply(exact_ids,paste,collapse="|",FUN.VALUE=character(1)),
  alias_entrez_candidates=vapply(alias_ids,paste,collapse="|",FUN.VALUE=character(1)),
  source_annotation="org.Mm.eg.db SYMBOL plus ALIAS collision screen",
  source_annotation_version=as.character(utils::packageVersion("org.Mm.eg.db")),
  stringsAsFactors=FALSE)
if(!identical(rows$source_entrez_ids[rows$gene_symbol=="FAS"],"14104") ||
   !identical(rows$source_canonical_symbol[rows$gene_symbol=="FAS"],"Fasn") ||
   !identical(rows$original_primary_refseq[rows$gene_symbol=="FAS"],"NM_007988"))
  stop("Original published probe annotation still ambiguous: refuse false Fas 14102 mapping")
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
  any_incomplete_source_mapping=any(rows$source_annotation_status!="unique"),
  unresolved_historical_FAS_acronym=FALSE,
  original_FAS_primary_refseq="NM_007988",
  original_FAS_source_gene_canonical="Fasn",
  original_FAS_gene_entrez="14104",
  primary_FAS_refseq_verified_against_orgMmEgDb=TRUE,
  primary_source_citation="https://doi.org/10.1002/hep.21510#Table-4-page-772",
  status="SOURCE_GENE_IDS_ONLY_NO_CROSS_PLATFORM_EXPRESSION_CALIBRATION",
  external_signature_scores=0L)
jsonlite::write_json(aud,file.path(out,"frozen_source_gene_identity_audit.json"),
                     auto_unbox=TRUE,pretty=TRUE)
cat("REAL frozen Nutrimouse 30/10 model genes + independent org.Mm.eg.db: ",
    paste(paste(genes,rows$source_entrez_ids,sep=":"),collapse=", "),
    "; mapped uniquely=",sum(rows$source_annotation_status=="unique"),
    "; no external signatures or PK claims\n",sep="")
