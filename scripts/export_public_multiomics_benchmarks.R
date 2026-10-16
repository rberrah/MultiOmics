args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "tmp/public-benchmarks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mixomics_root <- Sys.getenv("MIXOMICS_ROOT", "external/mixOmics")

write_matrix <- function(x, feature_ids, assay_ids, path) {
  stopifnot(nrow(x) == length(assay_ids))
  stopifnot(ncol(x) == length(feature_ids))
  out <- data.frame(feature_id = feature_ids, t(as.matrix(x)), check.names = FALSE)
  colnames(out) <- c("feature_id", assay_ids)
  write.csv(out, path, row.names = FALSE, quote = FALSE)
}

write_metadata <- function(subject_ids, assay_ids_by_omic, condition, out_path, extra = list()) {
  rows <- list()
  idx <- 1L
  for (i in seq_along(subject_ids)) {
    for (omic in names(assay_ids_by_omic)) {
      row <- data.frame(
        subject_id = subject_ids[[i]],
        sample_id = subject_ids[[i]],
        assay_id = assay_ids_by_omic[[omic]][[i]],
        omic = omic,
        condition = as.character(condition[[i]]),
        timepoint = "",
        batch = "",
        technical_replicate = "1",
        outcome = if ("outcome" %in% names(extra)) as.character(extra$outcome[[i]]) else "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      if (length(extra)) {
        for (nm in setdiff(names(extra), "outcome")) row[[nm]] <- as.character(extra[[nm]][[i]])
      }
      rows[[idx]] <- row
      idx <- idx + 1L
    }
  }
  write.csv(do.call(rbind, rows), out_path, row.names = FALSE, quote = FALSE)
}

# --- Nutrimouse ---
load(file.path(mixomics_root, "data", "nutrimouse.rda"))
nutri_dir <- file.path(out_dir, "nutrimouse")
dir.create(nutri_dir, recursive = TRUE, showWarnings = FALSE)

n <- nrow(nutrimouse$gene)
subjects <- sprintf("MOUSE%03d", seq_len(n))
rna_ids <- sprintf("RNA%03d", seq_len(n))
met_ids <- sprintf("MET%03d", seq_len(n))

write_matrix(nutrimouse$gene, colnames(nutrimouse$gene), rna_ids, file.path(nutri_dir, "transcriptomics.csv"))
write_matrix(nutrimouse$lipid, colnames(nutrimouse$lipid), met_ids, file.path(nutri_dir, "metabolomics.csv"))
write_metadata(
  subjects,
  list(transcriptomics = rna_ids, metabolomics = met_ids),
  nutrimouse$genotype,
  file.path(nutri_dir, "metadata.csv"),
  extra = list(diet = nutrimouse$diet)
)

cat("NUTRIMOUSE\n")
cat("subjects:", n, "\n")
cat("genotype:", paste(names(table(nutrimouse$genotype)), as.integer(table(nutrimouse$genotype)), collapse = "; "), "\n")
cat("diet:", paste(names(table(nutrimouse$diet)), as.integer(table(nutrimouse$diet)), collapse = "; "), "\n")
cat("gene_dim:", paste(dim(nutrimouse$gene), collapse = "x"), "\n")
cat("lipid_dim:", paste(dim(nutrimouse$lipid), collapse = "x"), "\n")
cat("gene_range:", paste(signif(range(as.matrix(nutrimouse$gene), na.rm = TRUE), 5), collapse = ".."), "\n")
cat("lipid_range:", paste(signif(range(as.matrix(nutrimouse$lipid), na.rm = TRUE), 5), collapse = ".."), "\n")
cat("gene_names:", paste(colnames(nutrimouse$gene), collapse = ","), "\n")
cat("lipid_names:", paste(colnames(nutrimouse$lipid), collapse = ","), "\n")

# --- breast.TCGA: HER2-enriched vs LumA ---
load(file.path(mixomics_root, "data", "breast.TCGA.rda"))
tcga <- breast.TCGA$data.train
keep <- as.character(tcga$subtype) %in% c("Her2", "LumA")
mrna <- tcga$mrna[keep, , drop = FALSE]
protein <- tcga$protein[keep, , drop = FALSE]
subtype <- droplevels(tcga$subtype[keep])

tcga_dir <- file.path(out_dir, "tcga-her2-luma")
dir.create(tcga_dir, recursive = TRUE, showWarnings = FALSE)

n2 <- nrow(mrna)
subjects2 <- sprintf("TCGA%03d", seq_len(n2))
rna2 <- sprintf("RNA%03d", seq_len(n2))
prot2 <- sprintf("PROT%03d", seq_len(n2))

write_matrix(mrna, colnames(mrna), rna2, file.path(tcga_dir, "transcriptomics.csv"))
write_matrix(protein, colnames(protein), prot2, file.path(tcga_dir, "proteomics.csv"))
write_metadata(
  subjects2,
  list(transcriptomics = rna2, proteomics = prot2),
  subtype,
  file.path(tcga_dir, "metadata.csv"),
  extra = list(outcome = subtype)
)

cat("\nTCGA_HER2_LUMA\n")
cat("subjects:", n2, "\n")
cat("subtype:", paste(names(table(subtype)), as.integer(table(subtype)), collapse = "; "), "\n")
cat("mrna_dim:", paste(dim(mrna), collapse = "x"), "\n")
cat("protein_dim:", paste(dim(protein), collapse = "x"), "\n")
cat("mrna_range:", paste(signif(range(as.matrix(mrna), na.rm = TRUE), 5), collapse = ".."), "\n")
cat("protein_range:", paste(signif(range(as.matrix(protein), na.rm = TRUE), 5), collapse = ".."), "\n")

marker_pattern <- "(?i)(ERBB2|HER2|GRB7|STARD3|PGAP3)"
cat("mrna_marker_features:", paste(grep(marker_pattern, colnames(mrna), value = TRUE, perl = TRUE), collapse = ","), "\n")
cat("protein_marker_features:", paste(grep(marker_pattern, colnames(protein), value = TRUE, perl = TRUE), collapse = ","), "\n")


# --- Published cross-omics truth pairs: NCI-60 IntLIM ---
nci_gm_root <- Sys.getenv("NCI60_GM_ROOT", "external/NCI60_GeneMetabolite_Data")
if (dir.exists(nci_gm_root)) {
  gene <- read.csv(file.path(nci_gm_root, "geneData.csv"), row.names = 1, check.names = FALSE)
  metab <- read.csv(file.path(nci_gm_root, "metabData.csv"), row.names = 1, check.names = FALSE)
  pdata <- read.csv(file.path(nci_gm_root, "pData.csv"), check.names = FALSE, stringsAsFactors = FALSE)

  export_pair <- function(gene_name, metab_name, filename) {
    stopifnot(gene_name %in% rownames(gene), metab_name %in% rownames(metab))
    p <- pdata[pdata$cancergroup %in% c("BPO", "Leukemia"), c("cell_line", "cancergroup")]
    common <- Reduce(intersect, list(p$cell_line, colnames(gene), colnames(metab)))
    p <- p[match(common, p$cell_line), , drop = FALSE]
    out <- data.frame(
      sample_id = common,
      condition = p$cancergroup,
      gene = as.numeric(gene[gene_name, common]),
      metabolite = as.numeric(metab[metab_name, common]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    write.csv(out, file.path(out_dir, filename), row.names = FALSE, quote = FALSE)
  }

  export_pair("FAM174B", "malic acid", "nci60_FAM174B_malic.csv")
  export_pair("DNER", "L-beta-imidazolelactic acid", "nci60_DNER_imidazole.csv")
  cat("\nNCI60_INTLIM_PAIRS\n")
  cat("gene_dim:", paste(dim(gene), collapse = "x"), "\n")
  cat("metab_dim:", paste(dim(metab), collapse = "x"), "\n")
}

# --- Published cross-omics truth pair: breast tumour vs normal ---
intlim_root <- Sys.getenv("INTLIM_ROOT", "external/IntLIMVignettes")
brca_root <- file.path(intlim_root, "BRCA_data")
if (dir.exists(brca_root)) {
  gene <- read.csv(file.path(brca_root, "geneData.csv"), row.names = 1, check.names = FALSE)
  metab <- read.csv(file.path(brca_root, "metabData.csv"), row.names = 1, check.names = FALSE)
  pdata <- read.csv(file.path(brca_root, "pData.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("GPT2" %in% rownames(gene), "2-hydroxyglutarate" %in% rownames(metab))

  p <- pdata[pdata$DIAG %in% c("NORMAL", "TUMOR"), c("id", "DIAG")]
  common <- Reduce(intersect, list(p$id, colnames(gene), colnames(metab)))
  p <- p[match(common, p$id), , drop = FALSE]
  out <- data.frame(
    sample_id = common,
    condition = p$DIAG,
    gene = as.numeric(gene["GPT2", common]),
    metabolite = as.numeric(metab["2-hydroxyglutarate", common]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  write.csv(out, file.path(out_dir, "brca_GPT2_2HG.csv"), row.names = FALSE, quote = FALSE)
  cat("\nBRCA_INTLIM_PAIR\n")
  cat("paired_samples:", nrow(out), "\n")
}


# --- NCI-60 incomplete multi-omics from Bioconductor missRows ---
missrows_root <- Sys.getenv("MISSROWS_ROOT", "external/missRows")
missrows_file <- file.path(missrows_root, "data", "NCI60.rda")
if (file.exists(missrows_file)) {
  suppressWarnings(load(missrows_file))
  trans_ids <- colnames(NCI60$dataTables$trans)
  prote_ids <- colnames(NCI60$dataTables$prote)
  all_ids <- union(trans_ids, prote_ids)
  rows <- list()
  k <- 1L
  for (id in trans_ids) {
    rows[[k]] <- data.frame(
      subject_id=id, sample_id=id, assay_id=paste0("RNA_", id),
      omic="transcriptomics", condition="NCI60", timepoint="", batch="",
      technical_replicate="1", outcome="", stringsAsFactors=FALSE
    )
    k <- k + 1L
  }
  for (id in prote_ids) {
    rows[[k]] <- data.frame(
      subject_id=id, sample_id=id, assay_id=paste0("PROT_", id),
      omic="proteomics", condition="NCI60", timepoint="", batch="",
      technical_replicate="1", outcome="", stringsAsFactors=FALSE
    )
    k <- k + 1L
  }

  miss_dir <- file.path(out_dir, "missrows-nci60-partial")
  dir.create(miss_dir, recursive=TRUE, showWarnings=FALSE)
  write.csv(do.call(rbind, rows), file.path(miss_dir, "metadata.csv"), row.names=FALSE, quote=FALSE)
  write.csv(do.call(rbind, rows), file.path(out_dir, "missrows_nci60_metadata.csv"), row.names=FALSE, quote=FALSE)

  trans <- t(as.matrix(NCI60$dataTables$trans))
  prote <- t(as.matrix(NCI60$dataTables$prote))
  top_variable_columns <- function(x, n = 120L) {
    vars <- apply(x, 2, stats::var, na.rm = TRUE)
    keep <- head(order(vars, decreasing = TRUE, na.last = NA), min(n, ncol(x)))
    x[, keep, drop = FALSE]
  }
  trans <- top_variable_columns(trans)
  prote <- top_variable_columns(prote)
  write_matrix(
    trans,
    colnames(trans),
    paste0("RNA_", rownames(trans)),
    file.path(miss_dir, "transcriptomics.csv")
  )
  write_matrix(
    prote,
    colnames(prote),
    paste0("PROT_", rownames(prote)),
    file.path(miss_dir, "proteomics.csv")
  )

  cat("\nMISSROWS_NCI60\n")
  cat("transcriptomic_subjects:", length(trans_ids), "\n")
  cat("proteomic_subjects:", length(prote_ids), "\n")
  cat("union_subjects:", length(all_ids), "\n")
  cat("matched_subjects:", length(intersect(trans_ids, prote_ids)), "\n")
}


# --- Inspect public Aging HF/CD 3-omics demo objects for deterministic export ---
xomics_root <- Sys.getenv("XOMICS_ROOT", "external/xOmicsShiny")
for (fname in c("AgingHFCD_RNAseq.RData","AgingHFCD_Proteomics.RData","AgingHFCD_Metabolomics.RData")) {
  path <- file.path(xomics_root, "data", fname)
  if (file.exists(path)) {
    env <- new.env(parent=emptyenv())
    load(path, envir=env)
    cat("\nXOMICS_OBJECT", fname, "\n")
    cat("objects:", paste(ls(env), collapse=","), "\n")
    for (nm in ls(env)) {
      obj <- get(nm, envir=env)
      cat("object:", nm, "class:", paste(class(obj), collapse="/"), "\n")
      if (!is.null(dim(obj))) cat("dim:", paste(dim(obj), collapse="x"), "\n")
      if (is.list(obj)) cat("names:", paste(names(obj), collapse=","), "\n")
      if (nm == "MetaData") {
        cat("groups:", paste(names(table(obj$group)), as.integer(table(obj$group)), collapse="; "), "\n")
      }
      if (nm == "data_wide") {
        cat("rownames_head:", paste(head(rownames(obj), 8), collapse=","), "\n")
      }
      if (nm == "data_results") {
        cat("uniqueid_head:", paste(head(obj$UniqueID, 8), collapse=","), "\n")
        hits <- obj[grepl("^(Ctsd|CTSD|St7|ST7)$", obj$Gene.Name, ignore.case=TRUE), , drop=FALSE]
        if (nrow(hits)) {
          cols <- intersect(c("UniqueID","Gene.Name","Old_CDvsYoung_CD_logFC","Old_CDvsYoung_CD_P.Value","Old_CDvsYoung_CD_Adj.P.Value"), colnames(hits))
          cat("AGING_TRUTH_ROWS\n")
          print(hits[, cols, drop=FALSE])
        }
      }
    }
  }
}


# --- Export an actual 3-omics fixed-diet age contrast from AgingHFCD ---
aging_dir <- file.path(out_dir, "aging-hfcd-old-vs-young-cd")
dir.create(aging_dir, recursive = TRUE, showWarnings = FALSE)

export_xomics_layer <- function(fname, layer, assay_prefix, n_features = 80L) {
  env <- new.env(parent = emptyenv())
  load(file.path(xomics_root, "data", fname), envir = env)
  meta <- env$MetaData[env$MetaData$group %in% c("Young_CD", "Old_CD"), , drop = FALSE]
  res <- env$data_results
  qcol <- "Old_CDvsYoung_CD_Adj.P.Value"
  lfccol <- "Old_CDvsYoung_CD_logFC"

  labels <- if (layer == "metabolomics") {
    as.character(res$UniqueID)
  } else {
    gene <- as.character(res$Gene.Name)
    ifelse(!is.na(gene) & nzchar(gene), gene, as.character(res$UniqueID))
  }

  ord <- order(res[[qcol]], abs(res[[lfccol]]), na.last = NA)
  forced <- if (layer == "metabolomics") integer(0) else which(tolower(labels) %in% c("ctsd", "st7"))
  candidate_idx <- unique(c(forced, ord))
  candidate_idx <- candidate_idx[!duplicated(labels[candidate_idx])]
  selected_idx <- head(candidate_idx, n_features)
  selected_ids <- as.character(res$UniqueID[selected_idx])
  selected_labels <- labels[selected_idx]

  long <- env$data_long[
    env$data_long$UniqueID %in% selected_ids &
      env$data_long$sampleid %in% meta$sampleid,
    c("UniqueID","sampleid","expr"),
    drop = FALSE
  ]

  assay_ids <- paste0(assay_prefix, "_", meta$sampleid)
  mat <- matrix(
    NA_real_,
    nrow = length(selected_ids),
    ncol = nrow(meta),
    dimnames = list(selected_labels, assay_ids)
  )
  row_index <- match(as.character(long$UniqueID), selected_ids)
  col_index <- match(as.character(long$sampleid), meta$sampleid)
  valid <- !is.na(row_index) & !is.na(col_index)
  mat[cbind(row_index[valid], col_index[valid])] <- as.numeric(long$expr[valid])

  out <- data.frame(feature_id = rownames(mat), mat, check.names = FALSE)
  write.csv(out, file.path(aging_dir, paste0(layer, ".csv")), row.names = FALSE, quote = FALSE, na = "")

  metadata <- data.frame(
    subject_id = as.character(meta$sampleid),
    sample_id = as.character(meta$sampleid),
    assay_id = assay_ids,
    omic = layer,
    condition = as.character(meta$group),
    timepoint = "",
    batch = "",
    technical_replicate = "1",
    outcome = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  reference <- data.frame(
    layer = layer,
    feature = selected_labels,
    original_id = selected_ids,
    official_old_vs_young_logFC = as.numeric(res[[lfccol]][selected_idx]),
    official_q = as.numeric(res[[qcol]][selected_idx]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(metadata = metadata, reference = reference, n_features = nrow(mat), n_samples = ncol(mat))
}

aging_rna <- export_xomics_layer("AgingHFCD_RNAseq.RData", "transcriptomics", "RNA")
aging_pro <- export_xomics_layer("AgingHFCD_Proteomics.RData", "proteomics", "PROT")
aging_met <- export_xomics_layer("AgingHFCD_Metabolomics.RData", "metabolomics", "MET")

write.csv(
  do.call(rbind, list(aging_rna$metadata, aging_pro$metadata, aging_met$metadata)),
  file.path(aging_dir, "metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  do.call(rbind, list(aging_rna$reference, aging_pro$reference, aging_met$reference)),
  file.path(aging_dir, "reference.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nAGING_HFCD_EXPORTED\n")
cat("rna:", aging_rna$n_features, "features,", aging_rna$n_samples, "samples\n")
cat("protein:", aging_pro$n_features, "features,", aging_pro$n_samples, "samples\n")
cat("metabolite:", aging_met$n_features, "features,", aging_met$n_samples, "samples\n")


# --- LRRK2 neuron RNA/proteome public benchmark inspection ---
for (fname in c("LRRK2_Neuron_RNA.RData","LRRK2_Neuron_Protein.RData")) {
  path <- file.path(xomics_root, "data", fname)
  if (file.exists(path)) {
    env <- new.env(parent=emptyenv())
    load(path, envir=env)
    cat("\nLRRK2_PUBLIC_INSPECTION", fname, "\n")
    cat("objects:", paste(ls(env), collapse=","), "\n")
    for (nm in ls(env)) {
      obj <- get(nm, envir=env)
      cat("object:", nm, "class:", paste(class(obj), collapse="/"), "\n")
      if (!is.null(dim(obj))) cat("dim:", paste(dim(obj), collapse="x"), "\n")
      if (is.list(obj)) cat("names:", paste(names(obj), collapse=","), "\n")
      if (nm == "MetaData") {
        cat("metadata_columns:", paste(colnames(obj), collapse=","), "\n")
        for (col in intersect(c("group","condition","Genotype","genotype"), colnames(obj))) {
          cat(col, ":", paste(names(table(obj[[col]])), as.integer(table(obj[[col]])), collapse="; "), "\n")
        }
      }
      if (nm == "data_results") {
        cat("result_columns:", paste(colnames(obj), collapse=","), "\n")
        hits <- obj[grepl("LRRK2|RAB|CLTC|DNM|SH3GL|SYNJ", apply(obj, 1, paste, collapse=" "), ignore.case=TRUE), , drop=FALSE]
        if (nrow(hits)) {
          cat("LRRK2_TRUTH_ROWS\n")
          print(head(hits, 20))
        }
      }
    }
  }
}


# --- Export LRRK2 G2019S vs control at day 35 (RNA + proteome) ---
lrrk_dir <- file.path(out_dir, "lrrk2-d35")
dir.create(lrrk_dir, recursive = TRUE, showWarnings = FALSE)

export_lrrk_layer <- function(fname, layer, assay_prefix, control_group, disease_group,
                              qcol, lfccol, n_features = 100L) {
  env <- new.env(parent = emptyenv())
  load(file.path(xomics_root, "data", fname), envir = env)

  meta <- env$MetaData[env$MetaData$group %in% c(control_group, disease_group), , drop = FALSE]
  res <- as.data.frame(env$data_results, stringsAsFactors = FALSE)

  symbol <- as.character(res$Gene.Name)
  symbol[is.na(symbol) | !nzchar(symbol)] <- as.character(res$UniqueID[is.na(symbol) | !nzchar(symbol)])

  ord <- order(res[[qcol]], -abs(res[[lfccol]]), na.last = NA)
  force_names <- c("RAB29","RAB25","RAB10","RAB3A","RAB3B","CLTC","SYNJ1","SYNJ2","DNM1L","SH3GLB1","SH3GLB2")
  forced <- which(toupper(symbol) %in% force_names)
  idx <- unique(c(forced, ord))
  idx <- idx[!duplicated(symbol[idx])]
  idx <- head(idx, n_features)

  ids <- as.character(res$UniqueID[idx])
  labels <- symbol[idx]

  long <- env$data_long[
    env$data_long$UniqueID %in% ids &
      env$data_long$sampleid %in% meta$sampleid,
    c("UniqueID","sampleid","expr","group"),
    drop = FALSE
  ]

  assay_ids <- paste0(assay_prefix, "_", meta$sampleid)
  mat <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = nrow(meta),
    dimnames = list(labels, assay_ids)
  )
  row_index <- match(as.character(long$UniqueID), ids)
  col_index <- match(as.character(long$sampleid), meta$sampleid)
  valid <- !is.na(row_index) & !is.na(col_index)
  mat[cbind(row_index[valid], col_index[valid])] <- as.numeric(long$expr[valid])

  out <- data.frame(feature_id = rownames(mat), mat, check.names = FALSE)
  write.csv(out, file.path(lrrk_dir, paste0(layer, ".csv")), row.names = FALSE, quote = FALSE, na = "")

  condition <- ifelse(meta$group == control_group, "Control", "G2019S")
  metadata <- data.frame(
    subject_id = paste0(layer, "_", meta$sampleid),
    sample_id = paste0(layer, "_", meta$sampleid),
    assay_id = assay_ids,
    omic = layer,
    condition = condition,
    timepoint = "D35",
    batch = "",
    technical_replicate = "1",
    outcome = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  reference <- data.frame(
    layer = layer,
    feature = labels,
    original_id = ids,
    official_logFC_G2019S_vs_Control = as.numeric(res[[lfccol]][idx]),
    official_q = as.numeric(res[[qcol]][idx]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(metadata=metadata, reference=reference, n_features=nrow(mat), n_samples=ncol(mat))
}

lrrk_rna <- export_lrrk_layer(
  "LRRK2_Neuron_RNA.RData", "transcriptomics", "RNA",
  "Control_D35", "G2019S_D35",
  "DiseaseD35vsCtlD35_Adj.P.value", "DiseaseD35vsCtlD35_logFC", 100L
)
lrrk_pro <- export_lrrk_layer(
  "LRRK2_Neuron_Protein.RData", "proteomics", "PROT",
  "CtlD35", "DiseaseD35",
  "DiseaseD35vsCtlD35_LIMMA.Adj.P.value", "DiseaseD35vsCtlD35_logFC", 100L
)

write.csv(
  do.call(rbind, list(lrrk_rna$metadata, lrrk_pro$metadata)),
  file.path(lrrk_dir, "metadata.csv"), row.names=FALSE, quote=FALSE
)
write.csv(
  do.call(rbind, list(lrrk_rna$reference, lrrk_pro$reference)),
  file.path(lrrk_dir, "reference.csv"), row.names=FALSE, quote=FALSE
)

cat("\nLRRK2_D35_EXPORTED\n")
cat("rna:", lrrk_rna$n_features, "features,", lrrk_rna$n_samples, "samples\n")
cat("protein:", lrrk_pro$n_features, "features,", lrrk_pro$n_samples, "samples\n")


# --- TCGA breast 3-subtype benchmark for multi-group inference ---
load(file.path(mixomics_root, "data", "breast.TCGA.rda"))
tcga3 <- breast.TCGA$data.train
keep3 <- as.character(tcga3$subtype) %in% c("Basal","Her2","LumA")
mrna3 <- tcga3$mrna[keep3,,drop=FALSE]
protein3 <- tcga3$protein[keep3,,drop=FALSE]
subtype3 <- droplevels(tcga3$subtype[keep3])

tcga3_dir <- file.path(out_dir, "tcga-three-subtypes")
dir.create(tcga3_dir, recursive=TRUE, showWarnings=FALSE)
n3 <- nrow(mrna3)
subjects3 <- sprintf("TCGA3_%03d", seq_len(n3))
rna3 <- sprintf("RNA3_%03d", seq_len(n3))
prot3 <- sprintf("PROT3_%03d", seq_len(n3))

write_matrix(mrna3, colnames(mrna3), rna3, file.path(tcga3_dir, "transcriptomics.csv"))
write_matrix(protein3, colnames(protein3), prot3, file.path(tcga3_dir, "proteomics.csv"))
write_metadata(
  subjects3,
  list(transcriptomics=rna3, proteomics=prot3),
  subtype3,
  file.path(tcga3_dir, "metadata.csv")
)
cat("\nTCGA_THREE_SUBTYPES_EXPORTED\n")
cat("subjects:", n3, "\n")
cat("subtypes:", paste(names(table(subtype3)), as.integer(table(subtype3)), collapse="; "), "\n")
