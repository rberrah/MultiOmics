#!/usr/bin/env Rscript
# Real, small-sample individual-assay score from the pinned public Nutrimouse
# study. This is NOT a predictor of PK, external-cohort validation, or an
# unsupervised MOFA factor. The ONLY fitted quantities are training-subject
# mean and SD of a predeclared, measured liver gene.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)
 stop("Usage: Rscript scripts/real_nutrimouse_frozen_score.R external/mixOmics output_dir")
source_file <- file.path(args[[1L]],"data","nutrimouse.rda")
dest <- args[[2L]]
if(!file.exists(source_file) || !requireNamespace("jsonlite",quietly=TRUE))
 stop("Pinned original Nutrimouse data and jsonlite are required")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
e <- new.env(parent=emptyenv())
load(source_file,envir=e)
if(!exists("nutrimouse",envir=e,inherits=FALSE))
 stop("Nutrimouse original object absent")
study <- e$nutrimouse
if(!all(c("gene","lipid","genotype","diet") %in% names(study)))
 stop("Expected real two-omic and experimental metadata are missing")
gene <- as.matrix(study$gene)
lipid <- as.matrix(study$lipid)
storage.mode(gene) <- "double"
storage.mode(lipid) <- "double"
if(!identical(dim(gene),c(40L,120L)) ||
   !identical(dim(lipid),c(40L,21L)) ||
   any(!is.finite(gene)) || any(!is.finite(lipid)) ||
   is.null(colnames(gene)) || is.null(colnames(lipid)) ||
   anyDuplicated(colnames(gene)) || anyDuplicated(colnames(lipid)))
 stop("The frozen genuine measured Nutrimouse assay contract has changed")

# This marker is chosen BEFORE inspecting phenotype association. It is not
# selected by a significant p-value or optimized to predict genotype.
marker <- "Cyp3a11"
location <- which(tolower(colnames(gene))==tolower(marker))
if(length(location)!=1L)
 stop("Prespecified Cyp3a11 marker absent or ambiguous: do not pick a post-hoc replacement")
marker <- colnames(gene)[location]
ids <- sprintf("NUTRIMOUSE_%03d",seq_len(40L))
rownames(gene) <- rownames(lipid) <- ids
colnames(gene) <- as.character(colnames(gene))
genotype <- as.character(study$genotype)
diet <- as.character(study$diet)
if(length(genotype)!=40L || length(diet)!=40L ||
   anyNA(genotype) || anyNA(diet) ||
   length(unique(genotype))!=2L || length(unique(diet))!=5L)
 stop("Experimental genotype/diet design changed")

# Index split is predeclared, fixed and independent of BOTH genotype and diet.
# One animal = one subject; two assays from that animal are NOT two subjects.
test_index <- seq(4L,40L,by=4L)
train_index <- setdiff(seq_len(40L),test_index)
stopifnot(length(train_index)==30L,length(test_index)==10L,
          !any(ids[train_index] %in% ids[test_index]))
if(min(table(factor(genotype[test_index],levels=unique(genotype))))<2L)
 stop("Both genuine genotypes need at least two heldout subjects for descriptive audit")
values <- gene[,location]
fit_on_train <- function(v) {
  c(mean=mean(v[train_index]),sd=stats::sd(v[train_index]))
}
parameters <- fit_on_train(values)
if(any(!is.finite(parameters)) || parameters[["sd"]]<=0)
 stop("Frozen training-only scaling parameters not estimable")
# Falsification: modifying ONLY held-out observations leaves fitted model
# unchanged, even if they are extreme.
counterfactual <- values
counterfactual[test_index] <- counterfactual[test_index]*1e6+999999
if(!identical(unname(parameters),unname(fit_on_train(counterfactual))))
 stop("Heldout samples leaked into fitted scoring parameters")
heldout_z <- (values[test_index]-parameters[["mean"]])/parameters[["sd"]]
if(any(!is.finite(heldout_z)))stop("Nonfinite genuine heldout scores")

# SECOND, deliberately modest 4-gene technical composite. The marker Cyp3a11
# is prespecified. Three more genes are selected by TRAINING-ONLY variance,
# with no genotype, diet, heldout value or significance screen. It demonstrates
# frozen multi-feature projection, NOT a validated pathway/activity score.
panel_candidate <- setdiff(seq_len(ncol(gene)),location)
train_sds <- apply(gene[train_index,panel_candidate,drop=FALSE],2,stats::sd)
valid <- which(is.finite(train_sds) & train_sds>1e-12)
if(length(valid)<3L)stop("Fewer than three nonconstant TRAINING-only genes")
ranked <- valid[order(-train_sds[valid],colnames(gene)[panel_candidate[valid]],method="radix")]
panel_indices <- c(location,panel_candidate[ranked[seq_len(3L)]])
panel_names <- colnames(gene)[panel_indices]
panel_train <- gene[train_index,panel_indices,drop=FALSE]
panel_means <- colMeans(panel_train)
panel_sds <- apply(panel_train,2,stats::sd)
if(any(!is.finite(panel_sds)) || any(panel_sds<=0))
 stop("Multigene training parameters are not finite")
panel_z <- sweep(sweep(gene[test_index,panel_indices,drop=FALSE],
                       2,panel_means,"-"),2,panel_sds,"/")
heldout_panel_z <- rowMeans(panel_z)
if(any(!is.finite(heldout_panel_z)))stop("Nonfinite heldout panel score")
# Negative-control: entire heldout matrix can change without affecting
# training-only feature selection, means or SDs.
counterfactual_matrix <- gene
counterfactual_matrix[test_index,] <- counterfactual_matrix[test_index,]*1e6+999999
refit_candidate_sds <- apply(counterfactual_matrix[train_index,panel_candidate,drop=FALSE],2,stats::sd)
if(!identical(unname(train_sds),unname(refit_candidate_sds)) ||
   !identical(unname(panel_means),
              unname(colMeans(counterfactual_matrix[train_index,panel_indices,drop=FALSE]))) ||
   !identical(unname(panel_sds),
              unname(apply(counterfactual_matrix[train_index,panel_indices,drop=FALSE],2,stats::sd))))
 stop("Heldout leakage into multigene selection, center or scale")

sha <- function(path) {
  result <- suppressWarnings(system2("sha256sum",shQuote(normalizePath(path)),
                                      stdout=TRUE,stderr=TRUE))
  if(length(result)!=1L || !grepl("^[0-9a-f]{64}  ",result))
    stop(paste("SHA-256 from real file unavailable:",path))
  substr(result,1L,64L)
}
write_json <- function(value,path) {
  jsonlite::write_json(value,path,pretty=TRUE,auto_unbox=TRUE,
                       null="null",na="null",digits=15L)
}
# Original measured matrices are preserved on their declared source scales:
# gene expression is NOT RNA-seq raw counts; lipids are percentages, not
# LC-MS ion intensities or plasma drug concentrations.
rna_assays <- paste0("RNA_",sprintf("%03d",seq_len(40L)))
lipid_assays <- paste0("LIPID_",sprintf("%03d",seq_len(40L)))
# Physical aliquots were not uniquely documented in the original source.
# Keep separate specimen IDs across RNA/lipid assays, paired by animal only.
rna_specimens <- paste0("SPEC_RNA_",sprintf("%03d",seq_len(40L)))
lipid_specimens <- paste0("SPEC_LIPID_",sprintf("%03d",seq_len(40L)))
metadata <- rbind(
  data.frame(assay_id=rna_assays,sample_id=rna_specimens,subject_id=ids,
             omic="transcriptomics",tissue="liver",timepoint="",
             condition=genotype,source_dataset_id="mixOmics_nutrimouse",
             stringsAsFactors=FALSE),
  data.frame(assay_id=lipid_assays,sample_id=lipid_specimens,subject_id=ids,
             omic="metabolomics",tissue="liver",timepoint="",
             condition=genotype,source_dataset_id="mixOmics_nutrimouse",
             stringsAsFactors=FALSE)
)
write.csv(metadata,file.path(dest,"metadata.csv"),row.names=FALSE,na="NA")
# Matrix columns are ASSAY IDs, not biological sample or row index aliases.
rna_values <- t(gene)
lipid_values <- t(lipid)
colnames(rna_values) <- rna_assays
colnames(lipid_values) <- lipid_assays
write.csv(data.frame(feature_id=colnames(gene),rna_values,check.names=FALSE),
          file.path(dest,"transcriptomics.csv"),row.names=FALSE,na="NA")
write.csv(data.frame(feature_id=colnames(lipid),lipid_values,check.names=FALSE),
          file.path(dest,"metabolomics.csv"),row.names=FALSE,na="NA")
train_measurements <- data.frame(subject_id=ids[train_index],
                                 feature_id=marker,value=values[train_index])
write.csv(train_measurements,file.path(dest,"training_marker_observations.csv"),
          row.names=FALSE)
training_input_sha <- sha(file.path(dest,"training_marker_observations.csv"))
feature_definition <- list(
  feature_id=marker,source_block="nutrimouse$gene",tissue="liver",
  original_measurement_type="preprocessed_macroarray_gene_expression",
  prespecified_before_phenotype_comparison=TRUE,
  score_scale="training_centered_z_score",
  formula="(measured_assay_expression - training_mean) / training_sample_sd",
  missing_value_policy="reject_missing_no_imputation",
  cohort="real_public_nutrimouse_one_study_not_external_validation"
)
definition_file <- file.path(dest,"feature_definition.json")
write_json(feature_definition,definition_file)
model <- list(
  method="one_feature_train_only_standardization",
  feature_id=marker,
  training_subject_ids=ids[train_index],
  heldout_subject_ids=ids[test_index],
  train_n=length(train_index),heldout_n=length(test_index),
  training_mean=unname(parameters[["mean"]]),
  training_sample_sd=unname(parameters[["sd"]]),
  training_measurements_sha256=training_input_sha,
  source_rda_sha256=sha(source_file),
  pinned_source_repository="mixOmics-org/mixOmics",
  pinned_source_commit="ef3e760526623d9e91e945e4b50be48d764efc40",
  phenotype_labels_used_for_fitting=FALSE,
  heldout_values_used_for_fitting=FALSE,
  statistical_model_not_clinical_predictor=TRUE
)
model_file <- file.path(dest,"frozen_model.json")
write_json(model,model_file)
model_sha <- sha(model_file)
feature_sha <- sha(definition_file)
# Save true bytes from the real fitted panel and selected training measurements.
panel_train_rows <- data.frame(
  subject_id=rep(ids[train_index],times=length(panel_indices)),
  feature_id=rep(panel_names,each=length(train_index)),
  value=as.vector(panel_train),stringsAsFactors=FALSE)
write.csv(panel_train_rows,file.path(dest,"training_panel_observations.csv"),row.names=FALSE)
panel_feature_definition <- list(
  anchor="Cyp3a11",source_block="nutrimouse$gene",tissue="liver",
  original_measurement_type="preprocessed_macroarray_gene_expression",
  method="anchor_plus_three_highest_train_only_standard_deviations",
  selected_gene_ids=as.list(panel_names),
  selection_and_all_preprocessing_training_only=TRUE,
  score_scale="mean_of_four_training_centered_gene_z_scores",
  formula="rowMeans((gene_expression - training_mean) / training_sample_sd)",
  feature_weights=c(.25,.25,.25,.25),
  phenotype_used_for_selection=FALSE,missing_policy="reject_missing_no_imputation",
  scientific_status="technical_four_gene_composite_not_biologically_validated")
panel_def_file <- file.path(dest,"panel_feature_definition.json")
write_json(panel_feature_definition,panel_def_file)
panel_model <- list(
  method="four_gene_train_only_variance_selection_and_scaling",
  feature_ids=as.list(panel_names),
  training_subject_ids=as.list(ids[train_index]),
  heldout_subject_ids=as.list(ids[test_index]),
  training_means=unname(as.numeric(panel_means)),
  training_sample_sds=unname(as.numeric(panel_sds)),
  training_measurements_sha256=sha(file.path(dest,"training_panel_observations.csv")),
  source_rda_sha256=sha(source_file),
  pinned_source_commit="ef3e760526623d9e91e945e4b50be48d764efc40",
  phenotype_labels_used_for_selection_or_fitting=FALSE,
  heldout_values_used_for_selection_or_fitting=FALSE,
  clinical_or_PK_prediction_validated=FALSE)
panel_model_file <- file.path(dest,"frozen_panel_model.json")
write_json(panel_model,panel_model_file)
panel_model_sha <- sha(panel_model_file)
panel_feature_sha <- sha(panel_def_file)
analysis_id <- paste0("source-",substr(sha(file.path(dest,"metadata.csv")),1L,16L))
scores <- data.frame(
  subject_id=ids[test_index],sample_id=rna_specimens[test_index],
  assay_id=rna_assays[test_index],omic="transcriptomics",
  signature_id="nutrimouse_liver_Cyp3a11_trainonly_z",
  score=as.numeric(heldout_z),score_scale="training_centered_z_score",
  fit_scope="training_frozen_test_projection",
  reference_model_sha256=model_sha,feature_definition_sha256=feature_sha,
  source_exchange_analysis_id=analysis_id,
  score_provenance="native_r_individual_score",stringsAsFactors=FALSE
)
panel_scores <- scores
panel_scores$signature_id <- "nutrimouse_four_gene_trainonly_composite"
panel_scores$score <- as.numeric(heldout_panel_z)
panel_scores$score_scale <- "mean_of_four_training_centered_gene_z_scores"
panel_scores$reference_model_sha256 <- panel_model_sha
panel_scores$feature_definition_sha256 <- panel_feature_sha
scores <- rbind(scores,panel_scores)
write.csv(scores,file.path(dest,"individual_signatures.csv"),row.names=FALSE)
if(nrow(scores)!=20L ||
   anyDuplicated(paste(scores$assay_id,scores$signature_id)) ||
   length(unique(scores$subject_id))!=10L ||
   length(unique(scores$signature_id))!=2L)
 stop("Actual ten heldout animals and both frozen per-assay scores incomplete")
# An exploratory phenotype check, strictly AFTER all scores were fixed:
# no conditional diet/genotype evidence or biological superiority claim.
levels <- sort(unique(genotype))
g1 <- heldout_z[genotype[test_index]==levels[[1L]]]
g2 <- heldout_z[genotype[test_index]==levels[[2L]]]
effect <- mean(g2)-mean(g1)
set.seed(20261017L)
null_effect <- replicate(199L,{
  perm <- sample(genotype[test_index],length(test_index),replace=FALSE)
  mean(heldout_z[perm==levels[[2L]]])-mean(heldout_z[perm==levels[[1L]]])
})
report <- list(
  dataset="actual measured public Nutrimouse, 40 animals, 120 macroarray genes, 21 liver fatty acids",
  experimental_source_sha256=sha(source_file),
  genotype_descriptive_levels=levels,
  train_n=30L,heldout_n=10L,
  no_subject_leakage=TRUE,heldout_perturbation_invariance_pass=TRUE,
  model_sha256=model_sha,feature_definition_sha256=feature_sha,
  panel_model_sha256=panel_model_sha,panel_feature_definition_sha256=panel_feature_sha,
  panel_selected_genes=panel_names,panel_selection_training_only=TRUE,
  panel_heldout_scores=as.numeric(heldout_panel_z),
  panel_performance_not_assessed_as_biomarker=TRUE,
  heldout_mean_score_by_genotype=c(mean(g1),mean(g2)),
  heldout_difference_second_minus_first=unname(effect),
  exploratory_permutation_n=199L,
  exploratory_plus_one_tail_frequency=(1+sum(abs(null_effect)>=abs(effect)))/200,
  permutation_effects=as.numeric(null_effect),
  limits=c(
    "One published cohort, no independently acquired external animals.",
    "Preselected 120 expression features, not untargeted transcriptomics.",
    "Marker preselected by name; 3 additional genes selected by training variance only, NOT biology.",
    "Genotype/diet may confound; 10 heldout animals give low precision.",
    "The empirical permutation is descriptive, not confirmatory FDR control.",
    "Neither 1-gene nor 4-gene score predicts human CYP activity, PK or therapeutic outcomes.",
    "No PK measures, dosing, pharmacodynamic response or clinical interpretation."
  ),
  biological_validation_status="exploratory_real_assay_projection_not_external_confirmatory",
  pk_efficacy_qualified=FALSE
)
write_json(report,file.path(dest,"biological_exploratory_audit.json"))
cat(sprintf("REAL Nutrimouse frozen-score PASS: train=%d test=%d marker=%s model=%s no PK claim\n",
            length(train_index),length(test_index),marker,substr(model_sha,1,12)))
