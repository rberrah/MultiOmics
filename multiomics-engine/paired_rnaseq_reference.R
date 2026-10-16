# Subject-blocked RNA-seq reference differential expression.
# "paired" = two unique biological samples per subject (same individual or
# same independently derived cell line), ONE in each condition. It does NOT
# describe technical replicates or two independent patients.
#
# Guard before DESeq2/limma-voom, after exact assay-input checks and aggregation.
# Never silently downgrade an invalid paired study to an independent design.
paired_rnaseq_design <- function(counts,metadata,protocol,covariates=character()) {
  if(!is.matrix(counts)||!is.numeric(counts)||ncol(counts)<20L||
     is.null(rownames(counts))||is.null(colnames(counts))||
     anyNA(rownames(counts))||anyNA(colnames(counts))||
     any(!nzchar(rownames(counts)))||any(!nzchar(colnames(counts)))||
     anyDuplicated(rownames(counts))||anyDuplicated(colnames(counts))||
     any(!is.finite(counts))||any(counts<0)||
     any(abs(counts-round(counts))>1e-8)||any(rowSums(counts)<=0))
    stop("Paired RNA-seq requires genuine finite nonnegative integer counts, unique measured samples and >=20 features.",call.=FALSE)
  if(!is.data.frame(metadata)||
     !all(c("subject_id","sample_id","condition","batch")%in%names(metadata))||
     is.null(rownames(metadata))||anyDuplicated(rownames(metadata))||
     nrow(metadata)!=nrow(counts)||!setequal(rownames(metadata),rownames(counts)))
    stop("Paired RNA-seq sample/metadata identifiers do not match exactly.",call.=FALSE)
  m<-metadata[match(rownames(counts),rownames(metadata)),,drop=FALSE]
  for(field in c("subject_id","sample_id","condition")) {
    a<-as.character(m[[field]])
    if(anyNA(a)||any(!nzchar(trimws(a))))
      stop(paste("Paired RNA-seq missing required",field,"identifier."),call.=FALSE)
  }
  if(anyDuplicated(as.character(m$sample_id))||
     !identical(as.character(m$sample_id),rownames(counts)))
    stop("Paired RNA-seq requires one unique biological sample identifier per measured assay after technical aggregation.",call.=FALSE)
  groups<-sort(unique(as.character(m$condition)))
  if(length(groups)!=2L)
    stop("Paired RNA-seq requires exactly two measured conditions.",call.=FALSE)
  subjects<-sort(unique(as.character(m$subject_id)))
  if(length(subjects)<4L)
    stop("Paired RNA-seq requires >=4 independent subjects or biological cell lines with both conditions; this is not a statistical power guarantee.",call.=FALSE)
  observed<-table(factor(as.character(m$subject_id),levels=subjects),
                  factor(as.character(m$condition),levels=groups))
  if(!identical(dim(observed),c(length(subjects),2L))||
     any(observed!=1L)||nrow(m)!=2L*length(subjects))
    stop("Paired RNA-seq requires exactly one measured biological sample per subject and condition; missing or duplicated pairs are not silently removed.",call.=FALSE)

  m$subject_id<-factor(as.character(m$subject_id),levels=subjects)
  m$condition<-factor(as.character(m$condition),levels=groups)
  # Paired identity is the first blocking effect, primary condition LAST. This
  # makes a limma-voom last coefficient the planned treatment contrast.
  nuisance<-setdiff(varying_terms(m,covariates),c("subject_id","condition"))
  if(length(nuisance)){
    for(k in nuisance){
      v<-as.character(m[[k]])
      if(anyNA(v)||any(!nzchar(trimws(v))))
        stop(paste("Paired RNA-seq cannot fit a missing or incomplete nuisance covariate:",k),call.=FALSE)
    }
  }
  formula<-safe_formula(c("subject_id",nuisance,"condition"))
  design<-stats::model.matrix(formula,data=m)
  if(nrow(design)!=nrow(counts)||qr(design)$rank!=ncol(design)||
     nrow(design)-ncol(design)<2L)
    stop("Paired RNA-seq design is confounded or has fewer than two residual degrees of freedom; no reference fit attempted.",call.=FALSE)
  # An exploratory result reports its sign explicitly (alphabetical group
  # reference). A confirmatory fit requires the researcher to declare BOTH
  # levels and their order; no implicit contrast is acceptable.
  declared<-identical(protocol$analysisIntent,"confirmatory")
  contrast<-if(declared){
    rnaseq_confirmatory_explicit_contrast(protocol,groups)
  } else c("condition",groups[[2L]],groups[[1L]])
  contract<-list(
    status="eligible_for_reference_fitting_only",
    pairing="same_subject_two_conditions",
    design=deparse(formula,width.cutoff=500L),
    model_samples=nrow(m),
    independent_biological_subjects=length(subjects),
    biological_pairs=length(subjects),
    measured_samples_per_condition=as.list(table(m$condition)),
    residual_degrees_of_freedom=nrow(design)-ncol(design),
    explicit_primary_contrast=as.list(unname(contrast)),
    contrast_was_prespecified=declared,
    batch_and_covariate_terms=as.list(nuisance),
    sample_identifiers_exact=TRUE,
    all_subject_condition_cells_observed_once=TRUE,
    certified=FALSE,
    limitations=c(
      "Two assays from one subject are not independent subjects.",
      "Four pairs is a software eligibility guard, not a power guarantee.",
      "Only an additive within-subject effect is tested; unmodeled condition-by-subject interactions remain a limitation.",
      "Independent biological and FDR validation are still required."
    )
  )
  list(metadata=m,design_formula=formula,contrast=unname(contrast),contract=contract)
}
