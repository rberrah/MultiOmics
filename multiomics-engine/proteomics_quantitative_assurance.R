# Proteomics measurement suitability, independent of LC-MS metabolite
# feature screens. No missing-value imputation, DE test or protein ID inference.
proteomics_quantitative_assurance <- function(measurements, manifest,
  platform=c("lfq_ms","tmt","affinity"), scale=c("raw","log2"),
  missing_threshold=.20, group_imbalance_threshold=.30) {
  platform<-match.arg(platform);scale<-match.arg(scale)
  x<-as.matrix(measurements)
  m<-as.data.frame(manifest,stringsAsFactors=FALSE)
  if(!is.numeric(x)||nrow(x)<6L||ncol(x)<2L||
     is.null(rownames(x))||is.null(colnames(x))||
     anyNA(rownames(x))||anyNA(colnames(x))||
     anyDuplicated(rownames(x))||anyDuplicated(colnames(x))||
     any(!nzchar(rownames(x)))||any(!nzchar(colnames(x)))||
     any(!is.finite(x)&!is.na(x))||
     (scale=="raw"&&any(x<0,na.rm=TRUE))||
     !all(c("assay_id","subject_id","condition")%in%names(m))||
     anyNA(m[,c("assay_id","subject_id","condition")])||
     anyDuplicated(m$assay_id)||anyDuplicated(m$subject_id)||
     !identical(as.character(m$assay_id),rownames(x))||
     any(!nzchar(trimws(as.character(m$subject_id))))||
     any(!nzchar(trimws(as.character(m$condition))))||
     any(!is.finite(c(missing_threshold,group_imbalance_threshold)))||
     missing_threshold<0||missing_threshold>1||
     group_imbalance_threshold<0||group_imbalance_threshold>1)
    stop("Proteomics QA requires non-duplicated independent subjects, exact assays, declared scale/platform and valid intensities.")
  group<-as.character(m$condition)
  levels<-sort(unique(group))
  if(length(levels)!=2L||any(base::table(group)<3L))
    stop("This protein QC route requires two groups with at least three independent subjects each.")
  r1<-colMeans(is.na(x[group==levels[1L],,drop=FALSE]))
  r2<-colMeans(is.na(x[group==levels[2L],,drop=FALSE]))
  overall<-colMeans(is.na(x))
  table<-data.frame(feature_id=colnames(x),
    missing_fraction=overall,
    missing_group_1=r1,missing_group_2=r2,
    absolute_missingness_difference=abs(r1-r2),
    flag_high_missingness=overall>=missing_threshold,
    flag_group_imbalance=abs(r1-r2)>=group_imbalance_threshold,
    stringsAsFactors=FALSE)
  batch_confounded<-NA
  if("batch"%in%names(m)) {
    b<-as.character(m$batch)
    if(anyNA(b)||any(!nzchar(trimws(b))))
      stop("Known proteomics batches cannot have missing labels.")
    if(length(unique(b))>1L){
      design<-stats::model.matrix(~factor(group)+factor(b))
      batch_confounded<-qr(design)$rank<ncol(design)
    } else batch_confounded<-FALSE
  }
  flags<-character()
  if(any(table$flag_high_missingness))
    flags<-c(flags,"protein_missing_fraction_high")
  if(any(table$flag_group_imbalance))
    flags<-c(flags,"protein_missingness_differs_by_group")
  if(isTRUE(batch_confounded))
    flags<-c(flags,"protein_batch_confounded_with_group")
  platform_checks<-switch(platform,
    lfq_ms=c("peptide_and_protein_inference_not_verified",
             "missing_not_at_random_left_censoring_not_excluded"),
    tmt=c("reference_channel_and_ratio_compression_not_verified",
          "isolation_interference_not_verified"),
    affinity=c("target_binding_specificity_and_panel_version_not_verified",
               "assay_limit_of_detection_not_verified"))
  flags<-unique(c(flags,platform_checks))
  list(
    status="method_specific_scientific_review_required",
    platform=platform,declared_scale=scale,
    sample_count=nrow(x),independent_subjects=length(unique(m$subject_id)),
    measured_features=ncol(x),features=table,
    completely_missing_proteins=sum(overall==1),
    batch_condition_fully_confounded=batch_confounded,
    alerts=flags,
    mnar_ruled_out=FALSE,
    feature_exclusions_automatically_applied=FALSE,
    statistical_significance_computed=FALSE,
    scientific_certified=FALSE,
    limitation="Platform-specific quality and biological replicate assumptions require human examination; no imputation or target identification is inferred."
  )
}
