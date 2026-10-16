# Conservative, instrument-agnostic quantitative MS QA.
# Summarises predeclared sample types and raw peak areas without imputing
# censored signals, removing unknown features or making pathway claims.
ms_quantitative_assurance <- function(intensities, manifest,
  qc_rsd_threshold=.30, blank_fraction_threshold=.20,
  max_missing_fraction=.20, imbalance_threshold=.30) {
  x <- as.matrix(intensities)
  meta <- as.data.frame(manifest,stringsAsFactors=FALSE)
  needed <- c("assay_id","sample_type","condition","injection_order")
  thresholds <- c(qc_rsd_threshold,blank_fraction_threshold,
                  max_missing_fraction,imbalance_threshold)
  if(!is.numeric(x)||!nrow(x)||ncol(x)<2L||
     is.null(rownames(x))||is.null(colnames(x))||
     anyNA(rownames(x))||anyNA(colnames(x))||
     anyDuplicated(rownames(x))||anyDuplicated(colnames(x))||
     any(!nzchar(rownames(x)))||any(!nzchar(colnames(x)))||
     any(is.infinite(x))||any(x<0,na.rm=TRUE)||
     any(!needed%in%names(meta))||
     anyNA(meta[,c('assay_id','sample_type','injection_order')])||
     anyDuplicated(meta$assay_id)||
     !identical(as.character(meta$assay_id),rownames(x))||
     any(!is.finite(thresholds))||
     any(thresholds<0)||any(thresholds>1))
    stop("MS QC needs aligned unique raw peak-area assays, known sample roles and valid thresholds.")
  roles <- tolower(trimws(as.character(meta$sample_type)))
  allowed <- c("biological","pooled_qc","blank")
  if(any(!roles%in%allowed))
    stop("Unknown MS sample_type: explicit biological, pooled_qc or blank required.")
  biological <- which(roles=="biological")
  qc <- which(roles=="pooled_qc")
  blanks <- which(roles=="blank")
  if(length(biological)<4L)stop("At least four biological assays required for descriptive MS QA.")
  if(anyNA(meta$condition[biological]))
    stop('Biological samples have missing condition annotations.')
  groups <- trimws(as.character(meta$condition[biological]))
  if(any(!nzchar(groups)))stop("Biological samples lack group declaration.")
  order <- suppressWarnings(as.numeric(meta$injection_order))
  if(any(!is.finite(order)) || anyDuplicated(order))
    stop("Unique finite injection orders are needed for drift review.")
  tested <- x[biological,,drop=FALSE]
  missing <- colMeans(is.na(tested))
  biological_finite <- colSums(is.finite(tested))
  pooled_qc_finite <- if(length(qc)) colSums(is.finite(x[qc,,drop=FALSE]))
    else rep(0L,ncol(x))
  blank_finite <- if(length(blanks)) colSums(is.finite(x[blanks,,drop=FALSE]))
    else rep(0L,ncol(x))

  # A pooled-QC series can exist while particular features are absent from
  # most/all injections. Never label that feature as QC-passing.
  qc_rsd <- rep(NA_real_,ncol(x))
  if(length(qc)>=3L) {
    qc_rsd <- vapply(seq_len(ncol(x)),function(j) {
      z<-x[qc,j]
      z<-z[is.finite(z)]
      if(length(z)<3L || mean(z)<=0) return(NA_real_)
      stats::sd(z)/mean(z)
    },numeric(1))
  }
  # The first branch is scalar. Use if/rep rather than ifelse(scalar,...):
  # scalar ifelse silently takes only the first feature's status and recycles
  # it to every molecule, hiding later missing QC series.
  qc_status <- if(length(qc)<3L)
    rep("too_few_pooled_qc_injections",ncol(x))
  else ifelse(pooled_qc_finite<3L,"too_few_finite_pooled_qc_values",
      ifelse(is.finite(qc_rsd),"rsd_evaluable",
             "pooled_qc_mean_zero_or_not_estimable"))

  blank_fraction <- rep(NA_real_,ncol(x))
  blank_status <- rep("not_evaluable",ncol(x))
  if(length(blanks)<2L) {
    blank_status[] <- "too_few_blank_injections"
  } else for(j in seq_len(ncol(x))) {
    if(blank_finite[[j]]<2L) {
      blank_status[[j]] <- "too_few_finite_blank_values"
      next
    }
    if(biological_finite[[j]]<2L) {
      blank_status[[j]] <- "too_few_finite_biological_values"
      next
    }
    numerator <- stats::median(x[blanks,j],na.rm=TRUE)
    denominator <- stats::median(tested[,j],na.rm=TRUE)
    if(denominator>0) {
      blank_fraction[[j]]<-numerator/denominator
      blank_status[[j]]<-"ratio_evaluable"
    } else if(numerator>0) {
      # This is positive solvent contamination against a biological median
      # of zero. NA would suppress the critical contamination flag.
      blank_status[[j]]<-"positive_blank_with_zero_biological_median"
    } else blank_status[[j]]<-"both_blank_and_biological_medians_zero"
  }

  strata <- sort(unique(groups))
  diff_missing <- rep(NA_real_,ncol(x))
  if(length(strata)==2L) {
    rates <- vapply(strata,function(g)
      colMeans(is.na(tested[groups==g,,drop=FALSE])),numeric(ncol(x)))
    diff_missing<-abs(rates[,1L]-rates[,2L])
  }

  drift <- rep(NA_real_,ncol(x))
  drift_status <- rep("too_few_pooled_qc_injections",ncol(x))
  if(length(qc)>=5L) for(j in seq_len(ncol(x))) {
    v<-x[qc,j]
    ok<-is.finite(v)
    if(sum(ok)<5L) {
      drift_status[[j]]<-"too_few_finite_pooled_qc_values"
    } else if(stats::sd(v[ok])==0) {
      drift_status[[j]]<-"constant_pooled_qc_signal_correlation_undefined"
    } else {
      drift[[j]]<-suppressWarnings(stats::cor(order[qc][ok],v[ok],
                                              method="spearman"))
      drift_status[[j]]<-"screen_evaluable"
    }
  }
  feature_table <- data.frame(
    feature_id=colnames(x),
    finite_biological_values=as.integer(biological_finite),
    finite_pooled_qc_values=as.integer(pooled_qc_finite),
    finite_blank_values=as.integer(blank_finite),
    biological_missing_fraction=missing,
    pooled_qc_rsd=qc_rsd,
    pooled_qc_precision_status=qc_status,
    blank_over_biological_ratio=blank_fraction,
    blank_background_status=blank_status,
    missingness_group_difference=diff_missing,
    pooled_qc_order_spearman=drift,
    pooled_qc_drift_status=drift_status,
    flag_missing=missing>max_missing_fraction,
    flag_qc_variation=is.finite(qc_rsd)&qc_rsd>qc_rsd_threshold,
    flag_qc_not_evaluable=(length(qc)>=3L)&!is.finite(qc_rsd),
    flag_blank_signal=(is.finite(blank_fraction)&
                         blank_fraction>=blank_fraction_threshold)|
      blank_status=="positive_blank_with_zero_biological_median",
    flag_blank_not_evaluable=(length(blanks)>=2L)&
      blank_status!="ratio_evaluable" &
      blank_status!="positive_blank_with_zero_biological_median",
    flag_group_missingness=is.finite(diff_missing)&
      diff_missing>=imbalance_threshold,
    stringsAsFactors=FALSE)
  batch_confounded <- NA
  if("batch"%in%names(meta)){
    batches <- trimws(as.character(meta$batch[biological]))
    if(anyNA(batches)||any(!nzchar(batches)))
      stop("Batch is declared but incomplete among biological MS assays.")
    if(length(unique(groups))==2L && length(unique(batches))>=2L) {
      # Perfect confounding prevents distinguishing biology from instrument.
      design <- stats::model.matrix(~factor(batches)+factor(groups))
      batch_confounded <- qr(design)$rank<ncol(design)
    } else batch_confounded <- FALSE
  }
  coverage <- list(
    biological_assays=length(biological),
    pooled_qc_assays=length(qc),
    blank_assays=length(blanks),
    declared_groups=as.list(base::table(groups)),
    batch_condition_fully_confounded=batch_confounded,
    biological_missingness_evaluable=TRUE,
    pooled_qc_precision_evaluable=length(qc)>=3L,
    blank_background_evaluable=length(blanks)>=2L,
    pooled_qc_drift_screen_evaluable=length(qc)>=5L,
    missingness_group_imbalance_evaluable=length(strata)==2L)
  alerts <- character()
  if(isTRUE(batch_confounded))alerts<-c(alerts,"batch_perfectly_confounded_with_biological_group")
  if(!coverage$pooled_qc_precision_evaluable)
    alerts <- c(alerts,"pooled_qc_not_sufficient_for_rsd")
  if(!coverage$blank_background_evaluable)
    alerts <- c(alerts,"blanks_not_sufficient_for_background")
  if(!coverage$pooled_qc_drift_screen_evaluable)
    alerts <- c(alerts,"pooled_qc_not_sufficient_for_drift")
  if(!coverage$missingness_group_imbalance_evaluable)
    alerts <- c(alerts,"exactly_two_biological_groups_needed_for_missingness_comparison")
  for(flag in c("flag_missing","flag_qc_variation","flag_qc_not_evaluable",
                "flag_blank_signal","flag_blank_not_evaluable",
                "flag_group_missingness")){
    if(any(feature_table[[flag]]))
      alerts <- c(alerts,paste0(flag,"_",sum(feature_table[[flag]])))
  }
  list(
    status=if(length(alerts))"scientific_review_required" else
      "measured_QC_screens_pass_not_certified",
    coverage=coverage,features=feature_table,
    alerts=alerts,
    parameters=list(qc_rsd_threshold=qc_rsd_threshold,
      blank_fraction_threshold=blank_fraction_threshold,
      max_missing_fraction=max_missing_fraction,
      imbalance_threshold=imbalance_threshold),
    feature_exclusions_automatically_applied=FALSE,
    per_feature_qc_observability_recorded=TRUE,
    mnar_diagnosis_possible_from_this_table=FALSE,
    pathway_or_compound_identity_certified=FALSE,
    scientific_certified=FALSE,
    limitation=paste("Peak-area QC is descriptive; zero and NA differ;",
      "an absence of group missingness imbalance does not rule out MNAR.",
      "A measured QC/blank injection is not evidence for a particular missing or zero-valued feature.",
      "No drift correction or biological feature exclusion is selected after seeing p-values.")
  )
}
