# Independent feature-coordinate matching for prespecified xcms sensitivity.
# Feature IDs are NEVER assumed stable across independent peak-picking runs.
# Reciprocal and UNIQUE neighbours only: ambiguous chromatographic evidence
# must not be counted as reproducible.
ms_feature_coordinates <- function(definitions, expected_ids) {
  if(!is.data.frame(definitions) || !length(expected_ids) ||
     anyNA(expected_ids) || anyDuplicated(expected_ids) ||
     !all(c("feature_id","mzmed","rtmed")%in%names(definitions)) ||
     nrow(definitions)!=length(expected_ids) ||
     anyDuplicated(definitions$feature_id) ||
     !setequal(as.character(definitions$feature_id),as.character(expected_ids)))
    stop("Native xcms feature coordinates do not agree with complete exported IDs.")
  x <- definitions[match(expected_ids,definitions$feature_id),,
                   drop=FALSE]
  x$feature_id <- as.character(x$feature_id)
  if(!is.numeric(x$mzmed) || !is.numeric(x$rtmed) ||
     anyNA(x[,c("mzmed","rtmed")]) ||
     any(!is.finite(x$mzmed)) || any(x$mzmed<=0) ||
     any(!is.finite(x$rtmed)) || any(x$rtmed<0))
    stop("Non-finite or invalid native xcms median m/z or RT; do not invent coordinates.")
  x[,c("feature_id","mzmed","rtmed"),drop=FALSE]
}

ms_reciprocal_unique_matches <- function(reference,alternative,
                                         ppm_tolerance=20,
                                         rt_seconds_tolerance=30) {
  validate <- function(x) {
    if(!is.data.frame(x)||!nrow(x)||
       !all(c("feature_id","mzmed","rtmed")%in%names(x)) ||
       anyNA(x[,c("feature_id","mzmed","rtmed")]) ||
       anyDuplicated(x$feature_id) ||
       any(!nzchar(as.character(x$feature_id))) ||
       !is.numeric(x$mzmed) || !is.numeric(x$rtmed) ||
       any(!is.finite(x$mzmed)) || any(x$mzmed<=0) ||
       any(!is.finite(x$rtmed)) || any(x$rtmed<0))
      stop("Feature-matching input must contain genuine finite unique coordinates.")
  }
  validate(reference)
  validate(alternative)
  if(length(ppm_tolerance)!=1L||!is.finite(ppm_tolerance)||
     ppm_tolerance<=0||ppm_tolerance>100 ||
     length(rt_seconds_tolerance)!=1L||
     !is.finite(rt_seconds_tolerance)||
     rt_seconds_tolerance<=0||rt_seconds_tolerance>300)
    stop("Prespecified reciprocal match tolerances are invalid.")
  candidate <- function(from,to) {
    # Use the query feature's mass in ppm denominator (correct query direction).
    # A sole unique match in *both* directions is needed to retain a pair.
    lapply(seq_len(nrow(from)),function(i) {
      ppm <- abs(to$mzmed-from$mzmed[[i]])/from$mzmed[[i]]*1e6
      rt <- abs(to$rtmed-from$rtmed[[i]])
      which(ppm<=ppm_tolerance & rt<=rt_seconds_tolerance)
    })
  }
  fwd <- candidate(reference,alternative)
  rev <- candidate(alternative,reference)
  result <- data.frame(
    reference_feature=as.character(reference$feature_id),
    alternative_feature=rep(NA_character_,nrow(reference)),
    candidate_count=vapply(fwd,length,integer(1)),
    match_status=rep("no_coordinate_candidate",nrow(reference)),
    matched_mass_error_ppm=rep(NA_real_,nrow(reference)),
    matched_rt_error_seconds=rep(NA_real_,nrow(reference)),
    stringsAsFactors=FALSE
  )
  for(i in seq_len(nrow(reference))) {
    hit <- fwd[[i]]
    if(length(hit)>1L) {
      result$match_status[[i]] <- "ambiguous_multiple_candidates"
      next
    }
    if(!length(hit)) next
    j <- hit[[1L]]
    if(length(rev[[j]])!=1L || rev[[j]][[1L]]!=i) {
      result$match_status[[i]] <- "not_reciprocally_unique"
      next
    }
    result$alternative_feature[[i]] <- as.character(alternative$feature_id[[j]])
    result$match_status[[i]] <- "reciprocal_unique"
    result$matched_mass_error_ppm[[i]] <-
      1e6*(alternative$mzmed[[j]]-reference$mzmed[[i]])/
      reference$mzmed[[i]]
    result$matched_rt_error_seconds[[i]] <-
      alternative$rtmed[[j]]-reference$rtmed[[i]]
  }
  if(anyDuplicated(result$alternative_feature[!is.na(result$alternative_feature)]))
    stop("Ambiguous duplicate alternative mapping detected.")
  result
}

ms_paired_feature_sensitivity <- function(matches,ref_features,alt_features) {
  if(!is.data.frame(matches)||!nrow(matches)||
     !is.data.frame(ref_features)||!is.data.frame(alt_features) ||
     !all(c("feature_id","pooled_qc_rsd",
            "blank_over_biological_ratio")%in%names(ref_features)) ||
     !all(c("feature_id","pooled_qc_rsd",
            "blank_over_biological_ratio")%in%names(alt_features)) ||
     !setequal(as.character(matches$reference_feature),
                   as.character(ref_features$feature_id)) ||
     anyDuplicated(ref_features$feature_id) ||
     anyDuplicated(alt_features$feature_id))
    stop("Paired QC sensitivity requires complete valid feature provenance.")
  ix <- which(matches$match_status=="reciprocal_unique")
  if(!length(ix)) {
    return(list(pairs=0L,both_qc_evaluable=0L,
                both_blank_evaluable=0L,qc_screen_agreement=NA_real_,
                blank_screen_agreement=NA_real_))
  }
  id <- matches$reference_feature[ix]
  jd <- matches$alternative_feature[ix]
  a <- ref_features[match(id,ref_features$feature_id),,drop=FALSE]
  z <- alt_features[match(jd,alt_features$feature_id),,drop=FALSE]
  if(anyNA(a$feature_id)||anyNA(z$feature_id))
    stop("QC feature ID pairing is not exact.")
  get_screen <- function(key,cutoff) {
    av <- a[[key]]
    zv <- z[[key]]
    if(!is.numeric(av)||!is.numeric(zv)||
       any(is.infinite(av))||any(is.infinite(zv))||
       any(av<0,na.rm=TRUE)||any(zv<0,na.rm=TRUE))
      stop("Invalid instrument QA metric or missing QC provenance.")
    k <- is.finite(av)&is.finite(zv)
    list(n=sum(k),agreement=if(sum(k)>0L)
      mean((av[k]>cutoff)==(zv[k]>cutoff)) else NA_real_)
  }
  qc <- get_screen("pooled_qc_rsd",.30)
  blank <- get_screen("blank_over_biological_ratio",.20)
  list(pairs=length(ix),both_qc_evaluable=qc$n,
       both_blank_evaluable=blank$n,
       qc_screen_agreement=qc$agreement,
       blank_screen_agreement=blank$agreement)
}
