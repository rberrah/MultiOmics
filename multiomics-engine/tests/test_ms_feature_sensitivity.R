#!/usr/bin/env Rscript
# Package-free adversarial scientific feature-matching test. Synthetic coordinates
# only; NOT external biological or instrument validation.
source("multiomics-engine/ms_feature_sensitivity.R")
ref <- data.frame(feature_id=c("F1","F2","F3","F4"),
                  mzmed=c(100,150,200,250),
                  rtmed=c(60,120,180,240))
alt <- data.frame(feature_id=c("X9","X2","X3","X4"),
                  mzmed=c(100.001,150.001,202,250.001),
                  rtmed=c(61,121,181,241))
matched <- ms_reciprocal_unique_matches(ref,alt)
stopifnot(nrow(matched)==4L,
  sum(matched$match_status=="reciprocal_unique")==3L,
  is.na(matched$alternative_feature[matched$reference_feature=="F3"]),
  !anyDuplicated(na.omit(matched$alternative_feature)),
  all(abs(matched$matched_mass_error_ppm[matched$match_status==
    "reciprocal_unique"])<20L))
# Splitting a reference feature into two apparent peaks must NOT be counted as
# biologically stable simply by taking the first nearest coordinate.
amb <- rbind(alt, data.frame(feature_id="X5",mzmed=100.0011,rtmed=60.8))
result <- ms_reciprocal_unique_matches(ref,amb)
stopifnot(result$match_status[[1L]]=="ambiguous_multiple_candidates",
  is.na(result$alternative_feature[[1L]]))
# Test reverse-only nonunique: two independent ref peaks mapping to one alt
# peak cannot both be counted as a robust feature.
twins <- rbind(ref,data.frame(feature_id="F5",mzmed=100.0015,rtmed=60.5))
result <- ms_reciprocal_unique_matches(twins,alt)
stopifnot(result$match_status[[1L]]=="not_reciprocally_unique",
          result$match_status[[5L]]=="not_reciprocally_unique")
# Reject corrupt coordinates, missing IDs, invalid ppm/RT.
reject <- function(a,b=alt,...) inherits(try(
  ms_reciprocal_unique_matches(a,b,...),silent=TRUE),"try-error")
x <- ref;x$mzmed[[1L]] <- NA_real_;stopifnot(reject(x))
x <- ref;x$mzmed[[1L]] <- 0;stopifnot(reject(x))
x <- ref;x$rtmed[[1L]] <- -1;stopifnot(reject(x))
x <- ref;x$feature_id[[1L]] <- "F2";stopifnot(reject(x))
stopifnot(reject(ref,ppm_tolerance=-1),
          reject(ref,rt_seconds_tolerance=0),
          reject(ref,rt_seconds_tolerance=Inf))
# Reordering rows with their coordinates preserves feature identity.
wrong <- ref[c(2L,1L,3L,4L),,drop=FALSE]
coordinates <- ms_feature_coordinates(wrong,ref$feature_id)
stopifnot(identical(coordinates$feature_id,ref$feature_id),
  identical(coordinates$mzmed,ref$mzmed))
wrong$feature_id[[1L]] <- "F1"
stopifnot(inherits(try(ms_feature_coordinates(
  wrong,ref$feature_id),silent=TRUE),"try-error"))
# Paired QC flags are DESCRIPTIVE, and a missing blank cannot be a clean blank.
refq <- data.frame(
  feature_id=ref$feature_id,
  pooled_qc_rsd=c(.1,.5,NA,.2),
  blank_over_biological_ratio=c(.1,.4,NA,NA))
altq <- data.frame(
  feature_id=alt$feature_id,
  pooled_qc_rsd=c(.2,.2,.3,.2),
  blank_over_biological_ratio=c(.1,.5,.1,NA))
compare <- ms_paired_feature_sensitivity(matched,refq,altq)
stopifnot(compare$pairs==3L,compare$both_qc_evaluable==3L,
  compare$both_blank_evaluable==2L,
  identical(compare$qc_screen_agreement,2/3),
  compare$blank_screen_agreement==1)
fake <- refq;fake$pooled_qc_rsd[[1L]] <- Inf
stopifnot(inherits(try(ms_paired_feature_sensitivity(matched,
  fake,altq),silent=TRUE),"try-error"))
missing <- matched
missing$match_status[] <- "no_coordinate_candidate"
empty <- ms_paired_feature_sensitivity(missing,refq,altq)
stopifnot(empty$pairs==0L,is.na(empty$qc_screen_agreement))
cat("PASS: strict reciprocal unique xcms m/z-RT matching, ambiguous peaks and paired RSD/blank audit; scientific claims remain descriptive\n")
