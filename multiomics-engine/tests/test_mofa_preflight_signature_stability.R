#!/usr/bin/env Rscript
# Pure R negative-control contracts; no MOFA2/mixOmics dependency required.
source("multiomics-engine/advanced_methods.R")
set.seed(20261009)
ids <- sprintf("S%03d",seq_len(48L))
a <- matrix(stats::rnorm(48*14),48,14,
  dimnames=list(ids,paste0("G",seq_len(14))))
b <- matrix(stats::rnorm(48*11),48,11,
  dimnames=list(rev(ids),paste0("P",seq_len(11))))
# Identical subject sets in arbitrary order are legal.
p <- mofa2_input_preflight(list(transcriptomics=a,proteomics=b))
stopifnot(
  identical(p$status,"eligible_for_exploration"),
  identical(p$nSharedSubjects,48L),
  length(p$views)==2L,
  isTRUE(p$views$transcriptomics$missingFraction==0)
)
miss <- b
rownames(miss)[1] <- "UNKNOWN_SAMPLE"
err <- try(mofa2_input_preflight(list(transcriptomics=a,proteomics=miss)),silent=TRUE)
stopifnot(inherits(err,"try-error"),
  grepl("different sample sets",as.character(err),fixed=TRUE))
duplicate <- b
rownames(duplicate)[1] <- rownames(duplicate)[2]
stopifnot(inherits(try(mofa2_input_preflight(
  list(transcriptomics=a,proteomics=duplicate)),silent=TRUE),"try-error"))
too_missing <- a
too_missing[1:32,1:10] <- NA_real_
stopifnot(inherits(try(mofa2_input_preflight(
  list(transcriptomics=too_missing,proteomics=b)),silent=TRUE),"try-error"))
bad_inf <- a
bad_inf[2,2] <- Inf
stopifnot(inherits(try(mofa2_input_preflight(
  list(transcriptomics=bad_inf,proteomics=b)),silent=TRUE),"try-error"))
bad_features <- a
colnames(bad_features)[1] <- colnames(bad_features)[2]
stopifnot(inherits(try(mofa2_input_preflight(
  list(transcriptomics=bad_features,proteomics=b)),silent=TRUE),"try-error"))
# Run-MOFA2 must fail on eligibility BEFORE attempting package installation.
err <- try(run_mofa2_blocks(list(transcriptomics=a,proteomics=miss),
  tempfile("mofa_preflight_")),silent=TRUE)
stopifnot(inherits(err,"try-error"),
  grepl("different sample sets",as.character(err),fixed=TRUE))

# Signature stability is a descriptive property of selected sets, not
# an automatically valid biomarker inference.
sel <- list(
  list(rna=c("G1","G2","G3"),prot=c("P1","P2")),
  list(rna=c("G1","G2","G4"),prot=c("P1","P3")),
  list(rna=c("G1","G5","G3"),prot=c("P1","P2"))
)
stability <- diablo_selection_stability(sel)
stopifnot(
  isTRUE(all.equal(stability$rna$meanPairwiseJaccard,
    mean(c(2/4,2/4,1/5)),tolerance=1e-12)),
  stability$prot$frequency$P1==1,
  stability$rna$frequency$G1==1,
  setequal(stability$rna$repeatedlySelected,c("G1","G2","G3")),
  identical(stability$rna$nPairwiseComparisons,3L)
)
stopifnot(inherits(try(diablo_selection_stability(
 list(sel[[1]],list(rna="G1",wrong="P1"))),silent=TRUE),"try-error"))
stopifnot(inherits(try(diablo_selection_stability(
 list(sel[[1]],list(rna=character(),prot="P1"))),silent=TRUE),"try-error"))
cat("MOFA2 input integrity and no hidden sample intersection; DIABLO stability Jaccard/frequency math and negative controls: PASS\n")
