#!/usr/bin/env Rscript
# Cross-language scientific contract: the JS browser's OLS-HC3 group estimates
# must agree with a fresh calculation from base R on the SAME fixed samples.
# Does not require limma/DESeq2 packages and does not assert broad validation.
stopifnot(nzchar(Sys.which("node")))
fetch_csv <- function(mode) {
  stdout <- suppressWarnings(system2("node", c("scripts/test_multiomics_statistical_calibration.mjs", mode), stdout=TRUE))
  if (!length(stdout) || any(grepl("^Error:", stdout))) stop("Node reference fixture failed")
  utils::read.csv(text=paste(stdout,collapse="\n"), stringsAsFactors=FALSE, check.names=FALSE)
}
dat <- fetch_csv("--reference-input")
got <- fetch_csv("--reference-results")
stopifnot(nrow(dat) == 24L * 24L, nrow(got)==24L)
fit_reference <- function(df) {
  design <- stats::model.matrix(~ condition, data=df)
  y <- df$value
  fit <- stats::lm.fit(design, y)
  xtxinv <- solve(crossprod(design))
  residual <- as.vector(y - design %*% fit$coefficients)
  hat <- rowSums((design %*% xtxinv) * design)
  stopifnot(all(hat < 0.999))
  meat <- crossprod(design, design * (residual/(1-hat))^2)
  vcovHC3 <- xtxinv %*% meat %*% xtxinv
  standardError <- sqrt(vcovHC3[2,2])
  statistic <- fit$coefficients[2]/standardError
  c(effect=fit$coefficients[2], se=standardError,
    p_value=2*stats::pt(-abs(statistic), df=nrow(design)-ncol(design)))
}
reference <- t(vapply(split(dat, dat$feature), fit_reference, numeric(3)))
reference <- reference[match(got$feature,rownames(reference)),,drop=FALSE]
stopifnot(all(!is.na(reference)))
reference_q <- stats::p.adjust(reference[,3],method="BH")
stopifnot(max(abs(got$effect-reference[,1]))<1e-6)
stopifnot(max(abs(got$se-reference[,2]))<1e-5)
stopifnot(max(abs(got$p_value-reference[,3]))<1e-5)
stopifnot(max(abs(got$q_value-reference_q))<1e-5)
cat("Browser adjusted 2-group OLS-HC3 vs independent base-R lm / sandwich / BH: PASS,",
    nrow(got),"features\n")
