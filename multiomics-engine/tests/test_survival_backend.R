#!/usr/bin/env Rscript
# Reference survival routing guardrails and feature-wise full-family BH check.
source("multiomics-engine/advanced_methods.R")
if(!requireNamespace("survival",quietly=TRUE))stop("survival R package required")
set.seed(20261009)
n <- 90L
id <- sprintf("S%03d",seq_len(n))
batch <- rep(c("B1","B2"),length.out=n)
age <- stats::runif(n,25,70)
a <- stats::rnorm(n)
b <- stats::rnorm(n)
c <- stats::rnorm(n)
intensity <- exp(0.55*a+0.012*(age-48)+0.30*(batch=="B2"))
event_time <- stats::rexp(n,rate=0.045*intensity)
censor_time <- stats::rexp(n,rate=0.012)
event <- as.integer(event_time<=censor_time)
time <- pmin(event_time,censor_time)
meta <- data.frame(sample_id=id,subject_id=id,batch=batch,age=age,
  survival_time=time,survival_event=event,stringsAsFactors=FALSE)
x <- cbind(A=a,B=b,C=c,CONSTANT=rep(1,n))
rownames(x) <- id
out_dir <- tempfile("multiomics_survival_reference_")
fit <- run_cox_survival_reference(x,meta,out_dir,covariates="age")
stopifnot(nrow(fit$results)==4L,
  fit$summary$attemptedFeatures==4L,
  fit$summary$estimable==3L,
  identical(fit$summary$tieMethod,"efron"),
  "CONSTANT"%in%fit$results$feature)
tab <- fit$results
ok <- is.finite(tab$pValue)
stopifnot(sum(ok)==3L,
  all(is.finite(tab$proportionalHazardsP[ok])),
  all(tab$hazardRatio[ok]>0),
  all.equal(tab$qValue[ok],
    stats::p.adjust(ifelse(is.finite(tab$pValue),tab$pValue,1),"BH")[ok],
    tolerance=1e-10)
)
eff <- tab$effect[tab$feature=="A"]
stopifnot(length(eff)==1L,is.finite(eff),eff>0)
# Cases deliberately refused before fitting any p-values.
duplicate <- meta
duplicate$subject_id[2] <- duplicate$subject_id[1]
err <- try(run_cox_survival_reference(x,duplicate,tempfile("survival_duplicate_"),covariates="age"),silent=TRUE)
stopifnot(inherits(err,"try-error"))
invalid <- meta
invalid$survival_event[2] <- 2
err <- try(run_cox_survival_reference(x,invalid,tempfile("survival_event_"),covariates="age"),silent=TRUE)
stopifnot(inherits(err,"try-error"))
few <- meta
few$survival_event <- as.integer(seq_len(n)<=4)
err <- try(run_cox_survival_reference(x,few,tempfile("survival_few_"),covariates="age"),silent=TRUE)
stopifnot(inherits(err,"try-error"))
missing_covariate <- meta
missing_covariate$age[3] <- NA_real_
err <- try(run_cox_survival_reference(x,missing_covariate,tempfile("survival_covariate_"),covariates="age"),silent=TRUE)
stopifnot(inherits(err,"try-error"))

# Validate the actual /reference-analysis router, not just the helper.
source("multiomics-engine/server.R")
meta$omic <- "proteomics"
meta$assay_id <- meta$sample_id
meta$sample_type <- "biological"
meta$condition <- "cohort"
meta$timepoint <- "T0"
m_csv <- paste(capture.output(utils::write.csv(meta,row.names=FALSE,na="")),collapse="\n")
mat <- data.frame(feature_id=colnames(x),t(x),check.names=FALSE)
x_csv <- paste(capture.output(utils::write.csv(mat,row.names=FALSE,na="")),collapse="\n")
payload <- list(
  metadataCsv=m_csv,
  matrices=list(proteomics=x_csv),
  dataTypes=list(proteomics="log_intensity"),
  columnMapping=list(
    sample_id="sample_id",subject_id="subject_id",assay_id="assay_id",
    omic="omic",sample_type="sample_type",condition="condition",
    timepoint="timepoint",batch="batch",survival_time="survival_time",
    survival_event="survival_event",age="age"
  ),
  protocol=list(objective="outcome",outcomeType="survival",
    designType="independent",longitudinal=FALSE,covariateColumns=list("age"))
)
reference <- run_backend_analysis(payload)
stopifnot(
  identical(reference$methods$proteomics_survival$status,"ok"),
  reference$methods$proteomics_survival$summary$estimable==3L,
  length(reference$methods$proteomics_survival$top)>0L,
  identical(reference$runtime$RVersion,as.character(getRversion())),
  identical(reference$runtime$RPlatform,R.version$platform),
  identical(reference$runtime$packageVersions$survival,
    as.character(utils::packageVersion("survival"))),
  !is.null(reference$runtime$limitation)
)

unlink(out_dir,recursive=TRUE)
cat("Reference survival backend: valid Cox fitted + full-feature-family BH; repeated subjects, bad events, few events, missing covariate refused | PASS\n")
