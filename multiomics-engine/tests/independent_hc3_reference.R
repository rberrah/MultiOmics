#!/usr/bin/env Rscript
# Independent numerical reference for browser feature ~ condition + batch + age.
# Uses only R base lm/QR, with HC3 sandwich constructed directly.
# This is NOT a statistical guarantee under MNAR and does not validate the
# upstream feature-selection or omics preprocessing from the JS pipeline.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript independent_hc3_reference.R observations.csv browser-estimates.csv r-estimates.csv")
data <- utils::read.csv(args[1],stringsAsFactors=FALSE,check.names=FALSE,na.strings=c("","NA"))
browser <- utils::read.csv(args[2],stringsAsFactors=FALSE,check.names=FALSE)
if(anyDuplicated(paste(browser$case_id,browser$feature)))stop("Duplicate feature in browser reference contract")
data$condition <- factor(data$condition,levels=c("control","treated"))
data$batch <- factor(data$batch,levels=c("B1","B2"))
data$age <- as.numeric(data$age)
data$value <- as.numeric(data$value)
if(any(!is.finite(data$age)))stop("Invalid age covariate in simulation")
source_keys <- paste(data$case_id,data$feature,sep="|")
observations_by_key <- split(seq_len(nrow(data)),source_keys)
ref <- vector("list",nrow(browser))
for(i in seq_len(nrow(browser))) {
  key <- paste(browser$case_id[i],browser$feature[i],sep="|")
  idx <- observations_by_key[[key]]
  if(is.null(idx))stop("Missing independent source observations for ",key)
  dat <- data[idx,,drop=FALSE]
  dat <- dat[is.finite(dat$value),,drop=FALSE]
  if(nrow(dat)<8L || length(unique(dat$condition))!=2L)stop("Nonestimable reference ",key)
  fit <- stats::lm(value ~ condition+batch+age,data=dat)
  X <- stats::model.matrix(fit)
  e <- stats::residuals(fit)
  coef <- stats::coef(fit)
  if(any(!is.finite(coef)) || qr(X)$rank != ncol(X))stop("Aliased R design ",key)
  ci <- which(colnames(X)=="conditiontreated")
  if(length(ci)!=1L)stop("Expected exactly one treatment coefficient in ",key)
  # Reference variance HC3 = (X'X)^-1 X' diag((e/(1-h))²) X (X'X)^-1
  bread <- solve(crossprod(X))
  h <- stats::hatvalues(fit)
  if(any(!is.finite(h))||any(h>=1-1e-7))stop("HC3 leverage is unidentified: ",key)
  meat <- crossprod(X * as.numeric(e/(1-h)))
  vcov_hc3 <- bread %*% meat %*% bread
  se <- sqrt(vcov_hc3[ci,ci])
  beta <- as.numeric(coef[ci])
  df <- stats::df.residual(fit)
  if(!is.finite(se)||se<=0||df<1L)stop("Invalid reference uncertainty ",key)
  p <- 2*stats::pt(-abs(beta/se),df=df)
  t95 <- stats::qt(0.975,df)
  ref[[i]] <- data.frame(
    case_id=browser$case_id[i],scenario=browser$scenario[i],
    n=browser$n[i],replicate=browser$replicate[i],
    effect_planted=browser$effect_planted[i],
    feature=browser$feature[i],effect=beta,
    standardError=se,pValue=p,qValue=NA_real_,
    ciLow=beta-t95*se,ciHigh=beta+t95*se,
    stringsAsFactors=FALSE
  )
}
r <- do.call(rbind,ref)
for(c in unique(r$case_id)) {
  idx <- which(r$case_id==c)
  r$qValue[idx] <- stats::p.adjust(r$pValue[idx],method="BH")
}
if(any(!is.finite(as.matrix(r[,c("effect","standardError","pValue","qValue","ciLow","ciHigh")]))))
  stop("Non-finite independent R reference statistics")
utils::write.csv(r,file=args[3],row.names=FALSE,na="",quote=FALSE)
cat(sprintf("Independent R lm/HC3: %d matched feature-level models, %d datasets | PASS\n",
  nrow(r),length(unique(r$case_id))))
