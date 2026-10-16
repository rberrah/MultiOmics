#!/usr/bin/env Rscript
# Independent survival::coxph R reference: unstratified Cox, Breslow ties.
# cox.zph is a diagnostic, not a certificate of proportional hazards.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript independent_cox_reference.R observations.csv browser-cox.csv r-cox.csv")
if(!requireNamespace("survival",quietly=TRUE))stop("R survival package not installed")
data <- utils::read.csv(args[1],stringsAsFactors=FALSE)
browser <- utils::read.csv(args[2],stringsAsFactors=FALSE)
if(nrow(browser)!=80L||anyDuplicated(browser$case_id))stop("Expected exactly 80 unique Cox cases")
fits <- vector("list",nrow(browser))
for(i in seq_len(nrow(browser))) {
  identifier <- browser$case_id[i]
  dat <- data[data$case_id==identifier,,drop=FALSE]
  if(nrow(dat)!=browser$n[i]||any(dat$time<=0)||any(!dat$event%in%c(0L,1L)))
    stop("Bad survival input ",identifier)
  model <- survival::coxph(
    survival::Surv(time,event)~feature+batch,
    data=dat,ties="breslow",x=TRUE,
    control=survival::coxph.control(iter.max=80,timefix=FALSE)
  )
  coefficient <- unname(stats::coef(model)["feature"])
  se <- sqrt(stats::vcov(model)["feature","feature"])
  if(!is.finite(coefficient)||!is.finite(se)||!(se>0))
    stop("Undefined independent Cox estimate ",identifier)
  p <- 2*stats::pnorm(-abs(coefficient/se))
  zph <- try(survival::cox.zph(model,transform="km"),silent=TRUE)
  ph_p <- if(inherits(zph,"try-error")) NA_real_ else unname(zph$table["feature","p"])
  fits[[i]] <- data.frame(
    case_id=identifier,
    beta=coefficient,se=se,p=p,
    ph_p=ph_p,converged=TRUE,
    stringsAsFactors=FALSE
  )
}
output <- do.call(rbind,fits)
utils::write.csv(output,args[3],row.names=FALSE,quote=FALSE,na="")
cat(sprintf("Independent R survival::coxph Breslow: %d fits; feature cox.zph diagnostic recorded | PASS\n",nrow(output)))
