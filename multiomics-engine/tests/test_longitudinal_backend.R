#!/usr/bin/env Rscript
# Guarded 2-group longitudinal reference and reproducible negative controls.
# Differences of slopes are tested explicitly, NEVER an arbitrary coefficient.
source("multiomics-engine/advanced_methods.R")
if(!requireNamespace("lmerTest",quietly=TRUE))stop("lmerTest package needed")
set.seed(20261009)
subjects <- sprintf("S%03d",seq_len(28L))
times <- c(0,1,2)
ids <- rep(subjects,each=length(times))
t <- rep(times,length(subjects))
arm <- rep(ifelse(seq_along(subjects)<=14L,"control","treatment"),each=length(times))
random_intercepts <- rep(stats::rnorm(28,sd=1.0),each=length(times))
base <- 10 + 0.22*t + random_intercepts
data <- data.frame(
  subject_id=ids,condition=arm,timepoint=paste0("T",t),
  time=t,batch=rep(c("B1","B2"),length.out=length(t)),
  stringsAsFactors=FALSE
)
rownames(data) <- paste0(ids,"_T",t)
x <- cbind(
  SIGNAL=base + 0.9*t*(arm=="treatment")+stats::rnorm(length(t),sd=0.32),
  NULL_1=base+stats::rnorm(length(t),sd=0.4),
  NULL_2=base+stats::rnorm(length(t),sd=0.4),
  CONSTANT=rep(5,length(t)),
  MISSING=base+stats::rnorm(length(t),sd=0.5)
)
rownames(x) <- rownames(data)
x[1:4,"MISSING"] <- NA_real_
rownames(x) <- rownames(data)

# Direct helper and actual backend must preserve the entire feature family.
outdir <- tempfile("lmer_reference_")
fit <- run_lmer_matrix(
  x,data,outdir,fixed_formula="condition * time + batch",
  interaction_term="conditiontreatment:time"
)
stopifnot(
  nrow(fit)==ncol(x),
  setequal(fit$feature,colnames(x)),
  identical(fit$term[1],"conditiontreatment:time"),
  identical(fit$status[fit$feature=="CONSTANT"],"constant_feature"),
  all(is.na(fit$q_bh[!is.finite(fit$p_value)])),
  all.equal(
    fit$q_bh[is.finite(fit$p_value)],
    stats::p.adjust(ifelse(is.finite(fit$p_value),fit$p_value,1),"BH")[is.finite(fit$p_value)],
    tolerance=1e-10
  )
)
signal <- fit[fit$feature=="SIGNAL",]
stopifnot(nrow(signal)==1L,is.finite(signal$p_value),
          signal$estimate>0.65)
# Compare the chosen random-effect structure to an independent direct fit.
direct_formula <- if(signal$random_effect_structure=="random_intercept_and_slope")
  value ~ condition * time + batch + (1+time|subject_id) else
  value ~ condition * time + batch + (1|subject_id)
direct_data <- data
direct_data$value <- x[,"SIGNAL"]
direct <- suppressWarnings(lmerTest::lmer(direct_formula,data=direct_data,REML=FALSE,
  control=lme4::lmerControl(optimizer="bobyqa")))
ref <- summary(direct)$coefficients["conditiontreatment:time",]
stopifnot(
  abs(ref["Estimate"]-signal$estimate)<1e-5,
  abs(ref["Std. Error"]-signal$std_error)<1e-5,
  abs(ref["Pr(>|t|)"]-signal$p_value)<1e-5
)

bad <- data
bad$time[2] <- bad$time[1]
error <- try(validate_longitudinal_reference_design(bad),silent=TRUE)
stopifnot(inherits(error,"try-error"))
crossover <- data
crossover$condition[2] <- "treatment"
error <- try(validate_longitudinal_reference_design(crossover),silent=TRUE)
stopifnot(inherits(error,"try-error"))
ambiguous <- try(parse_longitudinal_time(c("Baseline","Week 1")),silent=TRUE)
stopifnot(inherits(ambiguous,"try-error"))
mixed <- try(parse_longitudinal_time(c("T0","D1")),silent=TRUE)
stopifnot(inherits(mixed,"try-error"))
wrong <- try(run_lmer_matrix(x,data,tempfile(),
  fixed_formula="condition * time + batch",interaction_term="batchB2"),silent=TRUE)
stopifnot(inherits(wrong,"try-error"))

# Actual endpoint router, to avoid passing with an unconnected helper.
source("multiomics-engine/server.R")
m <- data
m$sample_id <- rownames(data)
m$assay_id <- rownames(data)
m$omic <- "proteomics"
m$sample_type <- "biological"
m$technical_replicate <- "1"
m_csv <- paste(capture.output(utils::write.csv(m,row.names=FALSE,na="")),collapse="\n")
features <- data.frame(feature_id=colnames(x),t(x),check.names=FALSE)
x_csv <- paste(capture.output(utils::write.csv(features,row.names=FALSE,na="")),collapse="\n")
payload <- list(
  metadataCsv=m_csv,
  matrices=list(proteomics=x_csv),
  dataTypes=list(proteomics="log_intensity"),
  columnMapping=list(
    subject_id="subject_id",sample_id="sample_id",assay_id="assay_id",
    omic="omic",condition="condition",timepoint="timepoint",
    batch="batch",sample_type="sample_type",technical_replicate="technical_replicate"),
  protocol=list(objective="time",designType="repeated",longitudinal=TRUE,
    covariateColumns=list())
)
reference <- run_backend_analysis(payload)
method <- reference$methods$proteomics_longitudinal
stopifnot(
  identical(method$status,"ok"),
  method$summary$attemptedFeatures==5L,
  method$summary$estimableFeatures>=2L,
  method$summary$blockedFeatures>=1L,
  identical(method$summary$timeScale,"T"),
  length(method$top)==5L
)
broken <- payload
broken$metadataCsv <- sub('"T0"','"Baseline"',m_csv,fixed=TRUE)
stopifnot(!identical(broken$metadataCsv,m_csv))
blocked <- run_backend_analysis(broken)$methods$proteomics_longitudinal
stopifnot(identical(blocked$status,"blocked"))

# Null-label exchange across subjects (not across visits), declared seeds.
# A simple smoke check for GROSS calibration problems, not a formal FDR study.
rates <- numeric(8)
for(replication in seq_len(8L)) {
  set.seed(20261009 + replication)
  u <- rep(stats::rnorm(28,sd=0.8),each=3L)
  y <- 10+u+0.22*t+stats::rnorm(length(t),sd=0.45)
  nullx <- cbind(NEGATIVE_CONTROL=y,
                 CONTROL_2=y+stats::rnorm(length(t),sd=0.6))
  rownames(nullx) <- rownames(data)
  nullfit <- run_lmer_matrix(nullx,data,tempfile(),
    fixed_formula="condition * time + batch",
    interaction_term="conditiontreatment:time",
    random_slope=FALSE)
  p <- nullfit$p_value[is.finite(nullfit$p_value)]
  stopifnot(length(p)>=1L)
  rates[replication] <- mean(p<0.05)
}
stopifnot(mean(rates)<=0.30) # loose regression threshold, NOT size equivalence
cat("Longitudinal reference: explicit slope, matched lmerTest, full BH family, blocked bad plans, routed backend, 8 null simulations | PASS\n")
unlink(outdir,recursive=TRUE)
