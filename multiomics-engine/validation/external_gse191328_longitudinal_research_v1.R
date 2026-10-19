#!/usr/bin/env Rscript
# Real public external execution qualification for MultiOmics native lmerTest path.
# Dataset: GEO GSE191328, Discovery-Infliximab RNA-seq.
# Prespecified design: week-14 remission group x baseline-to-week-2 time.
# Counts are transformed with TMM logCPM before lmerTest; raw counts are NOT
# passed into the confirmatory mixed-model route.
# This is numerical/model-contract qualification, NOT validation of remission
# biology, causal treatment effects, FDR, power or clinical use.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(edgeR)
  library(lmerTest)
  library(lme4)
  library(jsonlite)
  library(digest)
})
source("multiomics-engine/advanced_methods.R")

out <- Sys.getenv("MULTIOMICS_EXTERNAL_LONGITUDINAL_OUT",
                  "tmp/external-longitudinal-gse191328")
if (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Refusing to overwrite prior external longitudinal evidence.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw_dir <- file.path(out,"downloaded-counts")
dir.create(raw_dir,recursive=TRUE,showWarnings=FALSE)

gse_list <- GEOquery::getGEO("GSE191328",GSEMatrix=TRUE,getGPL=FALSE)
if (length(gse_list)!=1L) stop("Unexpected GSE191328 platform count.")
gse <- gse_list[[1L]]
pd <- Biobase::pData(gse)
if (nrow(pd)!=201L)
  stop("GSE191328 source sample count changed; review before qualification.")

normal <- function(x) gsub("[^a-z0-9]+","",tolower(as.character(x)))
# GEO pData also contains treatment_protocol_ch1. Prefer the semantic
# characteristic columns explicitly so a protocol paragraph can never be
# mistaken for the prespecified treatment label.
pick_characteristic <- function(label) {
  n <- names(pd)
  nn <- normal(n)
  target <- paste0(normal(label),"ch1")
  exact <- which(nn == target)
  if(length(exact)==1L)return(n[[exact]])
  # GEOquery versions may expose the parsed characteristic without :ch1.
  exact <- which(nn == normal(label))
  if(length(exact)==1L)return(n[[exact]])
  candidate <- which(grepl(normal(label),nn,fixed=TRUE) &
                     !grepl("protocol",nn,fixed=TRUE))
  if(length(candidate)!=1L)
    stop(paste("Ambiguous/missing GEO characteristic field:",label,
               "candidates=",paste(n[candidate],collapse="|")))
  n[[candidate]]
}
characteristic_value <- function(x,label) {
  raw <- trimws(as.character(x))
  # Some GEOquery representations retain "label: value"; others expose value.
  raw <- sub(paste0("^",label,"[[:space:]]*:[[:space:]]*"),"",
             raw,ignore.case=TRUE)
  trimws(raw)
}
patient_col <- pick_characteristic("patient")
cohort_col <- pick_characteristic("cohort")
treatment_col <- pick_characteristic("treatment")
remission_col <- pick_characteristic("remission at week 14")
time_col <- pick_characteristic("time point")

patient_raw <- characteristic_value(pd[[patient_col]],"patient")
cohort_raw <- characteristic_value(pd[[cohort_col]],"cohort")
treatment_raw <- characteristic_value(pd[[treatment_col]],"treatment")
remission_raw <- characteristic_value(pd[[remission_col]],"remission at week 14")
time_raw <- characteristic_value(pd[[time_col]],"time point")
cohort <- normal(cohort_raw)
treatment <- normal(treatment_raw)
tp <- normal(time_raw)
remission <- normal(remission_raw)
patient <- patient_raw
if(!all(c("infliximab") %in% unique(treatment)) ||
   !all(c("baseline","2w") %in% unique(tp)) ||
   !all(c("remission","nonremission") %in% unique(remission)) ||
   !any(grepl("discovery",cohort))) {
  compact <- function(x) paste(utils::head(sort(unique(as.character(x))),20L),collapse=" | ")
  stop(paste(
    "GSE191328 parsed characteristic levels do not match the prespecified public source contract.",
    "columns:",patient_col,cohort_col,treatment_col,remission_col,time_col,
    "cohort=",compact(cohort_raw),
    "treatment=",compact(treatment_raw),
    "remission=",compact(remission_raw),
    "time=",compact(time_raw)
  ))
}

keep <- grepl("discovery",cohort) & grepl("infliximab",cohort) &
        treatment=="infliximab" & tp %in% c("baseline","2w") &
        remission %in% c("remission","nonremission") &
        nzchar(patient)
sel <- pd[keep,,drop=FALSE]
if (!nrow(sel)) stop("Prespecified GSE191328 discovery-IFX baseline/2w subset absent.")

meta0 <- data.frame(
  sample_id=rownames(sel),
  subject_id=patient[keep],
  condition=ifelse(remission[keep]=="remission","remission","nonremission"),
  time=ifelse(tp[keep]=="baseline",0,2),
  time_label=time_raw[keep],
  cohort=cohort_raw[keep],
  treatment=treatment_raw[keep],
  stringsAsFactors=FALSE
)
if (anyDuplicated(paste(meta0$subject_id,meta0$time,sep="\r")))
  stop("Duplicate subject-time RNA-seq sample in prespecified source subset.")

# Eligibility is determined only from visit availability, before downloading
# expression or examining any feature-level result.
by_subject <- split(meta0,meta0$subject_id)
complete_ids <- names(by_subject)[vapply(by_subject,function(z)
  nrow(z)==2L && setequal(z$time,c(0,2)) && length(unique(z$condition))==1L,
  logical(1))]
meta <- meta0[meta0$subject_id %in% complete_ids,,drop=FALSE]
meta <- meta[order(meta$subject_id,meta$time),,drop=FALSE]
groups <- table(unique(meta[,c("subject_id","condition")])$condition)
if (length(unique(meta$subject_id)) < 10L || length(groups)!=2L || any(groups<5L))
  stop("GSE191328 baseline/week2 complete-pair design does not meet frozen lmer eligibility.")
# The clinical discovery cohort contains 14 patients (7/7), whereas the
# prespecified baseline+week-2 RNA-seq complete-pair rule yields 12 subjects,
# 6 remitters and 6 non-remitters in the current public GEO source. Exclusion
# is determined solely by visit availability before expression download or
# any feature-level fit; it is not selected from statistical results.
if (length(unique(meta$subject_id)) != 12L ||
    !identical(sort(as.integer(groups)),c(6L,6L)))
  stop(paste(
    "Unexpected GSE191328 complete-pair RNA-seq subset; review source.",
    "subjects=",length(unique(meta$subject_id)),
    "groups=",paste(names(groups),as.integer(groups),collapse="|")
  ))

# Resolve per-GSM supplementary count URL from GEO metadata. Prefer a URL in
# the series phenotype table; otherwise query the GSM record itself.
supp_cols <- grep("supplementary",names(pd),ignore.case=TRUE,value=TRUE)
resolve_url <- function(gsm) {
  row <- pd[gsm,,drop=FALSE]
  if(length(supp_cols)) {
    vals <- unlist(row[,supp_cols,drop=FALSE],use.names=FALSE)
    vals <- as.character(vals)
    vals <- vals[grepl("featureCounts.*txt[.]gz",vals,ignore.case=TRUE)]
    if(length(vals)) return(vals[[1L]])
  }
  obj <- GEOquery::getGEO(gsm,GSEMatrix=FALSE)
  m <- GEOquery::Meta(obj)
  vals <- unlist(m[grep("^supplementary_file",names(m))],use.names=FALSE)
  vals <- as.character(vals)
  vals <- vals[grepl("featureCounts.*txt[.]gz",vals,ignore.case=TRUE)]
  if(!length(vals)) stop(paste("No featureCounts supplement for",gsm))
  vals[[1L]]
}

read_counts <- function(gsm) {
  url <- resolve_url(gsm)
  url <- sub("^ftp://ftp[.]ncbi[.]nlm[.]nih[.]gov/","https://ftp.ncbi.nlm.nih.gov/",url)
  path <- file.path(raw_dir,paste0(gsm,".txt.gz"))
  utils::download.file(url,path,mode="wb",quiet=TRUE)
  if(!file.exists(path) || file.info(path)$size < 100000L)
    stop(paste("Incomplete GSE191328 count download:",gsm))
  tab <- utils::read.delim(gzfile(path),comment.char="#",check.names=FALSE,
                           stringsAsFactors=FALSE)
  if(ncol(tab)<2L || nrow(tab)<10000L)
    stop(paste("Unexpected featureCounts table:",gsm))
  ids <- as.character(tab[[1L]])
  count <- suppressWarnings(as.numeric(tab[[ncol(tab)]]))
  good <- nzchar(ids) & !grepl("^__",ids) & is.finite(count) & count>=0
  ids <- ids[good]; count <- count[good]
  if(length(ids)<10000L || anyDuplicated(ids))
    stop(paste("Invalid gene IDs/counts:",gsm))
  list(ids=ids,count=count,url=url,bytes=file.info(path)$size,
       sha256=digest::digest(file=path,algo="sha256"))
}

downloads <- lapply(meta$sample_id,read_counts)
ref_ids <- downloads[[1L]]$ids
if(any(vapply(downloads,function(z)!identical(z$ids,ref_ids),logical(1))))
  stop("GSE191328 featureCounts gene universe/order differs across selected samples.")
counts <- do.call(cbind,lapply(downloads,`[[`,"count"))
rownames(counts) <- ref_ids
colnames(counts) <- meta$sample_id

# Count-aware, outcome-independent preprocessing before the processed matrix
# enters the lmerTest contract. Filter and TMM normalization never inspect the
# remission interaction coefficient or p-values.
design_filter <- stats::model.matrix(~condition*time,data=meta)
y <- edgeR::DGEList(counts=counts)
filter <- edgeR::filterByExpr(y,design=design_filter)
if(sum(filter)<1000L) stop("Unexpectedly few expressed genes after edgeR filterByExpr.")
y <- y[filter,,keep.lib.sizes=FALSE]
y <- edgeR::calcNormFactors(y,method="TMM")
logcpm <- edgeR::cpm(y,log=TRUE,prior.count=2)

# Prespecified deterministic feature family, independent of model outcome.
candidate <- sort(rownames(logcpm)[vapply(seq_len(nrow(logcpm)),function(i)
  all(is.finite(logcpm[i,])) && stats::sd(logcpm[i,])>1e-10,logical(1))])
if(length(candidate)<64L) stop("Too few finite nonconstant GSE191328 genes.")
features <- candidate[seq_len(64L)]
matrix <- t(logcpm[features,meta$sample_id,drop=FALSE])
storage.mode(matrix) <- "double"
rownames(matrix) <- meta$sample_id
colnames(matrix) <- features
rownames(meta) <- meta$sample_id

wrapper <- run_lmer_matrix(
  matrix=matrix,metadata=meta,
  output_dir=file.path(out,"wrapper"),
  fixed_formula="condition * time",
  subject_column="subject_id",
  interaction_term=NULL,
  random_slope="auto"
)
if(nrow(wrapper)!=64L || !setequal(wrapper$feature,features))
  stop("Native longitudinal wrapper lost the prespecified feature family.")

# With exactly two visits per subject the frozen wrapper uses random intercept,
# not a random slope. Refit independently with the exact same statistical model.
groups_sorted <- sort(unique(meta$condition))
expected_term <- paste0("condition",groups_sorted[[2L]],":time")
meta$condition <- factor(meta$condition,levels=groups_sorted)
direct <- lapply(features,function(feature) {
  dat <- meta
  dat$value <- matrix[,feature]
  base_row <- data.frame(
    feature=feature,term=expected_term,
    estimate=NA_real_,std_error=NA_real_,df=NA_real_,
    statistic=NA_real_,p_value=NA_real_,q_bh=NA_real_,
    random_effect_structure=NA_character_,
    status="singular_or_nonconvergent_fit",
    stringsAsFactors=FALSE
  )
  fit <- try(suppressWarnings(lmerTest::lmer(
    value ~ condition*time + (1|subject_id),data=dat,REML=FALSE,
    control=lme4::lmerControl(optimizer="bobyqa")
  )),silent=TRUE)
  if(inherits(fit,"try-error") ||
     lme4::isSingular(fit,tol=1e-4) ||
     !is.null(fit@optinfo$conv$lme4$messages))
    return(base_row)
  co <- summary(fit)$coefficients
  if(!(expected_term %in% rownames(co)) ||
     !all(c("Estimate","Std. Error","df","t value","Pr(>|t|)") %in% colnames(co))) {
    base_row$status <- "missing_explicit_interaction"
    return(base_row)
  }
  z <- co[expected_term,,drop=FALSE]
  values <- c(z[1,"Estimate"],z[1,"Std. Error"],z[1,"df"],
              z[1,"t value"],z[1,"Pr(>|t|)"])
  if(any(!is.finite(values)) || values[[2L]]<=0 ||
     values[[5L]]<0 || values[[5L]]>1) {
    base_row$status <- "invalid_inference"
    return(base_row)
  }
  base_row$estimate <- unname(values[[1L]])
  base_row$std_error <- unname(values[[2L]])
  base_row$df <- unname(values[[3L]])
  base_row$statistic <- unname(values[[4L]])
  base_row$p_value <- unname(values[[5L]])
  base_row$random_effect_structure <- "random_intercept"
  base_row$status <- "ok"
  base_row
})
direct <- do.call(rbind,direct)
# Preserve the original 64-feature family in BH exactly as the wrapper does:
# failed fits contribute p=1 to the family denominator and retain q=NA.
direct_family_p <- ifelse(is.finite(direct$p_value),direct$p_value,1)
direct_all_q <- stats::p.adjust(direct_family_p,method="BH")
direct$q_bh[is.finite(direct$p_value)] <-
  direct_all_q[is.finite(direct$p_value)]
utils::write.csv(direct,file.path(out,"independent_direct_lmertest.csv"),row.names=FALSE)

w <- wrapper[match(features,wrapper$feature),,drop=FALSE]
d <- direct[match(features,direct$feature),,drop=FALSE]
if(anyNA(w$feature) || anyNA(d$feature))
  stop("External longitudinal feature matching failed.")
if(!identical(as.character(w$status),as.character(d$status)))
  stop("Wrapper/direct lmer fit-status decisions differ.")
estimable <- is.finite(w$p_value)
if(!any(estimable))
  stop("No estimable feature remained in the frozen external family.")
if(any(w$random_effect_structure[estimable]!="random_intercept"))
  stop("Two-visit external qualification unexpectedly used a random slope.")
if(any(!is.na(w$random_effect_structure[!estimable])))
  stop("Blocked longitudinal features must not advertise a fitted random-effects structure.")
fields <- c("estimate","std_error","df","statistic","p_value","q_bh")
max_abs <- setNames(numeric(length(fields)),fields)
for(name in fields) {
  a <- as.numeric(w[[name]][estimable])
  b <- as.numeric(d[[name]][estimable])
  delta <- abs(a-b)
  if(any(!is.finite(delta)))
    stop(paste("Non-finite external lmer comparison:",name))
  max_abs[[name]] <- max(delta)
  if(max_abs[[name]]>1e-8)
    stop(sprintf("GSE191328 lmer wrapper/direct mismatch %s: %.12g",name,max_abs[[name]]))
  if(any(is.na(w[[name]][estimable]) != is.na(d[[name]][estimable])))
    stop(paste("Wrapper/direct lmer missingness mismatch:",name))
}

feature_manifest <- data.frame(
  rank=seq_along(features),feature=features,
  selection_rule="lexicographic_first_64_after_edgeR_filterByExpr_TMM_logCPM_before_lmer_fit",
  stringsAsFactors=FALSE
)
utils::write.csv(feature_manifest,file.path(out,"prespecified_feature_family.csv"),row.names=FALSE)
utils::write.csv(meta,file.path(out,"longitudinal_metadata_used.csv"),row.names=FALSE)
utils::write.csv(data.frame(
  sample_id=meta$sample_id,
  source_url=vapply(downloads,`[[`,character(1),"url"),
  bytes=vapply(downloads,`[[`,numeric(1),"bytes"),
  sha256=vapply(downloads,`[[`,character(1),"sha256"),
  stringsAsFactors=FALSE
),file.path(out,"source_count_file_manifest.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

report <- list(
  schema="multiomics.external-longitudinal-gse191328-research-v1.v1",
  source=list(
    accession="GSE191328",
    pubmed=36153599L,
    doi="10.1186/s13073-022-01112-z",
    cohort="Discovery - Infliximab",
    endpointGrouping="remission at week 14",
    visits=c("baseline","2w")
  ),
  subjects=length(unique(meta$subject_id)),
  observations=nrow(meta),
  groups=as.list(table(unique(meta[,c("subject_id","condition")])$condition)),
  preprocessing="edgeR filterByExpr; TMM normalization; logCPM prior.count=2; then deterministic 64-feature family",
  prespecifiedSelection="first 64 lexicographically sorted filtered finite nonconstant Ensembl feature IDs; no lmer coefficient or p-value used",
  attemptedFeatures=length(features),
  estimableFeatures=sum(estimable),
  blockedSingularOrNonconvergent=sum(w$status=="singular_or_nonconvergent_fit"),
  blockedOther=sum(!estimable & w$status!="singular_or_nonconvergent_fit"),
  randomEffectStructure="random_intercept (two visits per subject by design); singular/nonconvergent features remain blocked in the complete family",
  expectedInteraction=expected_term,
  maxAbsoluteWrapperDirectDifference=as.list(max_abs),
  resultInterpretation="External real-data numerical execution qualification of the lmerTest contract. Remission grouping is observational; linear time with two visits is a slope/contrast parameterization, not proof of a biological trajectory or causal treatment effect.",
  softwareQualified=TRUE,
  statisticalCalibrationCertified=FALSE,
  biologicalValidation=FALSE,
  clinicalCertification=FALSE
)
jsonlite::write_json(report,file.path(out,"external_longitudinal_gse191328_evidence.json"),
                     pretty=TRUE,auto_unbox=TRUE,null="null")
cat(sprintf(
  "PASS GSE191328 external lmer qualification: %d subjects (%s), %d observations, %d prespecified genes, max numeric delta %.3g; NO FDR/biomarker/causal/clinical certification.\n",
  length(unique(meta$subject_id)),
  paste(names(groups),as.integer(groups),collapse=", "),
  nrow(meta),length(features),max(unlist(max_abs))
))
cat(sprintf(
  "Feature-fit status: %d estimable, %d singular/nonconvergent blocked, %d other blocked; complete 64-feature family retained.\n",
  sum(estimable),sum(w$status=="singular_or_nonconvergent_fit"),
  sum(!estimable & w$status!="singular_or_nonconvergent_fit")
))
