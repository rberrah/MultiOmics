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

normal <- function(x) tolower(gsub("[^a-z0-9]+","",as.character(x)))
pick <- function(patterns) {
  n <- names(pd)
  nn <- normal(n)
  hit <- which(Reduce(`|`,lapply(patterns,function(p)grepl(p,nn,fixed=TRUE))))
  if(!length(hit)) stop(paste("Missing GEO metadata field:",paste(patterns,collapse="/")))
  n[[hit[[1L]]]]
}
patient_col <- pick(c("patient"))
cohort_col <- pick(c("cohort"))
treatment_col <- pick(c("treatment"))
remission_col <- pick(c("remissionatweek14"))
time_col <- pick(c("timepoint"))

cohort <- normal(pd[[cohort_col]])
treatment <- normal(pd[[treatment_col]])
tp <- normal(pd[[time_col]])
remission <- normal(pd[[remission_col]])
patient <- trimws(as.character(pd[[patient_col]]))

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
  time_label=as.character(pd[[time_col]][keep]),
  cohort=as.character(pd[[cohort_col]][keep]),
  treatment=as.character(pd[[treatment_col]][keep]),
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
# Published discovery cohort is 14 with 7 remission / 7 non-remission. Do not
# silently accept a materially changed public source.
if (length(unique(meta$subject_id)) != 14L || !identical(sort(as.integer(groups)),c(7L,7L)))
  stop("Unexpected GSE191328 discovery cohort size/remission balance; review source.")

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
  fit <- lmerTest::lmer(
    value ~ condition*time + (1|subject_id),data=dat,REML=FALSE,
    control=lme4::lmerControl(optimizer="bobyqa")
  )
  if(lme4::isSingular(fit,tol=1e-4) || !is.null(fit@optinfo$conv$lme4$messages))
    stop(paste("Independent lmer singular/nonconvergent:",feature))
  co <- summary(fit)$coefficients
  if(!(expected_term %in% rownames(co))) stop("Missing direct group-by-time term.")
  z <- co[expected_term,,drop=FALSE]
  data.frame(
    feature=feature,term=expected_term,
    estimate=unname(z[1,"Estimate"]),
    std_error=unname(z[1,"Std. Error"]),
    df=unname(z[1,"df"]),
    statistic=unname(z[1,"t value"]),
    p_value=unname(z[1,"Pr(>|t|)"]),
    random_effect_structure="random_intercept",
    stringsAsFactors=FALSE
  )
})
direct <- do.call(rbind,direct)
direct$q_bh <- stats::p.adjust(direct$p_value,method="BH")
utils::write.csv(direct,file.path(out,"independent_direct_lmertest.csv"),row.names=FALSE)

w <- wrapper[match(features,wrapper$feature),,drop=FALSE]
d <- direct[match(features,direct$feature),,drop=FALSE]
if(any(!w$status %in% c("ok","ok_random_slope_fallback")))
  stop("Unexpected failed lmer wrapper fit in frozen GSE191328 family.")
if(any(w$random_effect_structure!="random_intercept"))
  stop("Two-visit external qualification unexpectedly used a random slope.")
fields <- c("estimate","std_error","df","statistic","p_value","q_bh")
max_abs <- setNames(numeric(length(fields)),fields)
for(name in fields) {
  delta <- abs(as.numeric(w[[name]])-as.numeric(d[[name]]))
  if(any(!is.finite(delta)))stop(paste("Non-finite external lmer comparison:",name))
  max_abs[[name]] <- max(delta)
  if(max_abs[[name]]>1e-8)
    stop(sprintf("GSE191328 lmer wrapper/direct mismatch %s: %.12g",name,max_abs[[name]]))
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
  estimableFeatures=sum(is.finite(w$p_value)),
  randomEffectStructure="random_intercept (two visits per subject by design)",
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
