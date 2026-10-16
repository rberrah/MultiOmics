#!/usr/bin/env Rscript
# Independent public processed-peak-area quality-audit replay.
# Dataset: Hubert & Phillips, D. melanogaster, two separate experiments,
# source commit pinned in CI. No artificial group, batch or QC run order.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop("Usage: Rscript external_hubert_peak_area.R <Exp1 peak> <Exp1 meta> <Exp2 peak> <Exp2 meta> <out>")
source("multiomics-engine/ms_quantitative_assurance.R")
out<-args[[5L]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("External audit output exists; refuse to reuse old results")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
all_summary<-list()
for(exp in 1:2){
 peakpath<-args[[(exp-1)*2+1]]
 metapath<-args[[(exp-1)*2+2]]
 original<-utils::read.csv(peakpath,check.names=FALSE,stringsAsFactors=FALSE,
                           na.strings=c("NA","N/A",""))
 meta<-utils::read.csv(metapath,check.names=FALSE,stringsAsFactors=FALSE)
 stopifnot(all(c("Metabolite","Column","Ion Mode")%in%names(original)),
           !anyNA(original$Metabolite),!anyNA(original$Column),
           !anyNA(original[["Ion Mode"]]))
 first_measurement<-match("Ion Mode",names(original))+1L
 sample_cols<-names(original)[seq.int(first_measurement,ncol(original))]
 trimmed<-trimws(sample_cols)
 if(anyDuplicated(trimmed)||any(!nzchar(trimmed)))
   stop("Study column IDs become ambiguous after whitespace trimming")
 names(original)[match(sample_cols,names(original))]<-trimmed
 qc_names<-trimmed[grepl("^QC[ _]",trimmed,ignore.case=TRUE)]
 if(length(qc_names)!=if(exp==1L)6L else 3L)
   stop("Public study QC count has changed or external data are incomplete")
 if(exp==1L) {
   stopifnot(all(c("Sample","Strain","Population","PoolSize","Age")%in%names(meta)))
   bio_names<-as.character(meta$Sample)
   groups<-as.character(meta$Strain)
   if(nrow(meta)!=35L||!setequal(bio_names,setdiff(trimmed,qc_names)))
     stop("Experiment 1 independent biological sample IDs do not match matrix")
 }else {
   stopifnot(all(c("RunOrder","Tube","Sample","Diet","PoolSize")%in%names(meta)))
   bio_names<-as.character(meta$Tube)
   groups<-as.character(meta$Diet)
   if(nrow(meta)!=48L||!setequal(bio_names,setdiff(trimmed,qc_names)) ||
      anyNA(meta$RunOrder)||anyDuplicated(meta$RunOrder))
     stop("Experiment 2 independently documented biological sample IDs / orders disagree")
 }
 if(anyNA(bio_names)||anyDuplicated(bio_names)||any(!nzchar(groups)))
   stop("Missing/mispaired biological sample or outcome metadata")
 if(length(unique(groups))!=2L)
   stop("Unexpected external phenotype groups; never infer missing classes")
 if(!all(c(bio_names,qc_names)%in%names(original)))
   stop("Some published data columns have no source metadata")
 # Two publications' tables have QC peak areas but do not provide full
 # per-QC injection orders; the ordered integer below is solely a local key.
 # injection_order_verified=FALSE makes the engine suppress drift estimates.
 audit_order<-c(bio_names,qc_names)
 manifest<-data.frame(assay_id=audit_order,
   sample_type=c(rep("biological",length(bio_names)),rep("pooled_qc",length(qc_names))),
   condition=c(groups,rep(NA_character_,length(qc_names))),
   injection_order=seq_along(audit_order),
   injection_order_verified=rep(FALSE,length(audit_order)),
   stringsAsFactors=FALSE)
 combinations<-unique(original[,c("Column","Ion Mode")])
 if(nrow(combinations)<2L)stop("Expected multiple independent LC-MS acquisition panels")
 for(panel in seq_len(nrow(combinations))){
   col<-combinations$Column[[panel]]
   polarity<-combinations[["Ion Mode"]][[panel]]
   ix<-which(original$Column==col & original[["Ion Mode"]]==polarity)
   if(length(ix)<2L)stop("MS QC panel cannot be estimated with fewer than two features")
   raw<-as.matrix(original[ix,audit_order,drop=FALSE])
   suppressWarnings(storage.mode(raw)<-"double")
   if(any(is.infinite(raw))||any(raw<0,na.rm=TRUE) ||
      all(!is.finite(raw)))stop("Invalid source-native measured peak areas")
   ids<-paste0("E",exp,"_",seq_along(ix),"_",
               make.names(original$Metabolite[ix],unique=FALSE))
   if(anyDuplicated(ids))stop("Ambiguous external feature provenance IDs")
   values<-t(raw)
   rownames(values)<-audit_order
   colnames(values)<-ids
   qa<-ms_quantitative_assurance(values,manifest)
   if(qa$coverage$biological_assays!=length(bio_names) ||
      qa$coverage$pooled_qc_assays!=length(qc_names) ||
      qa$coverage$blank_assays!=0L ||
      isTRUE(qa$coverage$pooled_qc_injection_order_verified)||
      qa$coverage$pooled_qc_drift_screen_evaluable ||
      any(is.finite(qa$features$pooled_qc_order_spearman))||
      any(qa$features$flag_qc_drift_screen)||
      any(is.finite(qa$features$blank_over_biological_ratio))||
      qa$scientific_certified||
      !identical(qa$feature_exclusions_automatically_applied,FALSE))
      stop("External study missing controls were incorrectly claimed as observed")
   if(!any(is.finite(qa$features$pooled_qc_rsd)))
      stop("No numerical QC precision is evaluable on this independently measured panel")
   flags<-qa$features
   flags$source_metabolite<-original$Metabolite[ix]
   flags$source_panel<-as.character(col)
   flags$source_ion_mode<-as.character(polarity)
   label<-paste0("experiment_",exp,"_panel_",panel)
   utils::write.csv(flags,file.path(out,paste0(label,"_real_qc_features.csv")),
                    row.names=FALSE,na="")
   all_summary[[length(all_summary)+1L]]<-data.frame(
     experiment=exp,
     panel=as.character(col),ion_mode=as.character(polarity),
     real_features=length(ix),real_biological_samples=length(bio_names),
     real_pooled_qc_samples=length(qc_names),real_blank_samples=0L,
     qc_rsd_evaluable=sum(is.finite(qa$features$pooled_qc_rsd)),
     flags_qc_high_variation=sum(qa$features$flag_qc_variation),
     blank_contamination_evaluable=FALSE,
     qc_injection_order_documented=FALSE,drift_evaluable=FALSE,
     batch_identity_documented=FALSE,batch_effect_evaluable=FALSE,
     scientific_certified=FALSE,stringsAsFactors=FALSE)
 }
 utils::write.csv(manifest,
                  file.path(out,paste0("experiment_",exp,"_real_sample_roles.csv")),
                  row.names=FALSE,na="")
}
summary<-do.call(rbind,all_summary)
if(nrow(summary)<4L || !all(sort(unique(summary$experiment))==1:2))
  stop("Both independent public experiments and their assay panels required")
utils::write.csv(summary,file.path(out,"external_real_public_QC_summary.csv"),row.names=FALSE)
writeLines(c(
 "External public measured LC-MS peak-area tables: two independent D. melanogaster experiments, one originating research source.",
 "PUBLISHED DATA SOURCE (not reproduced): https://github.com/hubertdl/Metabolomics_sampling_project_D.melanogaster",
 "Pinned upstream source commit: 7c1649ab8093e2175e1683d4ee88af6ab1aa4fea",
 "DOI dataset: https://zenodo.org/records/21461823",
 "Matrix and metadata are downloaded on GitHub Actions; upstream exact Git blob hashes are checked before R execution.",
 "QC is real, but blanks and run order for QC are UNDOCUMENTED in these public tables, so drift and blank contamination are NOT evaluable.",
 "Experiment 2 reports biological RunOrder but not the QC sequence. No QC drift is inferred.",
 "Do not mislabel Population/Strain/Diet/PoolSize as instrumental analytical batch.",
 "No changes to observed peak areas, imputation, molecule identification, selection, biological effect estimate, FDR or automatic exclusion.",
 paste("Source checkout SHA:",Sys.getenv("GITHUB_SHA","local"))
),file.path(out,"SOURCE_AND_LIMITATIONS.txt"))
cat("PASS independent public measured peak-area QC:",nrow(summary),
    "panels from 2 actual experiments,",
    sum(summary$real_features),"feature-panel rows; no fabricated blanks/drift\n")
