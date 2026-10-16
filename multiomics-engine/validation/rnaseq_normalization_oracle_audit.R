#!/usr/bin/env Rscript
# MultiOmics P0: synthetic DESeq2 size-factor oracle diagnostic.
# No production inference changes. All generated features have true null group effect.
suppressPackageStartupMessages(library(DESeq2))
if (!requireNamespace("jsonlite",quietly=TRUE)) stop("jsonlite required")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L||!nzchar(args[[1]])) stop("Expected output directory")
out<-args[[1]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,"per_dataset.csv")))stop("Cannot overwrite evidence")
reps<-as.integer(Sys.getenv("PMX_RNASEQ_ORACLE_REPS","40"))
if(is.na(reps)||reps<40L||reps>500L)stop("40-500 independent simulations required")
alpha<-0.05
nfeatures<-400L
features<-sprintf("SIM_GENE_%04d",seq_len(nfeatures))
scenarios<-c(n40_balanced=40L,n80_balanced=80L)
arms<-c("estimated_size_factors","known_true_size_factors")

make_null<-function(n,seed){
 set.seed(seed)
 group<-factor(rep(c("control","treated"),each=n/2L),levels=c("control","treated"))
 batch<-factor(rep(c("B1","B2"),length.out=n),levels=c("B1","B2"))
 names<-sprintf("S%03d",seq_len(n))
 meta<-data.frame(group=group,batch=batch,row.names=names)
 if(qr(model.matrix(~batch+group,data=meta))$rank!=3L ||
    any(table(group,batch)!=n/4L))stop("Nonidentifiable/confounded design")
 base<-exp(runif(nfeatures,log(45),log(300)))
 dispersion<-runif(nfeatures,0.08,0.32)
 library<-exp(rnorm(n,0,.22))
 counts<-vapply(seq_len(n),function(i){
   mu<-base*library[[i]]*exp(log(1.35)*(as.character(batch[i])=="B2"))
   rnbinom(nfeatures,mu=mu,size=1/dispersion)
 },numeric(nfeatures))
 storage.mode(counts)<-"integer"
 rownames(counts)<-features
 colnames(counts)<-names
 list(counts=counts,meta=meta,library=library,dispersion=dispersion)
}
fit<-function(study,use_oracle){
 dds<-DESeq2::DESeqDataSetFromMatrix(countData=study$counts,
          colData=study$meta,design=~batch+group)
 keep<-rowSums(DESeq2::counts(dds)>=10)>=max(2L,ceiling(.2*ncol(dds)))
 if(!any(keep))stop("No eligible genes")
 dds<-dds[keep,]
 true_sf<-study$library/exp(mean(log(study$library)))
 if(use_oracle)sizeFactors(dds)<-true_sf
 dds<-DESeq2::DESeq(dds,quiet=TRUE)
 obtained<-as.numeric(sizeFactors(dds))
 if(any(!is.finite(obtained))||any(obtained<=0)||length(obtained)!=length(true_sf))
    stop("Invalid inferred size factors")
 if(use_oracle && !isTRUE(all.equal(obtained,true_sf,tolerance=1e-9)))
    stop("Known factors changed during fitting")
 res<-DESeq2::results(dds,
    contrast=c("group","treated","control"),alpha=alpha)
 p<-as.numeric(res$pvalue)
 q<-as.numeric(res$padj)
 if(any(is.finite(p)&(p<0|p>1))||any(is.finite(q)&(q<0|q>1)))
    stop("Invalid probabilities")
 f<-sum(is.finite(q)&q<=alpha)
 logcenter<-function(x)as.numeric(log(x)-mean(log(x)))
 d<-as.numeric(dispersions(dds));valid<-is.finite(d)&d>0
 if(!any(valid))stop("No estimated dispersions")
 truth<-study$dispersion[keep]
 list(n_subjects=ncol(dds),n_generated=nfeatures,n_fitted=nrow(dds),
      n_finite_p=sum(is.finite(p)),n_finite_q=sum(is.finite(q)),
      any_false=as.integer(f>0L),false_discoveries=f,
      mean_raw_p_le_05=mean(is.finite(p)&p<=alpha),
      mean_missing_q=mean(!is.finite(q)),
      log_sf_rmse=sqrt(mean((logcenter(obtained)-logcenter(true_sf))^2)),
      median_dispersion_ratio=median(d[valid]/truth[valid]))
}
rows<-list();errors<-list();k<-0L
for(si in seq_along(scenarios)){
 label<-names(scenarios)[[si]]
 n<-unname(scenarios[[si]])
 for(r in seq_len(reps)){
   seed<-as.integer(31415000L+si*100000L+r)
   study<-make_null(n,seed)
   for(arm in arms){
     tryCatch({
       metrics<-fit(study,identical(arm,"known_true_size_factors"))
       k<-k+1L
       rows[[k]]<-data.frame(scenario=label,replicate=r,seed=seed,
                             arm=arm,as.data.frame(metrics),
                             stringsAsFactors=FALSE)
     },error=function(e){
       errors[[length(errors)+1L]]<<-data.frame(
         scenario=label,replicate=r,seed=seed,arm=arm,
         message=conditionMessage(e))
     })
   }
 }
 if(length(rows))write.csv(do.call(rbind,rows),
               file.path(out,"partial_per_dataset.csv"),row.names=FALSE)
 if(length(errors))write.csv(do.call(rbind,errors),
               file.path(out,"errors.csv"),row.names=FALSE)
 cat("Completed",label,reps,"genuine simulated datasets\n")
}
if(!length(rows))stop("No valid DESeq2 fits")
d<-do.call(rbind,rows)
write.csv(d,file.path(out,"per_dataset.csv"),row.names=FALSE)
expected<-length(scenarios)*reps*length(arms)
if(nrow(d)!=expected||length(errors))stop("Incomplete evaluation: inspect archived errors")
groups<-split(d,list(d$scenario,d$arm),drop=TRUE)
summary<-do.call(rbind,lapply(groups,function(x){
 if(nrow(x)!=reps)stop("Incomplete scenario/arm")
 hits<-sum(x$any_false)
 ci<-binom.test(hits,nrow(x))$conf.int
 data.frame(scenario=x$scenario[[1]],arm=x$arm[[1]],
   n_datasets=nrow(x),n_subjects=x$n_subjects[[1]],
   null_datasets_with_false_calls=hits,
   empirical_null_FWER=hits/nrow(x),
   exact_95_lower=ci[1],exact_95_upper=ci[2],
   mean_raw_p_le_05=mean(x$mean_raw_p_le_05),
   median_log_sf_rmse=median(x$log_sf_rmse),
   median_estimated_true_dispersion_ratio=median(x$median_dispersion_ratio),
   FDR_certified=FALSE)
}))
rownames(summary)<-NULL
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
wide<-reshape(d[,c("scenario","replicate","seed","arm","any_false")],
      idvar=c("scenario","replicate","seed"),timevar="arm",direction="wide")
a<-paste0("any_false.",arms)
if(nrow(wide)!=reps*length(scenarios)||anyNA(wide[,a,drop=FALSE]))
 stop("Paired comparisons lack datasets")
paired<-do.call(rbind,lapply(split(wide,wide$scenario),function(x){
 delta<-x[[a[1]]]-x[[a[2]]]
 set.seed(41410000L+sum(utf8ToInt(x$scenario[[1]])))
 boots<-replicate(2000L,mean(sample(delta,length(delta),replace=TRUE)))
 ci<-quantile(boots,c(.025,.975))
 data.frame(scenario=x$scenario[[1]],n_matched_datasets=nrow(x),
   estimated_minus_oracle_FWER=mean(delta),
   paired_bootstrap_lower=as.numeric(ci[1]),
   paired_bootstrap_upper=as.numeric(ci[2]),
   interpretation="descriptive_only")
}))
write.csv(paired,file.path(out,"paired_differences.csv"),row.names=FALSE)
jsonlite::write_json(list(
  schema="multiomics-oracle-normalization-audit/v1",
  datasets_per_scenario=reps,scenarios=as.list(scenarios),
  n_genes=nfeatures,n_genuine_fits=nrow(d),expected_fits=expected,
  true_condition_effects=0,known_library_factors="geometric_centered_DGP_factors",
  group_model="~batch+group",alpha=alpha,
  DESeq2_version=as.character(packageVersion("DESeq2")),
  R_version=as.character(getRversion()),
  inferential_result="diagnostic_only",
  certification=FALSE,clinical_validity=FALSE,
  limitations=c(
   "Only synthetic NB nulls with balanced group/batch designs",
   "No automatic replacement of standard DESeq2",
   "True library size factors cannot be observed in real data",
   "Changing size factors also changes estimated dispersion",
   "Forty Monte Carlo replicates imply wide uncertainty")),
   file.path(out,"manifest.json"),auto_unbox=TRUE,pretty=TRUE)
print(summary,row.names=FALSE)
print(paired,row.names=FALSE)
cat("ORACLE NORMALIZATION DIAGNOSTIC COMPLETE; NO FDR CERTIFICATION\n")
