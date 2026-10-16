#!/usr/bin/env Rscript
# Independent negative-binomial GLOBAL NULL diagnostic.
# NEVER change production DESeq2 based on this experiment.
suppressPackageStartupMessages(library(DESeq2))
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite required")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L||!nzchar(args[[1L]]))stop("Usage: script OUTPUT_DIR")
out<-args[[1L]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,"per_fit.csv")))stop("Refuse overwrite of calibration evidence")
nrep<-as.integer(Sys.getenv("PMX_RNASEQ_DISPERSION_ORACLE_REPS","40"))
if(length(nrep)!=1L||is.na(nrep)||nrep<40L||nrep>500L)stop("Use 40..500 independent replicates per scenario")
alpha<-.05
nfeatures<-400L
scenarios<-c(n40_balanced=40L,n80_balanced=80L)
arms<-c("estimated_sf_estimated_disp_wald_z",
        "true_sf_estimated_disp_wald_z",
        "true_sf_true_disp_wald_z",
        "true_sf_true_disp_wald_t")
sample_null<-function(n,seed){
 set.seed(seed)
 stopifnot(n%%4L==0L)
 group<-factor(rep(c("control","treated"),each=n/2L),levels=c("control","treated"))
 batch<-factor(rep(c("B1","B2"),length.out=n),levels=c("B1","B2"))
 smp<-sprintf("S%03d",seq_len(n))
 md<-data.frame(group=group,batch=batch,row.names=smp)
 stopifnot(all(table(group,batch)==n/4L),
           qr(stats::model.matrix(~batch+group,md))$rank==3L)
 gene_mu<-exp(runif(nfeatures,log(45),log(300)))
 dispersion<-runif(nfeatures,.08,.32)
 library<-exp(rnorm(n,0,.22))
 mat<-vapply(seq_len(n),function(j){
   mu<-gene_mu*library[[j]]*exp(log(1.35)*(as.character(batch[j])=="B2"))
   stats::rnbinom(nfeatures,mu=mu,size=1/dispersion)
 },numeric(nfeatures))
 storage.mode(mat)<-"integer"
 rownames(mat)<-sprintf("NB_%04d",seq_len(nfeatures))
 colnames(mat)<-smp
 list(counts=mat,meta=md,true_sf=library/exp(mean(log(library))),
      dispersion=dispersion)
}
run_fit<-function(study,arm){
 dds<-DESeq2::DESeqDataSetFromMatrix(
   countData=study$counts,colData=study$meta,design=~batch+group)
 keep<-rowSums(DESeq2::counts(dds)>=10)>=max(2L,ceiling(ncol(dds)*.2))
 if(!any(keep))stop("No eligible genes")
 dds<-dds[keep,]
 oracle_sf<-arm!="estimated_sf_estimated_disp_wald_z"
 oracle_disp<-grepl("_true_disp_",arm,fixed=TRUE)
 # grepl fixed means match the EXACT declared arm substring.
 if(oracle_sf)DESeq2::sizeFactors(dds)<-study$true_sf
 if(oracle_disp){
   DESeq2::dispersions(dds)<-study$dispersion[keep]
   if(!isTRUE(all.equal(as.numeric(DESeq2::dispersions(dds)),
                        study$dispersion[keep],tolerance=1e-12)))
      stop("True dispersion was not retained")
   dds<-DESeq2::nbinomWaldTest(dds,quiet=TRUE,
          useT=identical(arm,"true_sf_true_disp_wald_t"))
 }else{
   # Avoid outlier replacement which creates an uncontrolled extra comparison.
   dds<-DESeq2::DESeq(dds,quiet=TRUE,minReplicatesForReplace=Inf)
 }
 sf<-as.numeric(DESeq2::sizeFactors(dds))
 disp<-as.numeric(DESeq2::dispersions(dds))
 if(length(sf)!=length(study$true_sf)||any(!is.finite(sf))||
    any(!is.finite(disp))||any(disp<=0))stop("Invalid model parameters")
 if(oracle_sf&&!isTRUE(all.equal(sf,study$true_sf,tolerance=1e-9)))
    stop("Known size factors were not retained")
 truth<-study$dispersion[keep]
 if(oracle_disp&&!isTRUE(all.equal(disp,truth,tolerance=1e-10)))
    stop("Known dispersions were not retained")
 check<-function(independent_filtering){
   z<-DESeq2::results(dds,
      contrast=c("group","treated","control"),alpha=alpha,
      independentFiltering=independent_filtering)
   p<-as.numeric(z$pvalue);q<-as.numeric(z$padj)
   if(any(is.finite(p)&(p<0|p>1))||any(is.finite(q)&(q<0|q>1)))
     stop("Invalid p or q")
   c(n_finite_p=sum(is.finite(p)),n_finite_q=sum(is.finite(q)),
     raw_p_le_05=sum(is.finite(p)&p<=alpha),
     discoveries=sum(is.finite(q)&q<=alpha),
     any_false=as.integer(any(is.finite(q)&q<=alpha)),
     raw_p_lower_001=sum(is.finite(p)&p<.01))
 }
 m<-check(TRUE);no_filter<-check(FALSE)
 centered<-function(x)log(x)-mean(log(x))
 c(n_subjects=ncol(dds),n_features_fitted=nrow(dds),
   sf_log_rmse=sqrt(mean((centered(sf)-centered(study$true_sf))^2)),
   median_fitted_over_true_disp=median(disp/truth),
   m,stats::setNames(no_filter,paste0(names(no_filter),"_no_filter")))
}
result<-list();failures<-list();k<-0L
for(si in seq_along(scenarios)){
 label<-names(scenarios)[[si]]
 for(i in seq_len(nrep)){
   seed<-as.integer(41717000L+si*100000L+i)
   study<-sample_null(unname(scenarios[[si]]),seed)
   for(arm in arms)tryCatch({
     metrics<-run_fit(study,arm)
     k<-k+1L
     result[[k]]<-data.frame(scenario=label,replicate=i,seed=seed,
        arm=arm,as.list(metrics),check.names=FALSE,stringsAsFactors=FALSE)
   },error=function(e){
     failures[[length(failures)+1L]]<<-data.frame(
        scenario=label,replicate=i,seed=seed,arm=arm,
        error=conditionMessage(e),stringsAsFactors=FALSE)
   })
 }
 if(length(result))utils::write.csv(do.call(rbind,result),
    file.path(out,"partial_per_fit.csv"),row.names=FALSE)
 if(length(failures))utils::write.csv(do.call(rbind,failures),
    file.path(out,"errors.csv"),row.names=FALSE)
 cat("DISPERSION ORACLE SCENARIO COMPLETE:",label,nrep,"simulations\n")
}
if(!length(result))stop("No DESeq2 fits complete")
d<-do.call(rbind,result)
utils::write.csv(d,file.path(out,"per_fit.csv"),row.names=FALSE)
expected<-length(scenarios)*nrep*length(arms)
if(nrow(d)!=expected||length(failures))stop("Not all oracle comparison fits succeeded")
cols<-c("any_false","any_false_no_filter","raw_p_le_05","n_finite_p",
       "sf_log_rmse","median_fitted_over_true_disp")
for(x in cols)d[[x]]<-as.numeric(d[[x]])
cell<-split(d,list(d$scenario,d$arm),drop=TRUE)
summary<-do.call(rbind,lapply(cell,function(z){
 if(nrow(z)!=nrep)stop("Incomplete cell")
 x<-sum(z$any_false)
 y<-sum(z$any_false_no_filter)
 ci<-stats::binom.test(x,nrep)$conf.int
 ci_nf<-stats::binom.test(y,nrep)$conf.int
 data.frame(scenario=z$scenario[[1L]],arm=z$arm[[1L]],
   n_datasets=nrow(z),null_false_call_studies=x,
   global_null_fdr_fwer=x/nrep,exact95_lower=ci[1],exact95_upper=ci[2],
   null_fwer_without_independent_filter=y/nrep,
   nf_exact95_lower=ci_nf[1],nf_exact95_upper=ci_nf[2],
   mean_fraction_raw_p_le_05=mean(z$raw_p_le_05/z$n_finite_p),
   median_sf_log_rmse=stats::median(z$sf_log_rmse),
   median_dispersion_ratio=stats::median(z$median_fitted_over_true_disp),
   FDR_certified=FALSE,stringsAsFactors=FALSE)
}))
rownames(summary)<-NULL
utils::write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
w<-reshape(d[,c("scenario","replicate","seed","arm","any_false")],
    idvar=c("scenario","replicate","seed"),timevar="arm",direction="wide")
if(nrow(w)!=length(scenarios)*nrep)stop("Paired denominator mismatch")
pairs<-list();j<-0L
for(s in names(scenarios))for(i in seq_len(length(arms)-1L)){
 a<-w[w$scenario==s,,drop=FALSE]
 left<-as.numeric(a[[paste0("any_false.",arms[[i]])]])
 right<-as.numeric(a[[paste0("any_false.",arms[[i+1L]])]])
 if(anyNA(left)||anyNA(right))stop("Missing matched comparison")
 change<-left-right
 set.seed(51100000L+sum(utf8ToInt(s))+i)
 bootstrap<-replicate(2000L,mean(sample(change,length(change),replace=TRUE)))
 ci<-stats::quantile(bootstrap,c(.025,.975))
 j<-j+1L
 pairs[[j]]<-data.frame(scenario=s,left=arms[[i]],right=arms[[i+1L]],
   n_matched=nrow(a),mean_FWER_difference=mean(change),
   descriptive_bootstrap_lower=as.numeric(ci[1]),
   descriptive_bootstrap_upper=as.numeric(ci[2]))
}
utils::write.csv(do.call(rbind,pairs),file.path(out,"paired_differences.csv"),row.names=FALSE)
jsonlite::write_json(list(
 schema="multiomics-deseq2-true-dispersion-audit/v1",
 independent_simulations_per_scenario=nrep,designs=as.list(scenarios),
 arms=arms,n_genuine_fits=nrow(d),expected_fits=expected,
 all_genes_global_null=TRUE,n_features=nfeatures,
 model="~ batch + group",contrast=c("group","treated","control"),
 no_automatic_method_selection=TRUE,
 simulated_dispersion="Uniform(0.08,0.32); exactly known in oracle arms",
 oracle_size_factors="geometric_centered_DGP_factors",
 outlier_replacement_disabled_in_all_arms=TRUE,
 R=as.character(getRversion()),
 DESeq2=as.character(utils::packageVersion("DESeq2")),
 final_scientific_certification=FALSE,clinical_validation=FALSE,
 limitations=c(
  "Only balanced negative-binomial global-null generator",
  "True size factors and dispersions cannot be measured directly in patient samples",
  "Only 40 replicates per arm, confidence intervals are wide",
  "Changing dispersions also changes fitted GLM variances and coefficients",
  "This is diagnosis, not independent clinical or biological validation")),
 file.path(out,"manifest.json"),pretty=TRUE,auto_unbox=TRUE)
print(summary,row.names=FALSE)
print(do.call(rbind,pairs),row.names=FALSE)
cat("ORACLE DISPERSION DIAGNOSTIC COMPLETE; no FDR certification\n")
