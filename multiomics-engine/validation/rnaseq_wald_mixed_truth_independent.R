#!/usr/bin/env Rscript
# Independent follow-up to PR75; non-overlapping RNG seeds.
# This deliberately measures BOTH false discoveries and recovered signal
# for the two DESeq2 Wald references. Not a license to auto-select the winner.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Expected OUTPUT_DIRECTORY")
if(!requireNamespace("DESeq2",quietly=TRUE)||
   !requireNamespace("jsonlite",quietly=TRUE))
 stop("DESeq2 and jsonlite must actually be installed")
out<-args[[1L]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
nrep<-as.integer(Sys.getenv("PMX_RNASEQ_WALD_PAIRS","60"))
if(length(nrep)!=1L||is.na(nrep)||nrep<40L||nrep>500L)
 stop("Predeclared 40–500 independent Monte Carlo replicates required")
if(file.exists(file.path(out,"per_dataset.csv")))
 stop("Do not overwrite original results")
n<-40L;p<-400L
ids<-sprintf("BIO%03d",seq_len(n))
genes<-sprintf("ENSG_SIM%04d",seq_len(p))
truth<-c(rep(1.2,14L),rep(-1.2,14L),rep(0.65,12L),rep(0,p-40L))
cases<-list(
 null_balanced=list(nctrl=20L,mixed=FALSE),
 null_unbalanced=list(nctrl=12L,mixed=FALSE),
 mixed_balanced=list(nctrl=20L,mixed=TRUE),
 mixed_unbalanced=list(nctrl=12L,mixed=TRUE)
)
methods<-c("wald_normal","wald_t_reference")
rows<-list();errors<-list();z<-0L
for(ci in seq_along(cases)) {
 scenario<-names(cases)[[ci]]
 setting<-cases[[ci]]
 for(rep in seq_len(nrep)) {
   seed<-as.integer(30361013L+ci*100000L+rep)
   set.seed(seed)
   nctrl<-setting$nctrl
   group<-factor(c(rep("control",nctrl),rep("treated",n-nctrl)),
                 levels=c("control","treated"))
   batch<-if(nctrl==20L)
      factor(rep(c("B1","B2"),n/2L))
    else factor(c(rep("B1",9L),rep("B2",3L),
                  rep("B1",12L),rep("B2",16L)))
   md<-data.frame(group=group,batch=batch,row.names=ids)
   stopifnot(qr(model.matrix(~batch+group,data=md))$rank==3L)
   mu0<-exp(runif(p,log(45),log(300)))
   dispersion<-runif(p,0.08,0.32)
   library<-exp(rnorm(n,0,0.22))
   planted<-if(setting$mixed)truth else rep(0,p)
   X<-vapply(seq_len(n),function(i) {
     mu<-mu0*library[[i]]*exp(log(1.35)*(as.character(batch[i])=="B2"))*
           2^(planted*(as.character(group[i])=="treated"))
     rnbinom(p,mu=mu,size=1/dispersion)
   },numeric(p))
   rownames(X)<-genes;colnames(X)<-ids
   tryCatch({
     dds<-DESeq2::DESeqDataSetFromMatrix(countData=X,
                      colData=md,design=~batch+group)
     keep<-rowSums(DESeq2::counts(dds)>=10)>=max(2L,ceiling(ncol(dds)*0.2))
     if(!any(keep))stop("All genes filtered")
     dds<-dds[keep,]
     dds<-DESeq2::DESeq(dds,quiet=TRUE)
     control<-DESeq2::results(dds,contrast=c("group","treated","control"),alpha=.05)
     dds_t<-DESeq2::nbinomWaldTest(dds,useT=TRUE,quiet=TRUE)
     alternative<-DESeq2::results(dds_t,contrast=c("group","treated","control"),alpha=.05)
     for(method in methods){
       res<-if(method=="wald_normal")control else alternative
       ids_present<-rownames(res)
       idx<-match(ids_present,genes)
       if(anyNA(idx))stop("Unknown fitted features")
       selected<-is.finite(res$padj)&res$padj<=.05
       actual<-planted[idx]!=0
       R<-sum(selected);V<-sum(selected&!actual)
       T<-sum(selected&actual)
       if(V>R || T+V!=R)stop("Invalid truth accounting")
       z<-z+1L
       rows[[z]]<-data.frame(
        scenario=scenario,replicate=rep,seed=seed,method=method,
        n_subjects=n,n_features=p,n_planted=sum(planted!=0),
        n_tested=nrow(res),n_nonfinite_q=sum(!is.finite(res$padj)),
        discoveries=R,false_discoveries=V,true_discoveries=T,
        false_discovery_proportion=V/max(1L,R),
        power=if(setting$mixed)T/sum(planted!=0) else NA_real_,
        effect_mae_recovered=if(T>0L)
          mean(abs(as.numeric(res$log2FoldChange[selected&actual])-
                   planted[idx][selected&actual])) else NA_real_,
        any_false=as.integer(V>0L),stringsAsFactors=FALSE)
     }
   },error=function(e){
     errors[[length(errors)+1L]]<<-data.frame(scenario=scenario,
       replicate=rep,seed=seed,error=conditionMessage(e))
   })
 }
 if(length(rows))write.csv(do.call(rbind,rows),
       file.path(out,"partial_per_dataset.csv"),row.names=FALSE)
 if(length(errors))write.csv(do.call(rbind,errors),
       file.path(out,"errors.csv"),row.names=FALSE)
 cat("RNAseq WALD pair evaluation complete:",scenario,"(",nrep,")\n")
}
if(!length(rows))stop("No valid genuine fits")
d<-do.call(rbind,rows)
write.csv(d,file.path(out,"per_dataset.csv"),row.names=FALSE)
cells<-split(d,list(d$scenario,d$method),drop=TRUE)
bootstrap<-function(v,seed,nboot=2000L){
 set.seed(seed)
 means<-replicate(nboot,mean(sample(v,length(v),replace=TRUE)))
 unname(stats::quantile(means,c(.025,.975)))
}
summary<-do.call(rbind,lapply(cells,function(a){
  ncell<-nrow(a)
  if(ncell!=nrep)stop("Incomplete scenario-method cell")
  vals<-a$false_discovery_proportion
  b<-bootstrap(vals,30361013L+sum(utf8ToInt(a$scenario[1]))+
                sum(utf8ToInt(a$method[1])))
  fullnull<-all(a$n_planted==0)
  if(fullnull){
    exact<-stats::binom.test(sum(a$any_false),ncell)$conf.int
    b<-as.numeric(exact)
  }
  data.frame(scenario=a$scenario[[1]],method=a$method[[1]],
   n_independent_studies=ncell,mean_fdr=mean(vals),
   ci95_lower=b[1],ci95_upper=b[2],
   mean_power=if(fullnull)NA_real_ else mean(a$power),
   mean_true_discoveries=mean(a$true_discoveries),
   mean_false_discoveries=mean(a$false_discoveries),
   mean_effect_mae=if(fullnull)NA_real_ else mean(a$effect_mae_recovered,na.rm=TRUE),
   CI_method=if(fullnull)"Clopper-Pearson FWER=FDR under global null"
     else "2000 bootstrap dataset-level FDP draws",
   statistical_certification=FALSE)
}))
rownames(summary)<-NULL
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
manifest<-list(
 schema="pmx-rnaseq-wald-mixed-truth-independent/v1",
 seed_base=30361013L,independent_datasets=length(cases)*nrep,
 authentic_fits=nrow(d),expected_fits=length(cases)*nrep*2L,
 cases=names(cases),methods=methods,
 true_effects="40 of 400 planted log2FC values when mixed; all zero otherwise",
 r_version=as.character(getRversion()),
 deseq2_version=as.character(utils::packageVersion("DESeq2")),
 no_automatic_method_selection=TRUE,clinical_validation=FALSE,
 confirmatory_certification=FALSE,
 limitations=c(
  "Artificial negative-binomial generator; not an independent biological cohort",
  "Same cohort generators used for cases; external-model misspecification untested",
  "Repeated simulated datasets independent within cells, test arms paired within datasets",
  "Exploratory comparisons cannot establish general calibration or justify an automatic algorithm switch"))
jsonlite::write_json(manifest,file.path(out,"manifest.json"),pretty=TRUE,auto_unbox=TRUE)
print(summary,row.names=FALSE)
if(nrow(d)!=manifest$expected_fits||length(errors))
 stop("Incomplete mixed-truth study, retain ALL negative/error reports")
cat("MIXED TRUTH WALD FDR AND POWER COMPLETE: no automatic certification\n")
