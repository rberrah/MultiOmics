#!/usr/bin/env Rscript
source("multiomics-engine/rnaseq_deseq2_null_tail_metrics.R")
p<-c(seq(.001,.999,length.out=398),.0001,.9999)
z<-stats::qnorm(p/2,lower.tail=FALSE)
z[seq_along(z)%%2L==0L]<- -z[seq_along(z)%%2L==0L]
x<-rnaseq_null_tail_metrics(p,z)
stopifnot(nrow(x)==1L,x$finite_p==400L,x$finite_z==400L,
          x$p_le_005>0,x$absz_q99>0)
scenarios<-c("balanced_base","balanced_depth_dispersion","unequal_partial_batch")
rows<-do.call(rbind,lapply(seq_along(scenarios),function(k) {
  do.call(rbind,lapply(6001:6032,function(i) {
    seed<-780000L+k*10000L+i
    data.frame(scenario=scenarios[[k]],
      method=c("native_default","oracle_true_dispersion"),
      replicate=i,seed=seed,finite_p=400L,finite_z=400L,bh_discoveries=0L,
      p_le_010=c(.11,.10),p_le_005=c(.06,.05),p_le_001=c(.015,.01),
      p_le_0001=c(.002,.001),absz_ge_1644854=c(.11,.10),
      absz_ge_1959964=c(.06,.05),absz_ge_2575829=c(.015,.01),
      absz_ge_3290527=c(.002,.001),wald_squared_mean=c(1.08,1.01),
      absz_q95=c(2.05,1.96),absz_q99=c(2.70,2.58),absz_q999=c(3.5,3.3))
  }))
}))
ans<-rnaseq_null_tail_summary(rows)
stopifnot(nrow(ans$summary)==72L,nrow(ans$paired)==36L,
  all(ans$paired$independent_paired_studies==32L),
  all(!ans$summary$biological_gene_independence_assumed),
  all(!ans$summary$universal_calibration_certified))
reject<-function(z)inherits(try(rnaseq_null_tail_summary(z),silent=TRUE),"try-error")
stopifnot(reject(rows[-1,,drop=FALSE]))
bad<-rows;bad$seed[[1]]<-1L;stopifnot(reject(bad))
bad<-rows;bad$p_le_005[[1]]<-NA_real_;stopifnot(reject(bad))
cat("PASS: whole-study null p/Wald-tail metrics; genes never treated as independent calibration replicates\n")
