# Whole-study p-value/Wald tail summaries under synthetic global null.
# Genes are NOT independent replicates; arm comparisons happen across
# independently simulated whole studies.
rnaseq_null_tail_metrics <- function(p,stat) {
  p<-as.numeric(p);stat<-as.numeric(stat)
  okp<-is.finite(p)&p>=0&p<=1
  okz<-is.finite(stat)
  if(sum(okp)<300L||sum(okz)<300L)stop("Too few finite native null statistics.")
  ap<-abs(stat[okz])
  data.frame(
    finite_p=sum(okp),finite_z=sum(okz),
    p_le_010=mean(p[okp]<=.10),
    p_le_005=mean(p[okp]<=.05),
    p_le_001=mean(p[okp]<=.01),
    p_le_0001=mean(p[okp]<=.001),
    absz_ge_1644854=mean(ap>=1.644854),
    absz_ge_1959964=mean(ap>=1.959964),
    absz_ge_2575829=mean(ap>=2.575829),
    absz_ge_3290527=mean(ap>=3.290527),
    wald_mean=mean(stat[okz]),
    wald_squared_mean=mean(stat[okz]^2),
    absz_q90=unname(stats::quantile(ap,.90,names=FALSE,type=8)),
    absz_q95=unname(stats::quantile(ap,.95,names=FALSE,type=8)),
    absz_q99=unname(stats::quantile(ap,.99,names=FALSE,type=8)),
    absz_q999=unname(stats::quantile(ap,.999,names=FALSE,type=8)),
    stringsAsFactors=FALSE)
}

rnaseq_null_tail_summary <- function(rows,replicates=32L,seed_offset=6000L) {
  scenarios<-c("balanced_base","balanced_depth_dispersion","unequal_partial_batch")
  methods<-c("native_default","oracle_true_dispersion")
  metrics<-c("p_le_010","p_le_005","p_le_001","p_le_0001",
    "absz_ge_1644854","absz_ge_1959964","absz_ge_2575829","absz_ge_3290527",
    "wald_squared_mean","absz_q95","absz_q99","absz_q999")
  needed<-c("scenario","method","replicate","seed","finite_p","finite_z",
    "bh_discoveries",metrics)
  if(!is.data.frame(rows)||!nrow(rows)||!all(needed%in%names(rows))||
     anyNA(rows[,needed])||
     !identical(sort(unique(as.character(rows$scenario))),sort(scenarios))||
     !identical(sort(unique(as.character(rows$method))),sort(methods))||
     any(rows$finite_p<300L|rows$finite_z<300L)||
     any(rows$bh_discoveries<0L)||
     any(!is.finite(as.matrix(rows[,metrics]))))
    stop("Invalid/incomplete whole-study null-tail evidence.")
  cells<-split(rows,list(rows$scenario,rows$method),drop=TRUE)
  if(length(cells)!=6L||any(vapply(cells,nrow,integer(1))!=replicates))
    stop("Missing prespecified scenario/method cells.")
  out<-list();paired<-list()
  for(s in scenarios){
    local<-rows[rows$scenario==s,,drop=FALSE]
    blocks<-lapply(methods,function(method){
      z<-local[local$method==method,,drop=FALSE]
      z<-z[order(z$replicate),,drop=FALSE]
      if(nrow(z)!=replicates||anyDuplicated(z$replicate)||
         !identical(as.integer(z$replicate),
           seq.int(as.integer(seed_offset+1L),as.integer(seed_offset+replicates))))
        stop("Missing/duplicated or post-selected replicate IDs.")
      z
    })
    if(!identical(blocks[[1]]$seed,blocks[[2]]$seed))
      stop("Paired methods were not fitted on identical generated cohorts.")
    for(j in seq_along(methods)){
      z<-blocks[[j]]
      for(metric in metrics){
        v<-z[[metric]]
        out[[length(out)+1L]]<-data.frame(
          scenario=s,method=methods[[j]],metric=metric,
          independent_studies=nrow(z),mean=mean(v),median=stats::median(v),
          q25=unname(stats::quantile(v,.25,type=8)),
          q75=unname(stats::quantile(v,.75,type=8)),
          biological_gene_independence_assumed=FALSE,
          universal_calibration_certified=FALSE,
          stringsAsFactors=FALSE)
      }
    }
    for(metric in metrics){
      d<-blocks[[1]][[metric]]-blocks[[2]][[metric]]
      paired[[length(paired)+1L]]<-data.frame(
        scenario=s,metric=metric,
        independent_paired_studies=length(d),
        mean_default_minus_true_dispersion=mean(d),
        median_default_minus_true_dispersion=stats::median(d),
        q25=unname(stats::quantile(d,.25,type=8)),
        q75=unname(stats::quantile(d,.75,type=8)),
        mechanism_causally_identified=FALSE,
        stringsAsFactors=FALSE)
    }
  }
  list(summary=do.call(rbind,out),paired=do.call(rbind,paired))
}
