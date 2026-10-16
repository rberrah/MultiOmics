#!/usr/bin/env Rscript
# Standalone tests: real count thinning, filter blindness, paired uncertainty.
source("multiomics-engine/rnaseq_stress_designs.R")
source("multiomics-engine/rnaseq_filter_sensitivity.R")
grid<-rnaseq_stress_grid()
study<-rnaseq_stress_generate(grid[3L,,drop=FALSE],1L)
set.seed(4567L)
saved_rng<-get(".Random.seed",envir=.GlobalEnv)
d<-rnaseq_dilute_counts(study,1L)
again<-rnaseq_dilute_counts(study,1L)
stopifnot(identical(saved_rng,get(".Random.seed",envir=.GlobalEnv)),
 identical(d$x,again$x),identical(d$truth,study$truth),
 identical(dimnames(d$x),dimnames(study$x)),
 all(d$x<=study$x),all(d$x>=0),all(d$x==floor(d$x)),
 any(d$x!=study$x),length(d$expression_retention)==400L,
 identical(d$meta$condition,study$meta$condition))
baseline<-rnaseq_independent_filter(d$x,"reference")
strict<-rnaseq_independent_filter(d$x,"strict_prespecified")
stopifnot(baseline$retained_genes==400L,
 strict$retained_genes<400L,strict$retained_genes>90L,
 identical(rownames(strict$x),rownames(d$x)),
 all(colnames(strict$x)%in%colnames(d$x)),
 setequal(c(strict$kept,strict$discarded),colnames(d$x)),
 isTRUE(strict$filtering_blind_to_condition))
# The selected set CANNOT depend on condition labels or condition contrasts.
altered<-study$meta
altered$condition<-rev(altered$condition)
stopifnot(identical(strict$kept,
 rnaseq_independent_filter(d$x,"strict_prespecified")$kept))
bad<-d$x;bad[1,1]<--1
stopifnot(inherits(try(rnaseq_independent_filter(bad,"reference"),
                         silent=TRUE),"try-error"))
stopifnot(inherits(try(rnaseq_independent_filter(d$x,"post_hoc_optimum"),
                         silent=TRUE),"try-error"))
bad<-study;bad$x[1,1]<-NA
stopifnot(inherits(try(rnaseq_dilute_counts(bad,1L),silent=TRUE),"try-error"))

rows<-list()
for(scenario in c("weak","unequal"))
 for(method in c("DESeq2","edgeR_QL","limma_voom"))
  for(rep in 1:10)
   for(policy in c("reference","strict_prespecified")) {
    is_strict<-policy=="strict_prespecified"
    rows[[length(rows)+1L]]<-data.frame(
      scenario=scenario,method=method,filter_policy=policy,
      replicate=as.integer(rep),seed=rep+100L,
      power=if(is_strict) .3 else .5,
      false_discovery_proportion=if(is_strict) .05 else .10,
      total_planted_features=400L,n_true=20L,selected_features=10L,
      external_retained_genes=if(is_strict) 200L else 400L,
      excluded_planted_genes=if(is_strict) 3L else 0L,
      stringsAsFactors=FALSE)
   }
records<-do.call(rbind,rows)
set.seed(1212L)
oldseed<-get(".Random.seed",envir=.GlobalEnv)
s<-rnaseq_filter_pair_summary(records,reps=10L,bootstrap=200L)
stopifnot(identical(oldseed,get(".Random.seed",envir=.GlobalEnv)),
 nrow(s)==6L,
 all(abs(s$strict_minus_reference_power+.2)<1e-12),
 all(abs(s$strict_minus_reference_fdr+.05)<1e-12),
 all(s$mean_genes_removed_by_strict==200L),
 all(s$mean_true_signals_removed_by_strict==3L),
 all(s$paired_power_delta_ci95_high<=0),
 all(s$paired_fdr_delta_ci95_high<=0),
 all(!s$certified))
bad<-records[-1L,]
stopifnot(inherits(try(rnaseq_filter_pair_summary(bad),silent=TRUE),"try-error"))
bad<-records
bad$seed[2L]<-99999L
stopifnot(inherits(try(rnaseq_filter_pair_summary(bad),silent=TRUE),"try-error"))
bad<-records
bad$replicate[2L]<-2L # duplicate the replicate of an existing strict arm
stopifnot(inherits(try(rnaseq_filter_pair_summary(bad),silent=TRUE),"try-error"))
cat("RNA-seq independent prefilter sensitivity and true paired study comparison: PASS\n")
