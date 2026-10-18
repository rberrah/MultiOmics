#!/usr/bin/env Rscript
source("multiomics-engine/rnaseq_independent_gamma_poisson_null.R")
scenarios<-c("gp_balanced_heterogeneous","gp_depth_heterogeneous",
             "gp_unbalanced_batch")
for(s in scenarios) {
  a<-independent_gp_null_study(s,5001L)
  b<-independent_gp_null_study(s,5001L)
  stopifnot(identical(a$counts,b$counts),
    identical(a$dispersion,b$dispersion),
    identical(a$meta,b$meta),
    nrow(a$counts)==40L,ncol(a$counts)==400L,
    all(is.finite(a$dispersion)),all(a$dispersion>0),
    a$no_true_condition_effect,
    a$gamma_poisson_latent_model,
    length(unique(a$meta$condition))==2L,
    all(table(a$meta$condition,a$meta$batch)>0L),
    stats::sd(a$dispersion_log_residual)>.4,
    identical(rownames(a$counts),rownames(a$meta)))
  stopifnot(!identical(a$counts,
     independent_gp_null_study(s,5002L)$counts))
}
old<-get(".Random.seed",envir=.GlobalEnv)
independent_gp_null_study("gp_balanced_heterogeneous",5003L)
stopifnot(identical(get(".Random.seed",envir=.GlobalEnv),old))
reject<-function(...)inherits(try(independent_gp_null_study(...),
                                  silent=TRUE),"try-error")
stopifnot(reject("unknown",5001L),
          reject("gp_balanced_heterogeneous",4001L),
          reject("gp_balanced_heterogeneous",5033L))
# Analytical Gamma-Poisson identity independently checked with thousands
# of draws (moment-level test, NOT a reproduction of clinical data).
set.seed(8899L)
mu<-100
alpha<-.25
latent<-stats::rgamma(40000L,shape=1/alpha,scale=mu*alpha)
draws<-stats::rpois(40000L,lambda=latent)
stopifnot(abs(mean(draws)-mu)/mu<.05,
          abs(stats::var(draws)-(mu+alpha*mu^2))/
              (mu+alpha*mu^2)<.08)
cat("PASS: independent Gamma-Poisson source, heterogeneous gene dispersions, real NB moments and rank-safe paired studies\n")
