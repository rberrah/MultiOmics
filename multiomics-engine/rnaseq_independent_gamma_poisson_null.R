# Independent Gamma-Poisson generator: no sourcing/calling the original
# rnbinom() study builder. Genuine independent gene-specific random NB
# dispersion residuals, NOT forced to lie on a smooth dispersion-mean curve.
# Negative-binomial variance: Var(Y|mu,alpha)=mu+alpha*mu^2.
independent_gp_null_study <- function(scenario,replicate) {
  scenarios<-c("gp_balanced_heterogeneous","gp_depth_heterogeneous",
               "gp_unbalanced_batch")
  k<-match(scenario,scenarios)
  if(length(scenario)!=1L||is.na(k)||
     length(replicate)!=1L||!is.finite(replicate)||
     replicate!=floor(replicate)||replicate<5001L||replicate>5032L)
    stop("Only prespecified novel Gamma-Poisson study scenarios and seed block 5001–5032 allowed.")
  seed<-as.integer(910000L+10000L*k+as.integer(replicate))
  existed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(existed)previous<-get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(existed)assign(".Random.seed",previous,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(seed)
  g<-400L
  n_control<-if(k==3L)12L else 20L
  n_treated<-40L-n_control
  samples<-sprintf("GP_S%02d",seq_len(40L))
  group<-factor(c(rep("control",n_control),rep("treated",n_treated)),
                levels=c("control","treated"))
  if(k==3L) {
    # Both condition/batch cells remain nonempty (no rank deficiency).
    batch<-c(sample(c(rep("B1",10L),rep("B2",2L))),
             sample(c(rep("B1",7L),rep("B2",21L))))
  } else {
    batch<-c(rep(c("B1","B2"),n_control/2L),
             rep(c("B1","B2"),n_treated/2L))
  }
  meta<-data.frame(condition=group,
                   batch=factor(batch,levels=c("B1","B2")),
                   row.names=samples)
  if(qr(stats::model.matrix(~batch+condition,meta))$rank!=3L ||
     any(table(meta$condition,meta$batch)==0L))
    stop("Independent alternative source has invalid rank-deficient batch design.")
  means<-exp(stats::runif(g,log(30),log(500)))
  # Unlike first generator, the true alpha has substantial independent
  # gene-by-gene variability even at a fixed mean (log-sd 0.65).
  smooth_alpha<-0.04+9/means
  residual<-stats::rnorm(g,mean=0,sd=.65)
  alpha<-smooth_alpha*exp(residual)
  if(k==2L)alpha<-alpha*2.5
  if(any(!is.finite(alpha)|alpha<=0)||
     stats::sd(residual)<.4)
    stop("Independent gene-level alpha heterogeneity was not generated.")
  depth<-if(k==2L)exp(stats::runif(40L,log(.25),log(3.5)))
         else stats::runif(40L,.85,1.15)
  gene_batch_log_effect<-stats::rnorm(g,mean=0,sd=.32)
  b2<-as.numeric(meta$batch=="B2")
  # No effect of treatment group in this conditional mean model.
  mu<-outer(means,depth)*exp(outer(gene_batch_log_effect,b2))
  shape<-matrix(rep(1/alpha,40L),nrow=g,ncol=40L)
  scale<-mu*matrix(rep(alpha,40L),nrow=g,ncol=40L)
  if(any(!is.finite(mu)|mu<=0)||
     any(!is.finite(shape)|shape<=0)||
     any(!is.finite(scale)|scale<=0))
    stop("Gamma–Poisson generator has invalid latent parameters.")
  latent<-stats::rgamma(g*40L,shape=as.vector(shape),scale=as.vector(scale))
  sampled<-stats::rpois(g*40L,lambda=latent)
  counts<-matrix(sampled,nrow=g,ncol=40L,
                 dimnames=list(sprintf("GP_N%04d",seq_len(g)),samples))
  if(anyNA(counts)||any(counts<0)||
     any(counts!=floor(counts))||
     !identical(dim(counts),c(400L,40L)))
    stop("Corrupted independently simulated whole-study RNA-seq counts.")
  list(
    counts=t(counts),meta=meta,
    dispersion=alpha,depth=depth,
    gene_means=means,
    dispersion_log_residual=residual,
    gamma_poisson_latent_model=TRUE,
    no_true_condition_effect=TRUE,
    batch_log_effect=gene_batch_log_effect,
    scenario=scenario,
    replicate=as.integer(replicate),seed=seed)
}
