# Frozen independent RNA-seq STRESS designs for REAL reference-model audit.
# Pure base R: easier to review/test separately from Bioconductor fitting.
#
# The normal mixed-truth suite uses strong log2FC=1.5; these stress designs
# deliberately reduce signal sizes or effective independent sample sizes.
rnaseq_stress_grid <- function() {
  data.frame(
    scenario=c("small_effect_20v20","small_cohort_6v6",
               "unequal_batch_12v28","sparse_signals_20v20"),
    n_control=c(20L,6L,12L,20L),
    n_treated=c(20L,6L,28L,20L),
    n_alternatives=c(20L,20L,20L,4L),
    abs_log2_fc=c(.35,.75,.75,.75),
    dispersion_low=c(.12,.22,.18,.12),
    dispersion_high=c(.28,.42,.36,.30),
    depth_low=c(.55,.70,.45,.55),
    depth_high=c(1.8,1.4,2.1,1.8),
    batch_mode=c("balanced","balanced","unequal","balanced"),
    stringsAsFactors=FALSE)
}
rnaseq_stress_generate <- function(scenario,replicate,n_genes=400L) {
  required <- c("scenario","n_control","n_treated","n_alternatives",
                "abs_log2_fc","dispersion_low","dispersion_high",
                "depth_low","depth_high","batch_mode")
  if (!is.data.frame(scenario) || nrow(scenario)!=1L ||
      !all(required %in% names(scenario)) || anyNA(scenario[,required]))
    stop("Exactly one complete, prespecified RNA-seq scenario is required.")
  if(length(replicate)!=1L || !is.finite(replicate) ||
     replicate<1L || replicate!=floor(replicate) || replicate>10000L ||
     length(n_genes)!=1L || !is.finite(n_genes) ||
     n_genes<100L || n_genes!=floor(n_genes))
    stop("Invalid independent replicate number or number of tested genes.")
  p <- scenario[1L,]
  if(!p$scenario %in% rnaseq_stress_grid()$scenario ||
     !identical(unname(as.character(p$batch_mode)),if (p$scenario=="unequal_batch_12v28") "unequal" else "balanced") ||
     any(!is.finite(as.numeric(p[,c("n_control","n_treated","n_alternatives",
                "abs_log2_fc","dispersion_low","dispersion_high",
                "depth_low","depth_high")]))) ||
     p$n_control<6L || p$n_treated<6L ||
     any(c(p$n_control,p$n_treated,p$n_alternatives)!=
         floor(c(p$n_control,p$n_treated,p$n_alternatives))) ||
     p$n_alternatives<2L || p$n_alternatives>=n_genes ||
     p$n_alternatives%%2L!=0L ||
     p$abs_log2_fc<=0 || p$dispersion_low<=0 ||
     p$dispersion_high<p$dispersion_low ||
     p$depth_low<=0 || p$depth_high<p$depth_low)
    stop("Malformed RNA-seq stress design.")
  # Explicit stable scenario ID avoids changing old RNG streams if grid order
  # changes. No data-derived reseeding or optional feature selection.
  s <- match(p$scenario,rnaseq_stress_grid()$scenario)
  seed <- as.integer(910000L+100000L*s+as.integer(replicate))
  set.seed(seed)
  n <- p$n_control+p$n_treated
  ids <- sprintf("S%03d",seq_len(n))
  groups <- factor(c(rep("control",p$n_control),rep("treated",p$n_treated)),
                   levels=c("control","treated"))
  if(p$batch_mode=="unequal") {
    # 9/3 controls vs 9/19 treated; batch adjustment is identifiable but
    # strongly correlated with treatment assignment.
    if(p$n_control!=12L || p$n_treated!=28L)
      stop("The partially confounded design requires the frozen 12/28 split.")
    batch <- c(sample(c(rep("B1",9L),rep("B2",3L))),
               sample(c(rep("B1",9L),rep("B2",19L))))
  } else {
    batch <- c(rep(c("B1","B2"),length.out=p$n_control),
               rep(c("B1","B2"),length.out=p$n_treated))
  }
  metadata <- data.frame(subject_id=sprintf("SUBJ%03d",seq_len(n)),
    condition=groups,batch=factor(batch,levels=c("B1","B2")),
    row.names=ids,stringsAsFactors=FALSE)
  matrix <- stats::model.matrix(~batch+condition,metadata)
  if(qr(matrix)$rank!=ncol(matrix) ||
     any(table(metadata$condition,metadata$batch)<2L) ||
     anyDuplicated(metadata$subject_id))
    stop("Unidentifiable or insufficiently replicated batch/condition plan.")
  features <- sprintf("G%04d",seq_len(n_genes))
  chosen <- sample.int(n_genes,size=p$n_alternatives,replace=FALSE)
  truth <- stats::setNames(numeric(n_genes),features)
  truth[chosen] <- rep(c(p$abs_log2_fc,-p$abs_log2_fc),
                       length.out=p$n_alternatives)
  dispersion <- seq(p$dispersion_low,p$dispersion_high,length.out=n_genes)
  baseline <- exp(seq(log(30),log(300),length.out=n_genes))
  depth <- exp(stats::runif(n,log(p$depth_low),log(p$depth_high)))
  # This nuisance covariate varies by gene but is encoded in the design.
  batch_multiplier <- exp(seq(-.35,.35,length.out=n_genes))
  counts <- vapply(seq_len(n),function(i) {
    mu <- baseline*depth[[i]]*
      if(metadata$batch[[i]]=="B2") batch_multiplier else 1
    if(groups[[i]]=="treated") mu <- mu*2^unname(truth)
    stats::rnbinom(n_genes,mu=mu,size=1/dispersion)
  },numeric(n_genes))
  rownames(counts) <- features
  colnames(counts) <- ids
  list(x=t(counts),meta=metadata,truth=truth,seed=seed,
       scenario=as.character(p$scenario),
       observed_depth_ratio=max(depth)/min(depth),
       design_condition_number=kappa(matrix),
       n_control=p$n_control,n_treated=p$n_treated,
       planted_abs_log2fc=p$abs_log2_fc)
}
