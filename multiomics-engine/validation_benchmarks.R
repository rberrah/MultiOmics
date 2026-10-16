# Scientific benchmark helpers; deterministic and independent of MOFA2/mixOmics
# until a reference package is explicitly invoked. Not clinical certification.

simulate_mofa_known_truth <- function(seed=20261009L,
                                      scenario=c("biological","biological_plus_batch","batch_only"),
                                      n_subjects=90L) {
  scenario<-match.arg(scenario)
  if(n_subjects<60L || n_subjects %% 2L != 0L)stop("Synthetic MOFA needs even n>=60.")
  set.seed(as.integer(seed))
  ids<-sprintf("SIM%03d",seq_len(n_subjects))
  biological<-as.numeric(scale(stats::rnorm(n_subjects)))
  batch<-sample(rep(c("B1","B2"),each=n_subjects/2L))
  technical<-as.numeric(scale(as.integer(batch=="B2")))
  z<-if(scenario=="batch_only")rep(0,n_subjects) else biological
  t<-if(scenario=="biological")rep(0,n_subjects) else technical
  truth <- data.frame(sample=ids,biological=biological,
                      technical=technical,batch=batch)
  gen <- function(view,n_features) {
    v <- matrix(stats::rnorm(n_subjects*n_features,sd=0.55),
       nrow=n_subjects,ncol=n_features,
       dimnames=list(ids,paste0(view,"_",seq_len(n_features))))
    v[,seq_len(12L)]<-v[,seq_len(12L)]+
       outer(z,rep(2.0,12L))
    v[,13:24]<-v[,13:24]+outer(t,rep(2.5,12L))
    v
  }
  list(blocks=list(transcriptomics=gen("rna",32L),
                   proteomics=gen("protein",30L)),
       truth=truth,scenario=scenario,seed=as.integer(seed),
       expectedBiological=(scenario!="batch_only"),
       expectedTechnical=(scenario!="biological"))
}

mofa_factor_truth_diagnostics <- function(factors,truth,
                                          biological_threshold=0.65,
                                          batch_threshold=0.65) {
  # Evaluation uses truth ONLY AFTER fitting, not as a feature or covariate.
  if(!is.data.frame(factors) ||
     !all(c("sample","factor","value")%in%names(factors)) ||
     !all(c("sample","biological","technical","batch")%in%names(truth)))
    stop("MOFA factor score and simulation truth schemas differ.")
  if(anyNA(factors[,c("sample","factor","value")]) ||
     anyDuplicated(paste(factors$sample,factors$factor,sep="\r")))
    stop("Duplicated or missing factor/samples; do not silently average.")
  ids<-as.character(truth$sample)
  f_names<-unique(as.character(factors$factor))
  if(length(f_names)<1L || anyNA(f_names) ||
     any(!nzchar(f_names)) || anyDuplicated(ids))
    stop("Missing MOFA factors or duplicate truth IDs.")
  if(!setequal(unique(as.character(factors$sample)),ids))
    stop("MOFA score subjects do not match simulated truth; no silent deletion.")
  matrix_values<-vapply(f_names,function(name) {
    rows<-factors[as.character(factors$factor)==name,,drop=FALSE]
    if(nrow(rows)!=length(ids))stop("Factor lacks subjects:",name)
    v<-as.numeric(rows$value[match(ids,as.character(rows$sample))])
    if(any(!is.finite(v)))stop("Nonfinite MOFA scores:",name)
    v
  },numeric(length(ids)))
  if(is.null(dim(matrix_values)))matrix_values<-matrix(matrix_values,ncol=1L)
  colnames(matrix_values)<-f_names
  strength<-function(vector) {
    vapply(seq_len(ncol(matrix_values)),function(j) {
      x<-matrix_values[,j]
      if(stats::sd(x)<1e-10)return(0)
      abs(stats::cor(x,vector))
    },numeric(1))
  }
  biological<-strength(truth$biological)
  technical<-strength(truth$technical)
  if(any(!is.finite(biological))||any(!is.finite(technical)))
    stop("Unidentifiable truth/factor correlation.")
  # Sign and factor ordering are arbitrary in latent factor models.
  best_biological<-which.max(biological)
  best_technical<-which.max(technical)
  # Biological and technical truth must match DIFFERENT recovered factors.
  # Two high correlations against the same mixed factor do not establish
  # successful recovery of two independent latent processes.
  distinct_best <- list(score=-Inf,biological=NA_integer_,technical=NA_integer_)
  if(length(biological)>=2L)for(i in seq_along(biological))
    for(j in seq_along(technical))if(i!=j) {
      score <- min(biological[[i]],technical[[j]])
      if(score>distinct_best$score)
        distinct_best<-list(score=score,biological=i,technical=j)
    }
  list(
    factorNames=f_names,
    absCorrelationBiological=stats::setNames(biological,f_names),
    absCorrelationTechnical=stats::setNames(technical,f_names),
    maxBiological=unname(max(biological)),
    maxTechnical=unname(max(technical)),
    strongestBiologicalFactor=f_names[[best_biological]],
    strongestTechnicalFactor=f_names[[best_technical]],
    distinctTwoFactorRecovery=is.finite(distinct_best$score) &&
      biological[[distinct_best$biological]]>=biological_threshold &&
      technical[[distinct_best$technical]]>=batch_threshold,
    oneToOneMinAbsoluteCorrelation=if(is.finite(distinct_best$score))
      unname(distinct_best$score) else NA_real_,
    matchedIndependentFactorNames=if(is.finite(distinct_best$score))
      c(biological=f_names[[distinct_best$biological]],
        technical=f_names[[distinct_best$technical]]) else character(),
    biologicalRecovered=unname(max(biological)>=biological_threshold),
    batchSignalDetected=unname(max(technical)>=batch_threshold),
    batchDominated=unname(max(technical)>=batch_threshold &&
       max(technical)>max(biological)),
    limitation="Synthetic factor recovery only. Abs(correlation) is sign/permutation invariant for each truth; this does not establish factor identifiability or external biological meaning."
  )
}

# Descriptive null distribution from fully rerun training-only heldout DIABLO.
# Label permutation is at SUBJECT level before each outer split. No permutation
# result may reuse fitted training transforms or selected feature identities.
run_diablo_permutation_benchmark <- function(blocks,outcome,output_dir,
                 n_permutations=8L,permutation_seed=20261009L,
                 split_seeds=c(20261009L,20261010L),
                 tune=TRUE) {
  if(!is.numeric(n_permutations)||length(n_permutations)!=1L||
     n_permutations<8L || n_permutations>99L ||
     n_permutations!=as.integer(n_permutations))
    stop("Permutations must be a predeclared integer 8..99.")
  if(length(split_seeds)<2L || anyDuplicated(split_seeds))
    stop("Two or more prespecified independent outer split seeds required.")
  actual <- run_diablo_outer_holdout(blocks,outcome,
          file.path(output_dir,"observed"),seeds=split_seeds,tune=tune)
  observed <- vapply(actual$summary$signatureStability,
             function(x)x$meanPairwiseJaccard,numeric(1))
  original_ber<-actual$summary$meanBER
  if(!is.finite(original_ber)||any(!is.finite(observed)))
    stop("Observed heldout metrics unavailable.")
  # preserve names; use local PRNG seeds, not time/RNG state.
  shuffled <- outcome
  perms <- vector("list",n_permutations)
  for(i in seq_len(n_permutations)) {
    set.seed(as.integer(permutation_seed+i*97L))
    shuffled[]<-sample(as.character(outcome),length(outcome),replace=FALSE)
    permutation <- run_diablo_outer_holdout(blocks,shuffled,
       file.path(output_dir,sprintf("null_%03d",i)),
       seeds=split_seeds,tune=tune)
    null_j <- vapply(permutation$summary$signatureStability,
                   function(x)x$meanPairwiseJaccard,numeric(1))
    if(!identical(sort(names(null_j)),sort(names(observed))) ||
       any(!is.finite(null_j)) || !is.finite(permutation$summary$meanBER))
      stop("Null permutation omitted a block/fold or produced undefined metrics.")
    perms[[i]]<-list(index=i,meanBER=permutation$summary$meanBER,
                   jaccard=as.list(null_j))
  }
  null_ber<-vapply(perms,function(x)x$meanBER,numeric(1))
  jaccard_p <- lapply(names(observed),function(view) {
    values<-vapply(perms,function(x)x$jaccard[[view]],numeric(1))
    list(observed=unname(observed[[view]]),
         nullMean=mean(values),
         nullValues=as.numeric(values),
         upperTailEmpiricalP=(1+sum(values>=observed[[view]]))/
                              (1+n_permutations))
  })
  names(jaccard_p)<-names(observed)
  summary <- list(
    nPermutations=as.integer(n_permutations),
    permutationSeed=as.integer(permutation_seed),
    outerSplitSeeds=as.integer(split_seeds),
    tuningRepeatedForEveryPermutation=isTRUE(tune),
    observedMeanBER=original_ber,
    nullMeanBER=mean(null_ber),
    nullBERValues=as.numeric(null_ber),
    lowerTailBEREmpiricalP=(1+sum(null_ber<=original_ber))/
                           (1+n_permutations),
    featureStability=jaccard_p,
    minimumPossibleP=1/(1+n_permutations),
    status="exploratory_null_distribution_not_formal_calibration",
    note="All permutations refit train-only preprocessing, tuning and model. With 8 permutations, p cannot be below 1/9. Repeated splits share patients: no clinical, biomarker or confirmatory inference."
  )
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(summary=summary,permutations=perms,observed=actual$summary),
          file.path(output_dir,"diablo_permutation_null.rds"))
  summary
}
