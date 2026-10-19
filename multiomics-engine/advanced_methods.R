# Advanced multi-omics adapters for local/server execution.
# These functions intentionally use reference R implementations rather than
# relabelling browser-native PCA/PLS code as MOFA2 or DIABLO.
#
# Methodological policy:
# - deterministic seeds are recorded;
# - supervised multiblock models are tuned inside repeated CV before final fit;
# - unsupervised factor analysis refuses clearly underpowered sample counts;
# - longitudinal models attempt a random slope only when the design can support it;
# - sensitivity methods are compared, never merged by averaging p-values/q-values.

require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this advanced method.", pkg), call. = FALSE)
  }
}

validate_blocks <- function(blocks, min_samples = 3L) {
  if (!is.list(blocks) || length(blocks) < 2L) {
    stop("blocks must contain at least two omics matrices.", call. = FALSE)
  }
  if (is.null(names(blocks)) || any(names(blocks) == "")) {
    stop("Every block must have a name.", call. = FALSE)
  }
  blocks <- lapply(blocks, function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    if (is.null(rownames(x))) stop("Every block must use sample IDs as row names.", call. = FALSE)
    x
  })
  common <- Reduce(intersect, lapply(blocks, rownames))
  if (length(common) < as.integer(min_samples)) {
    stop(sprintf("At least %d samples must be shared across blocks.", as.integer(min_samples)), call. = FALSE)
  }
  lapply(blocks, function(x) x[common, , drop = FALSE])
}

write_matrix_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(sample_id = rownames(x), x, check.names = FALSE), path, row.names = FALSE)
}

filter_informative_features <- function(x, min_observed = 3L, max_features = NULL) {
  x <- as.matrix(x)
  observed <- colSums(is.finite(x))
  variances <- vapply(seq_len(ncol(x)), function(j) {
    values <- x[, j]
    values <- values[is.finite(values)]
    if (length(values) < min_observed) return(NA_real_)
    stats::var(values)
  }, numeric(1))
  keep <- observed >= min_observed & is.finite(variances) & variances > 0
  x <- x[, keep, drop = FALSE]
  variances <- variances[keep]
  if (!is.null(max_features) && is.finite(max_features) && ncol(x) > max_features) {
    order_idx <- order(variances, decreasing = TRUE, na.last = NA)
    x <- x[, order_idx[seq_len(as.integer(max_features))], drop = FALSE]
  }
  x
}

# Audit before create_mofa(). No silent intersection of unmatched subjects:
# otherwise the effective population can depend on which views were supplied.
mofa2_input_preflight <- function(blocks, min_samples=16L,
                                  max_missing_fraction=0.40) {
  if(!is.list(blocks) || length(blocks)<2L || is.null(names(blocks)) ||
     anyNA(names(blocks)) || any(!nzchar(names(blocks))) ||
     anyDuplicated(names(blocks)))
    stop("MOFA2 requires at least two uniquely named omics layers.")
  if(!is.finite(max_missing_fraction) || max_missing_fraction<0 ||
     max_missing_fraction>=1)
    stop("Invalid MOFA2 missingness safeguard.")
  views <- list()
  subject_sets <- list()
  for(layer in names(blocks)) {
    x <- as.matrix(blocks[[layer]])
    suppressWarnings(storage.mode(x)<-"double")
    ids <- rownames(x)
    features <- colnames(x)
    if(is.null(ids) || is.null(features) ||
       anyNA(ids) || anyNA(features) ||
       any(!nzchar(ids)) || any(!nzchar(features)) ||
       anyDuplicated(ids) || anyDuplicated(features))
      stop(sprintf("%s: unique nonempty sample and feature identifiers are required.",layer))
    if(nrow(x)<min_samples || ncol(x)<2L)
      stop(sprintf("%s: MOFA2 needs at least %s independent samples and two features.",
                   layer,min_samples))
    if(any(is.infinite(x)) || any(is.nan(x)))
      stop(sprintf("%s: infinite/NaN measurements are not valid missing values.",layer))
    missing <- mean(is.na(x))
    if(!is.finite(missing) || missing>max_missing_fraction)
      stop(sprintf("%s: excessive missing measurements (%.1f%%). Review the assay and missingness mechanism.",
                   layer,100*missing))
    subject_sets[[layer]] <- ids
    views[[layer]] <- list(
      samples=nrow(x),originalFeatures=ncol(x),
      missingFraction=missing,
      finiteFeatures=sum(colSums(is.finite(x))>=3L),
      note="Latent variance can reflect technical nuisance effects; MOFA2 factors are descriptive, not causal findings."
    )
  }
  reference <- subject_sets[[1]]
  if(any(vapply(subject_sets,function(ids)
       !setequal(ids,reference),logical(1))))
    stop("MOFA2 views have different sample sets. This backend does not model missing whole views: never silently discard unpaired subjects.")
  list(status="eligible_for_exploration",
       views=views,nSharedSubjects=length(reference),
       maxMissingFraction=max_missing_fraction,
       limitations="Same subject sets and low missingness do not establish MAR. Batch, biological identity, latent-factor number and stability need independent audit.")
}

# Input matrices throughout PMx are sample × feature; MOFA2's *list of
# matrices* interface is feature × sample. Preserve the explicit biological
# sample IDs as MOFA2 column names and align both views identically.
mofa2_orient_reference_matrices <- function(blocks) {
  mofa2_input_preflight(blocks)
  ids <- rownames(blocks[[1]])
  reference <- lapply(blocks,function(z) {
    z <- as.matrix(z)
    storage.mode(z) <- "double"
    t(z[ids,,drop=FALSE])
  })
  if(any(vapply(reference,function(z)
       !identical(colnames(z),ids),logical(1))))
    stop("MOFA2 biological subject orientation or order is inconsistent.")
  reference
}

run_mofa2_blocks <- function(
  blocks,
  output_dir,
  factors = 5L,
  seed = 20260924L,
  convergence_mode = "medium",
  min_samples = 16L,
  max_features_per_block = NULL
) {
  # Data validity must be established before the reference package is loaded,
  # so unsupported sample-set loss never appears as an installation error.
  input_audit <- mofa2_input_preflight(blocks,min_samples=min_samples)
  require_namespace("MOFA2")
  blocks <- validate_blocks(blocks, min_samples = min_samples)
  raw_feature_counts <- vapply(blocks, ncol, integer(1))
  blocks <- lapply(blocks, filter_informative_features, max_features = max_features_per_block)
  if (any(vapply(blocks, ncol, integer(1)) < 2L)) {
    stop("MOFA2 requires at least two informative non-constant features in every retained block.", call. = FALSE)
  }

  # MOFA2's own guidance stresses adequate normalisation, removal of known
  # technical effects and enough samples. The backend handles known technical
  # series upstream; here we additionally scale views so a high-variance assay
  # does not dominate merely because of measurement scale.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  n_samples <- nrow(blocks[[1]])
  factors <- max(1L, min(as.integer(factors), 10L, max(1L, n_samples - 2L)))

  reference_matrices <- mofa2_orient_reference_matrices(blocks)
  model <- MOFA2::create_mofa(reference_matrices)
  data_options <- MOFA2::get_default_data_options(model)
  model_options <- MOFA2::get_default_model_options(model)
  train_options <- MOFA2::get_default_training_options(model)
  if ("scale_views" %in% names(data_options)) data_options$scale_views <- TRUE
  model_options$num_factors <- factors
  train_options$seed <- as.integer(seed)
  train_options$convergence_mode <- convergence_mode
  train_options$verbose <- FALSE

  model <- MOFA2::prepare_mofa(
    model,
    data_options = data_options,
    model_options = model_options,
    training_options = train_options
  )
  trained <- MOFA2::run_mofa(
    model,
    outfile = file.path(output_dir, "mofa_model.hdf5"),
    use_basilisk = TRUE
  )

  factors_df <- MOFA2::get_factors(trained, factors = "all", as.data.frame = TRUE)
  weights_df <- MOFA2::get_weights(trained, views = "all", factors = "all", as.data.frame = TRUE)
  utils::write.csv(factors_df, file.path(output_dir, "mofa_factors.csv"), row.names = FALSE)
  utils::write.csv(weights_df, file.path(output_dir, "mofa_weights.csv"), row.names = FALSE)

  variance_explained <- try(MOFA2::calculate_variance_explained(trained), silent = TRUE)
  if (!inherits(variance_explained, "try-error")) {
    saveRDS(variance_explained, file.path(output_dir, "mofa_variance_explained.rds"))
  }

  summary <- list(
    method = "MOFA2",
    seed = as.integer(seed),
    samples = n_samples,
    factors_requested = as.integer(factors),
    inputAudit = input_audit,
    matrixOrientation = "feature_rows_subject_columns_for_MOFA2",
    implementationContract = "mofa2-feature-rows-sample-columns-v2",
    subjectIdentityPreserved = identical(colnames(reference_matrices[[1]]),rownames(blocks[[1]])),
    factorStatus = "descriptive_latent_covariance_not_validated_biomarker",
    blocks_before_filtering = raw_feature_counts,
    blocks_after_filtering = vapply(blocks, ncol, integer(1)),
    view_scaling = isTRUE(data_options$scale_views),
    guardrails = list(
      minimum_shared_samples = as.integer(min_samples),
      removed_constant_or_unobserved_features = TRUE,
      known_technical_effects_expected_to_be_handled_upstream = TRUE,
      interpretation = "Unsupervised latent factors describe covariance; they are not causal mechanisms or validated biomarkers."
    )
  )
  saveRDS(summary, file.path(output_dir, "mofa_summary.rds"))
  invisible(list(model = trained, summary = summary, variance_explained = variance_explained))
}

make_keepx_grid <- function(n_features) {
  n_features <- as.integer(n_features)
  if (n_features <= 1L) return(1L)
  candidates <- unique(pmin(n_features, c(2L, 5L, 10L, 20L, 50L, max(2L, round(n_features * 0.10)))))
  candidates <- sort(unique(as.integer(candidates[candidates >= 1L & candidates <= n_features])))
  if (length(candidates) < 2L) candidates <- sort(unique(c(1L, n_features)))
  candidates
}

normalise_keepx <- function(keepX, blocks, ncomp) {
  if (is.null(keepX)) {
    return(lapply(blocks, function(x) rep(min(20L, max(1L, ncol(x))), ncomp)))
  }
  out <- keepX
  for (view in names(blocks)) {
    values <- as.integer(out[[view]])
    if (!length(values)) values <- min(20L, max(1L, ncol(blocks[[view]])))
    values <- pmax(1L, pmin(values, ncol(blocks[[view]])))
    if (length(values) < ncomp) values <- c(values, rep(tail(values, 1), ncomp - length(values)))
    out[[view]] <- values[seq_len(ncomp)]
  }
  out[names(blocks)]
}

safe_perf_summary <- function(perf) {
  if (is.null(perf) || inherits(perf, "try-error")) return(NULL)
  error_rate <- perf$error.rate
  if (is.null(error_rate)) {
    return(list(
      status = "available",
      note = "Performance object saved; compact BER extraction unavailable for this mixOmics version. Internal cross-validation is not external validation."
    ))
  }
  ber <- NULL
  if (is.list(error_rate) && !is.null(error_rate$BER)) ber <- error_rate$BER
  if (is.null(ber) && is.list(error_rate) && !is.null(error_rate$WeightedVote) && !is.null(error_rate$WeightedVote$BER)) ber <- error_rate$WeightedVote$BER
  numeric_ber <- suppressWarnings(as.numeric(unlist(ber)))
  numeric_ber <- numeric_ber[is.finite(numeric_ber)]
  list(
    status = "available",
    ber_values_extracted = length(numeric_ber),
    ber_median = if (length(numeric_ber)) stats::median(numeric_ber) else NA_real_,
    ber_mean = if (length(numeric_ber)) mean(numeric_ber) else NA_real_,
    ber_min = if (length(numeric_ber)) min(numeric_ber) else NA_real_,
    ber_max = if (length(numeric_ber)) max(numeric_ber) else NA_real_,
    note = "These BER values summarise repeated internal cross-validation on the same cohort used for model tuning. They are descriptive internal-validation diagnostics, not an unbiased external performance estimate."
  )
}


# Strict subject-level outer evaluation, separate from tuning/perf() on the
# development cohort. Only accepts log-valued as-supplied matrices; all filters,
# imputation, centering and standardisation are fit on TRAINING samples only.
# This is a repeated internal held-out validation, NOT external replication.
prepare_diablo_holdout_split <- function(blocks, outcome, seed=20261009L,
                                         test_fraction=0.25,
                                         max_features_per_block=500L) {
  if(!is.list(blocks) || length(blocks)<2L || is.null(names(blocks)) ||
     any(!nzchar(names(blocks))) || anyDuplicated(names(blocks)))
    stop("Holdout DIABLO needs >=2 uniquely named omic blocks.")
  blocks <- lapply(blocks,function(z) {
    z <- as.matrix(z)
    storage.mode(z) <- "double"
    if(is.null(rownames(z))||anyNA(rownames(z))||anyDuplicated(rownames(z)) ||
       any(!nzchar(rownames(z))) || is.null(colnames(z))||
       anyDuplicated(colnames(z)) || any(!nzchar(colnames(z))))
      stop("Holdout DIABLO requires unique subject/sample and feature identifiers.")
    z
  })
  common <- Reduce(intersect,lapply(blocks,rownames))
  if(length(common)<32L)
    stop("Holdout DIABLO requires >=32 shared independent subjects.")
  blocks <- lapply(blocks,function(z)z[common,,drop=FALSE])
  if(is.null(names(outcome))||anyDuplicated(names(outcome)) ||
     anyNA(outcome[common])||any(!nzchar(as.character(outcome[common]))))
    stop("Holdout labels must be uniquely named, complete and match every shared subject.")
  y <- droplevels(factor(outcome[common]))
  if(nlevels(y)!=2L || min(table(y))<12L)
    stop("Holdout DIABLO currently supports 2 classes with >=12 independent subjects per class.")
  if(!is.finite(test_fraction)||test_fraction<0.15||test_fraction>0.40)
    stop("Holdout test fraction outside supported range.")
  set.seed(as.integer(seed))
  test <- unlist(lapply(split(seq_along(common),y),function(ix) {
    size <- max(3L,as.integer(round(length(ix)*test_fraction)))
    sample(ix,size=size,replace=FALSE)
  }),use.names=FALSE)
  train <- setdiff(seq_along(common),test)
  if(length(intersect(train,test))||length(unique(c(train,test)))!=length(common))
    stop("Train and holdout subjects are not disjoint and exhaustive.")
  if(min(table(y[train]))<8L || min(table(y[test]))<3L)
    stop("Insufficient class support after subject-level partition.")
  recipes <- list()
  train_blocks <- list()
  test_blocks <- list()
  for(layer in names(blocks)) {
    z <- blocks[[layer]]
    tr <- z[train,,drop=FALSE]
    te <- z[test,,drop=FALSE]
    # Fold-specific feature selection: no full-cohort variance or missingness.
    observed <- colSums(is.finite(tr))
    missing_fraction <- 1-observed/nrow(tr)
    variance <- vapply(seq_len(ncol(tr)),function(j) {
      v <- tr[,j];v <- v[is.finite(v)]
      if(length(v)<3L)return(NA_real_)
      stats::var(v)
    },numeric(1))
    available <- which(observed>=8L & missing_fraction<=0.25 &
                         is.finite(variance)&variance>1e-12)
    if(length(available)<2L)
      stop(paste("Insufficient train-only informative features:",layer))
    if(length(available)>max_features_per_block)
      available<-available[order(variance[available],decreasing=TRUE,
                               na.last=NA)[seq_len(max_features_per_block)]]
    tr <- tr[,available,drop=FALSE]
    te <- te[,available,drop=FALSE]
    medians <- apply(tr,2,stats::median,na.rm=TRUE)
    for(j in seq_len(ncol(tr))) {
      tr[!is.finite(tr[,j]),j] <- medians[j]
      te[!is.finite(te[,j]),j] <- medians[j]
    }
    means <- colMeans(tr)
    sd_train <- apply(tr,2,stats::sd)
    if(any(!is.finite(sd_train))||any(sd_train<=1e-12))
      stop("Non-identifiable training-only scaling.")
    tr <- sweep(sweep(tr,2,means,"-"),2,sd_train,"/")
    te <- sweep(sweep(te,2,means,"-"),2,sd_train,"/")
    if(any(!is.finite(tr))||any(!is.finite(te)))
      stop("Non-finite fold-local omics values.")
    train_blocks[[layer]] <- tr
    test_blocks[[layer]] <- te
    recipes[[layer]] <- list(
      featureNames=colnames(tr),
      medians=medians,
      means=means,
      standardDeviations=sd_train,
      fittedOn="training_subjects_only"
    )
  }
  names_y <- as.character(y)
  list(
    train=train_blocks,test=test_blocks,
    y_train=factor(names_y[train],levels=levels(y)),
    y_test=factor(names_y[test],levels=levels(y)),
    trainIds=common[train],testIds=common[test],
    trainIndices=train,testIndices=test,
    recipes=recipes,
    seed=as.integer(seed),
    boundary="Repeated cohort holdout with fold-local preprocessing, not independent external validation."
  )
}

# Compare selected feature identities, never their unstable sign/order or scores.
# Repeated train/test splits overlap; this is a descriptive repeatability measure
# and not an independent sampling-based confidence interval.
diablo_selection_stability <- function(selected_by_split) {
  if(!is.list(selected_by_split)||length(selected_by_split)<2L)
    stop("DIABLO stability requires >=2 successful predeclared folds.")
  blocks<-names(selected_by_split[[1]])
  if(is.null(blocks)||length(blocks)<2L||anyDuplicated(blocks)||
     any(!nzchar(blocks)) ||
     any(vapply(selected_by_split,function(fold)
       !setequal(names(fold),blocks),logical(1))))
    stop("Different omics blocks between DIABLO stability folds.")
  out<-list()
  for(view in blocks) {
    sets<-lapply(selected_by_split,function(fold) {
      vals<-as.character(fold[[view]])
      if(!length(vals)||anyNA(vals)||any(!nzchar(vals)))
        stop(paste("Missing selected features for",view))
      unique(vals)
    })
    pairs<-utils::combn(seq_along(sets),2L)
    jaccard<-vapply(seq_len(ncol(pairs)),function(j) {
      a<-sets[[pairs[1,j]]]; b<-sets[[pairs[2,j]]]
      length(intersect(a,b))/length(union(a,b))
    },numeric(1))
    tab<-sort(table(unlist(sets,use.names=FALSE)),decreasing=TRUE)
    out[[view]]<-list(
      meanPairwiseJaccard=unname(mean(jaccard)),
      minPairwiseJaccard=unname(min(jaccard)),
      maxPairwiseJaccard=unname(max(jaccard)),
      nPairwiseComparisons=length(jaccard),
      nDistinctSelected=length(tab),
      nSelectedPerFold=as.integer(lengths(sets)),
      frequency=as.list(stats::setNames(as.numeric(tab)/length(sets),names(tab))),
      repeatedlySelected=as.character(names(tab)[as.numeric(tab)>=2L]),
      limitation="Repeated folds share original subjects; selection frequencies and Jaccard are internal exploratory stability, not externally validated biomarkers."
    )
  }
  out
}

# mixOmics::selectVar() returns named variables in version-dependent slots.
# Fail rather than counting all feature loadings or unmatched identifiers.
selected_diablo_features <- function(model,view,training_features) {
  # For block.splsda, selectVar(model, comp=1) returns a list indexed
  # by VIEW. Passing block=view can still return that outer list; reading
  # $name on the outer object silently fails in affected mixOmics releases.
  selected <- mixOmics::selectVar(model,comp=1L)
  if(!is.list(selected) || is.null(selected[[view]]))
    stop(paste("mixOmics selectVar lacks block:",view))
  choice <- selected[[view]]
  ids <- as.character(unlist(choice$name,use.names=FALSE))
  if(!length(ids) && (is.matrix(choice$value)||is.data.frame(choice$value)))
    ids <- rownames(choice$value)
  ids <- unique(as.character(ids))
  if(!length(ids)||anyNA(ids)||any(!nzchar(ids)) ||
     any(!ids %in% training_features))
    stop(paste("mixOmics selected features cannot be matched to train-only names:",view))
  ids
}

# The model is refit inside each outer split; any tuning is restricted to
# training subjects. Invalid folds fail the WHOLE benchmark, never silently
# reduce the denominator and make performance look better.
run_diablo_outer_holdout <- function(blocks,outcome,output_dir,
                                     seeds=c(20261009L,20261010L,20261011L),
                                     tune=TRUE) {
  require_namespace("mixOmics")
  if(length(seeds)<2L||anyDuplicated(seeds))stop("Need >=2 distinct predeclared outer seeds.")
  all_rows <- list()
  summaries <- list()
  selections <- list()
  for(index in seq_along(seeds)) {
    split <- prepare_diablo_holdout_split(blocks,outcome,seed=seeds[index])
    ncomp <- 1L
    design <- matrix(0.1,length(split$train),length(split$train),
      dimnames=list(names(split$train),names(split$train)))
    diag(design) <- 0
    keep <- lapply(split$train,function(z)min(10L,ncol(z)))
    names(keep) <- names(split$train)
    if(isTRUE(tune)) {
      grid <- lapply(split$train,function(z)sort(unique(pmin(ncol(z),c(2L,5L,10L)))))
      if(all(lengths(grid)>=2L)) {
        set.seed(as.integer(seeds[index]))
        trained <- mixOmics::tune.block.splsda(
          X=split$train,Y=split$y_train,ncomp=ncomp,
          test.keepX=grid,design=design,
          validation="Mfold",folds=3L,nrepeat=2L,measure="BER",
          dist="max.dist",progressBar=FALSE,seed=as.integer(seeds[index])
        )
        if(!is.null(trained$choice.keepX)) {
          keep <- lapply(names(split$train),function(name)
            as.integer(trained$choice.keepX[[name]][1]))
          names(keep) <- names(split$train)
        } else stop("Training-only DIABLO tuning did not return a keepX choice.")
      }
    }
    keep <- lapply(keep,function(k)rep(as.integer(k[1]),ncomp))
    set.seed(as.integer(seeds[index]))
    model <- mixOmics::block.splsda(
      X=split$train,Y=split$y_train,ncomp=ncomp,
      keepX=keep,design=design,scale=FALSE
    )
    selections[[index]] <- lapply(names(split$train),function(view)
      selected_diablo_features(model,view,colnames(split$train[[view]])))
    names(selections[[index]]) <- names(split$train)
    predicted <- predict(model,newdata=split$test,dist="max.dist")
    vote <- predicted$WeightedVote[["max.dist"]]
    if(is.null(vote) || nrow(as.matrix(vote))!=length(split$testIds))
      stop("DIABLO weighted vote missing or wrong holdout subject count.")
    classes <- as.character(as.matrix(vote)[,ncomp])
    if(length(classes)!=length(split$testIds) ||
       anyNA(classes)||any(!classes%in%levels(split$y_test)))
      stop("Invalid/ambiguous DIABLO holdout predictions, no silent omissions.")
    actual <- as.character(split$y_test)
    counts <- table(split$y_test)
    recall <- vapply(levels(split$y_test),function(label)
      mean(classes[actual==label]==label),numeric(1))
    ber <- mean(1-recall)
    all_rows[[index]] <- data.frame(
      seed=seeds[index],sampleId=split$testIds,
      actual=actual,predicted=classes,correct=as.integer(actual==classes),
      stringsAsFactors=FALSE
    )
    summaries[[index]] <- list(seed=seeds[index],
      trainingSubjects=length(split$trainIds),
      heldOutSubjects=length(split$testIds),
      balancedErrorRate=ber,
      accuracy=mean(actual==classes),
      classRecall=as.list(recall),
      trainingOnlyTuning=isTRUE(tune))
  }
  predictions <- do.call(rbind,all_rows)
  if(nrow(predictions)!=sum(vapply(summaries,function(x)x$heldOutSubjects,integer(1))))
    stop("DIABLO holdout missing predictions.")
  stability <- diablo_selection_stability(selections)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(predictions,file.path(output_dir,"diablo_outer_holdout.csv"),row.names=FALSE)
  saveRDS(list(versions="internal-fold-feature-selections-v1",
      selectedBySplit=selections,stability=stability),
      file.path(output_dir,"diablo_signature_stability.rds"))
  list(status="ok",summary=list(
    signatureStability=stability,
    signatureStatus="internal_repeatability_only_not_validated_biomarker",
    evaluation="Repeated subject-level heldout, with training-only filtering/imputation/scaling and training-only keepX tuning",
    seeds=as.integer(seeds),folds=summaries,
    meanBER=mean(vapply(summaries,function(x)x$balancedErrorRate,numeric(1))),
    maxBER=max(vapply(summaries,function(x)x$balancedErrorRate,numeric(1))),
    trainingPipeline="as-supplied log-scale omics, no varying batch or nuisance covariates; observed samples must be truly independent",
    caveat="Held-out subsets reuse the same source cohort across splits: not external replication. Outcome-related pre-import processing, MS global QC and study selection remain outside this validation."
  ),predictions=predictions)
}

run_diablo_blocks <- function(
  blocks,
  outcome,
  output_dir,
  ncomp = 2L,
  keepX = NULL,
  seed = 20260924L,
  tune = TRUE,
  nrepeat = 5L,
  folds = NULL
) {
  require_namespace("mixOmics")
  blocks <- validate_blocks(blocks, min_samples = 6L)
  blocks <- lapply(blocks, filter_informative_features)
  if (any(vapply(blocks, ncol, integer(1)) < 2L)) {
    stop("DIABLO requires at least two informative features in every block.", call. = FALSE)
  }

  common <- rownames(blocks[[1]])
  if (is.null(names(outcome))) {
    if (length(outcome) != length(common)) stop("Unnamed outcome must have one value per shared sample.", call. = FALSE)
    y <- outcome
  } else {
    y <- outcome[common]
  }
  if (any(is.na(y)) || any(!nzchar(as.character(y)))) {
    stop("Outcome contains missing values for shared samples.", call. = FALSE)
  }
  y <- droplevels(factor(y))
  if (nlevels(y) < 2L) stop("DIABLO requires at least two outcome classes.", call. = FALSE)
  class_counts <- table(y)
  min_class <- min(class_counts)
  if (min_class < 3L) {
    stop("DIABLO tuning requires at least three samples in every class; use an exploratory non-supervised route instead.", call. = FALSE)
  }

  if (is.null(folds)) folds <- min(5L, as.integer(min_class))
  folds <- as.integer(max(3L, min(folds, min_class)))
  # The automatic backend previously requested three repeats. Enforce at least
  # five here so both direct calls and browser-triggered runs receive a more
  # stable repeated-CV diagnostic without relying on the caller.
  nrepeat <- as.integer(max(5L, nrepeat))
  max_components <- max(1L, min(3L, nlevels(y), length(y) - 2L))
  requested_components <- as.integer(max(1L, min(ncomp, max_components)))

  design <- matrix(
    0.1,
    nrow = length(blocks),
    ncol = length(blocks),
    dimnames = list(names(blocks), names(blocks))
  )
  diag(design) <- 0

  set.seed(seed)
  tuning <- NULL
  tuning_error <- NULL
  tuned <- FALSE
  final_ncomp <- requested_components
  final_keepX <- normalise_keepx(keepX, blocks, final_ncomp)
  test_keepX <- NULL

  if (isTRUE(tune) && is.null(keepX)) {
    test_keepX <- lapply(blocks, function(x) make_keepx_grid(ncol(x)))
    if (all(vapply(test_keepX, length, integer(1)) >= 2L)) {
      tuning <- try(
        mixOmics::tune.block.splsda(
          X = blocks,
          Y = y,
          ncomp = max_components,
          test.keepX = test_keepX,
          validation = "Mfold",
          folds = folds,
          dist = "max.dist",
          measure = "BER",
          nrepeat = nrepeat,
          design = design,
          progressBar = FALSE,
          seed = as.integer(seed)
        ),
        silent = TRUE
      )
      if (!inherits(tuning, "try-error")) {
        tuned_ncomp <- suppressWarnings(as.integer(tuning$choice.ncomp$ncomp))
        if (length(tuned_ncomp) && is.finite(tuned_ncomp) && tuned_ncomp >= 1L) {
          final_ncomp <- min(tuned_ncomp, max_components)
        }
        if (!is.null(tuning$choice.keepX)) {
          final_keepX <- normalise_keepx(tuning$choice.keepX, blocks, final_ncomp)
          tuned <- TRUE
        }
      } else {
        tuning_error <- as.character(tuning)
      }
    } else {
      tuning_error <- "Automatic sparsity tuning was not attempted because at least one block had fewer than two admissible keepX candidates."
    }
  }

  fit <- mixOmics::block.splsda(
    X = blocks,
    Y = y,
    ncomp = final_ncomp,
    keepX = final_keepX,
    design = design
  )

  performance_seed <- as.integer(seed) + 1L
  performance <- try(
    mixOmics::perf(
      fit,
      validation = "Mfold",
      folds = folds,
      nrepeat = nrepeat,
      progressBar = FALSE,
      seed = performance_seed
    ),
    silent = TRUE
  )
  performance_error <- if (inherits(performance, "try-error")) as.character(performance) else NULL

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, file.path(output_dir, "diablo_model.rds"))
  if (!is.null(tuning) && !inherits(tuning, "try-error")) saveRDS(tuning, file.path(output_dir, "diablo_tuning.rds"))
  if (!inherits(performance, "try-error")) saveRDS(performance, file.path(output_dir, "diablo_performance.rds"))

  variates <- do.call(cbind, lapply(names(fit$variates), function(view) {
    x <- fit$variates[[view]]
    colnames(x) <- paste0(view, "_comp", seq_len(ncol(x)))
    x
  }))
  write_matrix_csv(variates, file.path(output_dir, "diablo_variates.csv"))

  selected <- lapply(names(blocks), function(view) {
    lapply(seq_len(final_ncomp), function(comp) mixOmics::selectVar(fit, block = view, comp = comp)$value)
  })
  names(selected) <- names(blocks)
  saveRDS(selected, file.path(output_dir, "diablo_selected_variables.rds"))

  summary <- list(
    method = "mixOmics DIABLO / block.splsda",
    seed = as.integer(seed),
    performance_seed = performance_seed,
    samples = length(y),
    classes = as.list(class_counts),
    blocks = vapply(blocks, ncol, integer(1)),
    tuning = list(
      requested = isTRUE(tune),
      completed = tuned,
      validation = "repeated M-fold CV; folds bounded by the smallest class",
      folds = folds,
      repeats = nrepeat,
      measure = "BER",
      candidate_keepX = test_keepX,
      ncomp = final_ncomp,
      keepX = final_keepX,
      error = tuning_error,
      note = if (tuned) "Final model uses repeated cross-validated sparsity choices." else "Automatic tuning was unavailable; conservative keepX values were used and this model should be treated as exploratory."
    ),
    internal_performance = safe_perf_summary(performance),
    internal_performance_error = performance_error,
    validation_boundary = "Hyperparameter tuning and the reported perf() diagnostic use the same cohort. Repeated CV reduces split sensitivity but does not create an independent validation cohort; external validation remains required for generalisable biomarker or prediction claims.",
    interpretation = "Selected variables are a supervised multivariate signature. Their stability and external validity must be assessed before biomarker claims."
  )
  saveRDS(summary, file.path(output_dir, "diablo_summary.rds"))
  invisible(list(model = fit, summary = summary, tuning = tuning, performance = performance))
}

# Narrow R1 *eligibility* preflight: two independent groups, an additive
# full-rank design and raw nonnegative integer counts. This intentionally
# DOES NOT certify data provenance, statistical power, differential expression,
# false-discovery control or any biological conclusion. No silent sample
# intersection, pseudoreplication or feature-ID repair in this path.
rnaseq_confirmatory_preflight <- function(counts, metadata, design_formula,
                                          contrast,
                                          min_subjects_per_group=3L) {
  if(!is.matrix(counts)||!is.numeric(counts)||
     nrow(counts)<6L||ncol(counts)<20L||
     is.null(rownames(counts))||is.null(colnames(counts))||
     anyNA(rownames(counts))||anyNA(colnames(counts))||
     any(!nzchar(rownames(counts)))||any(!nzchar(colnames(counts)))||
     anyDuplicated(rownames(counts))||anyDuplicated(colnames(counts)))
    stop("R1 RNA-seq requires >=6 biological samples, >=20 features and unique nonempty input IDs")
  if(any(!is.finite(counts))||any(counts<0)||
     any(abs(counts-round(counts))>1e-8)||any(rowSums(counts)<=0))
    stop("R1 RNA-seq requires finite, non-negative, unnormalised integer counts and nonempty sample libraries")
  if(!is.data.frame(metadata)||is.null(rownames(metadata))||
     anyNA(rownames(metadata))||anyDuplicated(rownames(metadata))||
     any(!nzchar(rownames(metadata)))||
     !setequal(rownames(metadata),rownames(counts))||
     nrow(metadata)!=nrow(counts))
    stop("R1 RNA-seq sample IDs must match metadata exactly; never drop or multiply biological samples")
  if(!all(c("subject_id","condition")%in%names(metadata))||
     anyNA(metadata$subject_id)||any(!nzchar(trimws(as.character(metadata$subject_id))))||
     anyDuplicated(as.character(metadata$subject_id)))
    stop("R1 RNA-seq needs one unique, known biological subject per sample (no pseudoreplicates)")
  if(!inherits(design_formula,"formula")||
     length(design_formula)!=2L||
     attr(stats::terms(design_formula),"intercept")!=1L||
     any(attr(stats::terms(design_formula),"order")!=1L))
    stop("R1 RNA-seq only supports an additive formula with an intercept (no interactions, pairing or offsets)")
  fields<-all.vars(design_formula)
  labels<-attr(stats::terms(design_formula),"term.labels")
  if(!length(labels)||any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$",labels))||
     !"condition"%in%labels||anyDuplicated(labels)||
     any(!fields%in%names(metadata)))
    stop("R1 RNA-seq needs a named condition main effect and plain additive measured covariates")
  m<-metadata[match(rownames(counts),rownames(metadata)),,drop=FALSE]
  for(field in fields){
    value<-as.character(m[[field]])
    if(anyNA(value)||any(!nzchar(trimws(value))))
      stop(paste("R1 RNA-seq has missing or blank design covariate:",field))
  }
  groups<-unique(as.character(m$condition))
  if(length(groups)!=2L||any(!nzchar(groups)))
    stop("R1 RNA-seq requires exactly two nonempty independent conditions")
  n_by_group<-table(as.character(m$condition))
  if(length(min_subjects_per_group)!=1L||!is.finite(min_subjects_per_group)||
     min_subjects_per_group!=as.integer(min_subjects_per_group)||
     min_subjects_per_group<3L||
     any(n_by_group<min_subjects_per_group))
    stop("R1 RNA-seq has too few independent biological subjects per group (minimum 3, not a power claim)")
  if(!is.character(contrast)||length(contrast)!=3L||
     !identical(contrast[[1L]],"condition")||
     identical(contrast[[2L]],contrast[[3L]])||
     !setequal(contrast[2:3],groups))
    stop("R1 RNA-seq requires one explicit two-level condition contrast, not an implicit coefficient")
  batch<-if("batch"%in%names(m))trimws(as.character(m$batch)) else character()
  observed_batch<-unique(batch[!is.na(batch)&nzchar(batch)])
  if(length(observed_batch)>1L && !"batch"%in%labels)
    stop("R1 RNA-seq detects varying known batches excluded from the additive design")
  if("batch"%in%labels && (anyNA(batch)||any(!nzchar(batch))))
    stop("R1 RNA-seq cannot adjust a partially missing technical batch")
  if("condition"%in%labels && !"condition"%in%names(m))
    stop("Missing condition in confirmatory design")
  for(field in fields) {
    if(is.character(m[[field]]))m[[field]]<-factor(m[[field]])
  }
  matrix_design<-tryCatch(stats::model.matrix(design_formula,data=m,
                                               na.action=stats::na.fail),
                          error=function(e)stop("R1 RNA-seq design cannot be constructed: ",conditionMessage(e)))
  if(nrow(matrix_design)!=nrow(counts)||qr(matrix_design)$rank!=ncol(matrix_design)||
     nrow(matrix_design)-ncol(matrix_design)<2L)
    stop("R1 RNA-seq design is non-identifiable or has insufficient residual degrees of freedom")
  list(
    status="eligible_for_reference_fitting_only",
    cohort_samples=nrow(counts),measured_features=ncol(counts),
    independent_subjects=length(unique(m$subject_id)),
    group_sizes=as.list(n_by_group),
    rank=qr(matrix_design)$rank,model_columns=ncol(matrix_design),
    residual_df=nrow(matrix_design)-ncol(matrix_design),
    condition_contrast=unname(contrast),covariates=labels[labels!="condition"],
    missing_or_unverified_batch=length(observed_batch)==0L,
    rejected_silent_subject_loss=TRUE,
    certified=FALSE,
    limitations=c(
      "Integer-looking input is not proof of genuine raw sequencer counts.",
      "At least three subjects per group is a software guard, NOT a statistical power guarantee.",
      "Batch measurement, library and mapping QC, biological independence, protocol preregistration, multiple-testing family and external scientific review remain necessary.",
      "Passing preflight is not confirmatory inference, independent replication or clinical certification."
    )
  )
}

run_deseq2_counts <- function(
  counts,
  metadata,
  output_dir,
  design_formula = ~ condition,
  contrast = NULL,
  alpha = 0.05,
  confirmatory_contract = FALSE
) {
  input_contract <- if(isTRUE(confirmatory_contract))
    rnaseq_confirmatory_preflight(as.matrix(counts), as.data.frame(metadata),
                                 design_formula, contrast)
    else NULL
  require_namespace("DESeq2")
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  if (is.null(rownames(counts))) stop("counts must use sample IDs as row names.", call. = FALSE)
  metadata <- as.data.frame(metadata)
  if (is.null(rownames(metadata))) stop("metadata must use sample IDs as row names.", call. = FALSE)
  common <- intersect(rownames(counts), rownames(metadata))
  if (length(common) < 3L) stop("At least three matched samples are required.", call. = FALSE)
  counts <- counts[common, , drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(abs(counts - round(counts)) > 1e-8)) {
    stop("DESeq2 requires non-negative, finite, unnormalised integer counts. Do not silently round normalised or fractional abundances.", call. = FALSE)
  }
  design <- stats::model.matrix(design_formula, data=metadata)
  if (qr(design)$rank < ncol(design)) stop("DESeq2 design is rank-deficient.", call. = FALSE)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = t(counts),
    colData = metadata,
    design = design_formula
  )
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= max(2L, ceiling(ncol(dds) * 0.2))
  if (!any(keep)) stop("Expression filtering removed all genes.", call. = FALSE)
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  res <- if (is.null(contrast)) {
    DESeq2::results(dds, alpha = alpha)
  } else {
    DESeq2::results(dds, contrast = contrast, alpha = alpha)
  }
  result <- data.frame(feature = rownames(res), as.data.frame(res), row.names = NULL, check.names = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)
  saveRDS(dds, file.path(output_dir, "deseq2_model.rds"))
  invisible(list(model = dds, results = result,
                 confirmatory_input_contract = input_contract))
}

run_edger_ql_counts <- function(
  counts,
  metadata,
  output_dir,
  design_formula = ~ condition,
  coefficient = NULL,
  robust = TRUE
) {
  require_namespace("edgeR")
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  metadata <- as.data.frame(metadata)
  if (is.null(rownames(counts)) || is.null(rownames(metadata))) {
    stop("counts and metadata must use sample IDs as row names.", call. = FALSE)
  }
  common <- intersect(rownames(counts), rownames(metadata))
  if (length(common) < 3L) stop("At least three matched samples are required.", call. = FALSE)
  x <- counts[common, , drop = FALSE]
  meta <- metadata[common, , drop = FALSE]
  design <- stats::model.matrix(design_formula, data = meta)
  if (qr(design)$rank < ncol(design)) stop("The edgeR design matrix is rank-deficient.", call. = FALSE)

  if (any(!is.finite(x)) || any(x < 0)) {
    stop("RNA-seq counts must be finite and non-negative.", call. = FALSE)
  }
  # edgeR and voom can process estimated fractional counts directly; never
  # silently round them as that modifies the input and its mean/variance.
  y <- edgeR::DGEList(counts = t(x))
  keep <- edgeR::filterByExpr(y, design = design)
  if (!any(keep)) stop("edgeR::filterByExpr removed all genes.", call. = FALSE)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design = design, robust = robust)
  fit <- edgeR::glmQLFit(y, design = design, robust = robust)
  if (is.null(coefficient)) coefficient <- ncol(design)
  test <- edgeR::glmQLFTest(fit, coef = coefficient)
  table <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table
  result <- data.frame(feature = rownames(table), table, row.names = NULL, check.names = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(output_dir, "edger_ql_results.csv"), row.names = FALSE)
  saveRDS(list(dge = y, fit = fit, test = test, design = design), file.path(output_dir, "edger_ql_model.rds"))
  invisible(list(model = fit, test = test, results = result, design = design, dge = y))
}

run_voom_counts <- function(
  counts,
  metadata,
  output_dir,
  design_formula = ~ condition,
  coefficient = NULL,
  robust = TRUE,
  trend = FALSE
) {
  require_namespace("edgeR")
  require_namespace("limma")
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  metadata <- as.data.frame(metadata)
  if (is.null(rownames(counts)) || is.null(rownames(metadata))) {
    stop("counts and metadata must use sample IDs as row names.", call. = FALSE)
  }
  common <- intersect(rownames(counts), rownames(metadata))
  if (length(common) < 3L) stop("At least three matched samples are required.", call. = FALSE)
  x <- counts[common, , drop = FALSE]
  meta <- metadata[common, , drop = FALSE]
  design <- stats::model.matrix(design_formula, data = meta)
  if (qr(design)$rank < ncol(design)) stop("The voom design matrix is rank-deficient.", call. = FALSE)

  if (any(!is.finite(x)) || any(x < 0)) {
    stop("RNA-seq counts must be finite and non-negative.", call. = FALSE)
  }
  # edgeR and voom can process estimated fractional counts directly; never
  # silently round them as that modifies the input and its mean/variance.
  y <- edgeR::DGEList(counts = t(x))
  keep <- edgeR::filterByExpr(y, design = design)
  if (!any(keep)) stop("edgeR::filterByExpr removed all genes.", call. = FALSE)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  voom <- limma::voom(y, design = design, plot = FALSE)
  fit <- limma::lmFit(voom, design)
  fit <- limma::eBayes(fit, robust = robust, trend = trend)
  if (is.null(coefficient)) coefficient <- ncol(design)
  result <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "P")
  result <- data.frame(feature = rownames(result), result, row.names = NULL, check.names = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(output_dir, "voom_results.csv"), row.names = FALSE)
  utils::write.csv(data.frame(feature = rownames(voom$E), voom$E, check.names = FALSE), file.path(output_dir, "voom_logcpm.csv"), row.names = FALSE)
  saveRDS(list(dge = y, voom = voom, fit = fit, design = design), file.path(output_dir, "voom_model.rds"))
  invisible(list(model = fit, results = result, design = design, voom = voom, dge = y))
}

compare_rnaseq_methods <- function(deseq_results, voom_results, alpha = 0.05) {
  a <- as.data.frame(deseq_results, stringsAsFactors = FALSE)
  b <- as.data.frame(voom_results, stringsAsFactors = FALSE)
  if (!all(c("feature", "stat", "padj") %in% names(a)) ||
      !all(c("feature", "t", "adj.P.Val") %in% names(b))) {
    return(list(status = "unavailable", message = "Required DESeq2/voom result columns were not available."))
  }
  merged <- merge(
    a[, c("feature", "stat", "padj")],
    b[, c("feature", "t", "adj.P.Val")],
    by = "feature",
    all = FALSE
  )
  names(merged) <- c("feature", "deseq2_stat", "deseq2_q", "voom_t", "voom_q")
  valid <- is.finite(merged$deseq2_stat) & is.finite(merged$voom_t)
  rho <- if (sum(valid) >= 3L) suppressWarnings(stats::cor(merged$deseq2_stat[valid], merged$voom_t[valid], method = "spearman")) else NA_real_
  sig_a <- is.finite(merged$deseq2_q) & merged$deseq2_q <= alpha
  sig_b <- is.finite(merged$voom_q) & merged$voom_q <= alpha
  union_sig <- sig_a | sig_b
  sign_agreement <- if (any(union_sig)) mean(sign(merged$deseq2_stat[union_sig]) == sign(merged$voom_t[union_sig]), na.rm = TRUE) else NA_real_
  overlap <- sum(sig_a & sig_b)
  union_n <- sum(union_sig)
  jaccard <- if (union_n > 0L) overlap / union_n else NA_real_
  list(
    status = "ok",
    common_features = nrow(merged),
    spearman_statistics = unname(rho),
    deseq2_significant = sum(sig_a),
    voom_significant = sum(sig_b),
    significant_overlap = overlap,
    significant_jaccard = jaccard,
    sign_agreement_among_any_significant = sign_agreement,
    alpha = alpha,
    interpretation = "Concordance is a sensitivity analysis. q-values from distinct methods are not averaged or combined."
  )
}

run_limma_matrix <- function(
  matrix,
  metadata,
  output_dir,
  design_formula = ~ condition,
  coefficient = NULL,
  trend = TRUE,
  robust = TRUE
) {
  require_namespace("limma")
  matrix <- as.matrix(matrix)
  storage.mode(matrix) <- "double"
  metadata <- as.data.frame(metadata)
  if (is.null(rownames(matrix)) || is.null(rownames(metadata))) {
    stop("matrix and metadata must use sample IDs as row names.", call. = FALSE)
  }
  common <- intersect(rownames(matrix), rownames(metadata))
  if (length(common) < 3L) stop("At least three matched samples are required.", call. = FALSE)
  x <- matrix[common, , drop = FALSE]
  meta <- metadata[common, , drop = FALSE]
  design <- stats::model.matrix(design_formula, data = meta)
  if (qr(design)$rank < ncol(design)) stop("The limma design matrix is rank-deficient.", call. = FALSE)

  fit <- limma::lmFit(t(x), design)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  if (is.null(coefficient)) coefficient <- ncol(design)
  result <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "P")
  result <- data.frame(feature = rownames(result), result, row.names = NULL, check.names = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(output_dir, "limma_results.csv"), row.names = FALSE)
  saveRDS(fit, file.path(output_dir, "limma_model.rds"))
  invisible(list(model = fit, results = result, design = design))
}

# Strict time units: numeric and Tn/Dn/Wn only. A mixed unit or labels like
# baseline/followup cannot be converted deterministically into a rate.
parse_longitudinal_time <- function(timepoint) {
  raw <- trimws(as.character(timepoint))
  if (anyNA(raw) || any(!nzchar(raw))) stop("Longitudinal timepoint is missing.")
  kind <- ifelse(grepl("^[+]?[0-9]+([.][0-9]+)?$",raw),"numeric",
    ifelse(grepl("^[Tt][+]?[0-9]+([.][0-9]+)?$",raw),"T",
      ifelse(grepl("^[Dd][+]?[0-9]+([.][0-9]+)?$",raw),"D",
        ifelse(grepl("^[Ww][+]?[0-9]+([.][0-9]+)?$",raw),"W","invalid"))))
  if (any(kind=="invalid") || length(unique(kind))!=1L)
    stop("Longitudinal timepoint must use one consistent numeric or Tn/Dn/Wn scale (e.g. T0,T1,T2); no mixed units or ambiguous visit names.")
  values <- suppressWarnings(as.numeric(if(kind[1]=="numeric")raw else substring(raw,2L)))
  if(any(!is.finite(values)))stop("Unparseable longitudinal timepoint.")
  list(time=values,timeScale=kind[1])
}

validate_longitudinal_reference_design <- function(meta,subject_column="subject_id") {
  if(!subject_column%in%names(meta)||!"condition"%in%names(meta)||!"time"%in%names(meta))
    stop("Longitudinal design needs subject_id, condition and numeric time.")
  ids <- as.character(meta[[subject_column]])
  conditions <- as.character(meta$condition)
  if (anyNA(ids)||any(!nzchar(ids))||anyNA(conditions)||any(!nzchar(conditions)))
    stop("Longitudinal subject_id/condition cannot be missing.")
  if(any(!is.finite(meta$time)))stop("Longitudinal time contains non-finite entries.")
  if(anyDuplicated(paste(ids,format(meta$time,digits=14),sep="\r")))
    stop("Duplicated subject-time observation. Aggregate technical replicates first; do not count technical measurements as biological subjects.")
  subject_conditions <- split(conditions,ids)
  if(any(lengths(lapply(subject_conditions,unique)) != 1L))
    stop("Condition changes within subjects: use a crossover/time-dependent design, not fixed-group longitudinal lmer.")
  by_subject <- split(meta$time,ids)
  if(length(by_subject)<10L || any(lengths(lapply(by_subject,unique))<2L))
    stop("Longitudinal reference requires >=10 independent subjects and >=2 observed distinct times for every subject; sparse subjects need a prespecified missing-visits model.")
  groups <- sort(unique(conditions))
  if(length(groups)!=2L)stop("This reference longitudinal contrast requires exactly two fixed treatment groups.")
  subjects_per_group <- table(vapply(subject_conditions,function(x)x[1],character(1)))
  if(any(subjects_per_group<5L))
    stop("Longitudinal reference needs at least 5 independent subjects per group.")
  for(group in groups)if(length(unique(meta$time[conditions==group]))<2L)
    stop("Each group needs at least two distinct time points.")
  list(nSubjects=length(by_subject),nObservations=nrow(meta),
       subjectsPerGroup=as.list(subjects_per_group),
       contrast=paste0("condition",groups[2],":time"),
       groups=groups,
       note="Longitudinal time slope assumes linear effect and stable group membership; checked eligibility does not certify MAR censoring, residual assumptions or adequate power.")
}

subject_has_repeated_time <- function(meta, subject_column = "subject_id", time_column = "time") {
  if (!all(c(subject_column,time_column)%in%names(meta)))return(FALSE)
  counts <- tapply(meta[[time_column]],meta[[subject_column]],
    function(x)length(unique(x[is.finite(x)])))
  counts <- counts[is.finite(counts)]
  length(counts)>=10L && sum(counts>=3L)>=10L
}

run_lmer_matrix <- function(
  matrix,
  metadata,
  output_dir,
  fixed_formula = "condition * time",
  subject_column = "subject_id",
  interaction_term = NULL,
  random_slope = "auto"
) {
  require_namespace("lmerTest")
  require_namespace("lme4")
  matrix <- as.matrix(matrix)
  storage.mode(matrix) <- "double"
  metadata <- as.data.frame(metadata)
  if(is.null(rownames(matrix))||is.null(rownames(metadata)))
    stop("Longitudinal sample identifiers must be explicit matrix and metadata row names.")
  if(anyDuplicated(rownames(matrix))||anyDuplicated(rownames(metadata)))
    stop("Longitudinal sample identifiers must be unique after technical replicate aggregation.")
  common <- intersect(rownames(matrix),rownames(metadata))
  if(length(common)<20L)stop("At least 20 matched longitudinal observations are needed.")
  x <- matrix[common,,drop=FALSE]
  meta <- metadata[common,,drop=FALSE]
  if(!"time"%in%names(meta)||!"condition"%in%names(meta))
    stop("Missing numeric time or group for explicit longitudinal contrast.")
  meta$time <- as.numeric(meta$time)
  design <- validate_longitudinal_reference_design(meta,subject_column)
  meta$condition <- factor(as.character(meta$condition),levels=design$groups)
  expected_term <- design$contrast
  if(!is.null(interaction_term) && !identical(interaction_term,expected_term))
    stop(paste("Invalid longitudinal contrast; expected",expected_term))
  if(!grepl("condition[[:space:]]*\\*[[:space:]]*time",fixed_formula))
    stop("Longitudinal reference requires an explicit condition-by-time interaction.")
  fixed <- stats::as.formula(paste0("~",fixed_formula))
  if(anyNA(stats::model.frame(fixed,data=meta,na.action=stats::na.pass)))
    stop("Missing covariates in a longitudinal model. Resolve before fitting.")
  design_matrix <- stats::model.matrix(fixed,data=meta)
  if(qr(design_matrix)$rank<ncol(design_matrix) ||
      !(expected_term%in%colnames(design_matrix)))
    stop("Aliased longitudinal fixed effects or missing explicit condition-by-time coefficient.")
  try_slope <- identical(random_slope,TRUE)||identical(random_slope,"yes")||
    (identical(random_slope,"auto") && subject_has_repeated_time(meta,subject_column,"time"))
  intercept_formula <- stats::as.formula(paste0("value ~ ",fixed_formula," + (1|",subject_column,")"))
  slope_formula <- stats::as.formula(paste0("value ~ ",fixed_formula," + (1 + time|",subject_column,")"))
  results <- vector("list",ncol(x))
  for(j in seq_len(ncol(x))) {
    feature <- colnames(x)[j]
    if(is.null(feature)||!nzchar(feature)) feature <- paste0("F",j)
    dat <- meta
    dat$value <- x[,j]
    dat <- dat[is.finite(dat$value),,drop=FALSE]
    row <- data.frame(
      feature=feature,term=expected_term,estimate=NA_real_,
      std_error=NA_real_,df=NA_real_,statistic=NA_real_,
      p_value=NA_real_,q_bh=NA_real_,
      n_observations=nrow(dat),
      n_subjects=length(unique(dat[[subject_column]])),
      random_effect_structure=NA_character_,
      status="not_estimable",stringsAsFactors=FALSE
    )
    if(nrow(dat)<20L) {row$status <- "insufficient_observations";results[[j]]<-row;next}
    eligible <- try(validate_longitudinal_reference_design(dat,subject_column),silent=TRUE)
    if(inherits(eligible,"try-error")) {
      row$status <- "incomplete_feature_specific_subject_visits"
      results[[j]]<-row;next
    }
    if(length(unique(dat$value))<3L||stats::sd(dat$value)<1e-10) {
      row$status <- "constant_feature";results[[j]]<-row;next
    }
    fd <- stats::model.matrix(fixed,data=dat)
    if(qr(fd)$rank<ncol(fd)) {
      row$status <- "feature_specific_rank_deficient"
      results[[j]]<-row;next
    }
    fit <- NULL
    used <- "random_intercept"
    fallback <- FALSE
    if(try_slope && subject_has_repeated_time(dat,subject_column,"time")) {
      candidate <- try(suppressWarnings(lmerTest::lmer(
        slope_formula,data=dat,REML=FALSE,
        control=lme4::lmerControl(optimizer="bobyqa"))),silent=TRUE)
      if(!inherits(candidate,"try-error")&&!lme4::isSingular(candidate,tol=1e-4) &&
          is.null(candidate@optinfo$conv$lme4$messages)) {
        fit <- candidate;used <- "random_intercept_and_slope"
      } else fallback <- TRUE
    }
    if(is.null(fit)) {
      candidate <- try(suppressWarnings(lmerTest::lmer(
        intercept_formula,data=dat,REML=FALSE,
        control=lme4::lmerControl(optimizer="bobyqa"))),silent=TRUE)
      if(!inherits(candidate,"try-error")&&!lme4::isSingular(candidate,tol=1e-4) &&
          is.null(candidate@optinfo$conv$lme4$messages))fit <- candidate
    }
    if(is.null(fit)) {
      row$status <- "singular_or_nonconvergent_fit"
      results[[j]]<-row;next
    }
    coefs <- summary(fit)$coefficients
    if(!(expected_term%in%rownames(coefs)) ||
        !all(c("Estimate","Std. Error","Pr(>|t|)")%in%colnames(coefs))) {
      row$status <- "missing_explicit_interaction"
      results[[j]]<-row;next
    }
    xterm <- coefs[expected_term,,drop=FALSE]
    beta <- unname(xterm[1,"Estimate"])
    se <- unname(xterm[1,"Std. Error"])
    p <- unname(xterm[1,"Pr(>|t|)"])
    if(!all(is.finite(c(beta,se,p))) || se<=0 || p<0 || p>1) {
      row$status <- "invalid_inference"
      results[[j]]<-row;next
    }
    row$estimate <- beta
    row$std_error <- se
    row$df <- if("df"%in%colnames(xterm))unname(xterm[1,"df"])else NA_real_
    row$statistic <- unname(xterm[1,"t value"])
    row$p_value <- p
    row$random_effect_structure <- used
    row$status <- if(fallback)"ok_random_slope_fallback" else "ok"
    results[[j]] <- row
  }
  result <- do.call(rbind,results)
  if(is.null(result))result <- data.frame()
  if(nrow(result)) {
    # BH denominator is the ORIGINAL supplied feature family: failures p=1.
    complete_p <- ifelse(is.finite(result$p_value),result$p_value,1)
    all_q <- stats::p.adjust(complete_p,method="BH")
    result$q_bh <- ifelse(is.finite(result$p_value),all_q,NA_real_)
  }
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(result,file.path(output_dir,"lmer_results.csv"),row.names=FALSE)
  attr(result,"design_summary") <- c(design,
    list(attemptedFeatures=ncol(x),
         estimableFeatures=sum(is.finite(result$p_value)),
         blockedFeatures=sum(!is.finite(result$p_value)),
         randomSlopeFallbacks=sum(result$status=="ok_random_slope_fallback")))
  result
}

run_fgsea_ranked <- function(
  ranks,
  pathways,
  output_dir,
  min_size = 10L,
  max_size = 500L,
  seed = 20260924L
) {
  require_namespace("fgsea")
  if (is.null(names(ranks))) stop("ranks must be a named numeric vector.", call. = FALSE)
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  set.seed(seed)
  result <- fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = as.integer(min_size),
    maxSize = as.integer(max_size)
  )
  result <- as.data.frame(result)
  if ("leadingEdge" %in% names(result)) {
    result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(output_dir, "fgsea_results.csv"), row.names = FALSE)
  invisible(result)
}


# Feature-wise R survival reference (not automated clinical confirmation).
# Unit of analysis is the unique independent subject, never each omics assay.
# This intentionally declines longitudinal/time-dependent and competing-risk
# designs; such models require their own protocol and validation.
run_cox_survival_reference <- function(
  matrix,
  metadata,
  output_dir,
  covariates = character(),
  ties = "efron",
  min_events = 12L
) {
  require_namespace("survival")
  if (!identical(ties,"efron") && !identical(ties,"breslow"))
    stop("Unsupported Cox ties policy.")
  x <- as.matrix(matrix)
  storage.mode(x) <- "double"
  if (is.null(rownames(x)) || !length(rownames(x)) || anyDuplicated(rownames(x)))
    stop("Survival reference requires unique sample identifiers on each matrix row.")
  m <- metadata[match(rownames(x), metadata$sample_id),,drop=FALSE]
  if (anyNA(m$sample_id) || anyNA(m$subject_id) ||
      anyDuplicated(m$subject_id) || any(!nzchar(as.character(m$subject_id))))
    stop("Survival reference requires exactly one observation per independent subject. Repeated visits need a time-dependent or clustered model.")
  if (nrow(x) < 20L)stop("Survival reference requires at least 20 independent subjects.")
  times <- suppressWarnings(as.numeric(m$survival_time))
  events <- suppressWarnings(as.numeric(m$survival_event))
  if (any(!is.finite(times) | times<=0) ||
      any(!is.finite(events) | !events%in%c(0,1)))
    stop("Complete, positive survival times and binary 0/1 event coding required; never silently omit or recode outcomes.")
  covariates <- unique(as.character(covariates))
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_]*$",covariates)) ||
      any(covariates%in%c("assay_feature","survival_time","survival_event","subject_id","sample_id")))
    stop("Invalid survival covariate identifier.")
  if (any(!covariates%in%names(m)))stop("Mapped survival covariate not present in sample metadata.")
  nuisance <- character()
  if ("batch"%in%names(m)) {
    batch_values <- as.character(m$batch)
    batch_values[is.na(batch_values)] <- ""
    if (length(unique(batch_values[nzchar(batch_values)])) > 1L) {
      if(any(!nzchar(batch_values)))stop("Batch covariate has missing values; resolve before Cox reference.")
      m$batch <- factor(batch_values)
      nuisance <- c(nuisance,"batch")
    }
  }
  for (name in covariates) {
    values <- m[[name]]
    if (anyNA(values) || any(!nzchar(as.character(values))))
      stop(paste("Missing values in survival covariate",name))
    numeric <- suppressWarnings(as.numeric(as.character(values)))
    if(all(is.finite(numeric)))m[[name]] <- numeric
    else m[[name]] <- factor(as.character(values))
    if(length(unique(m[[name]]))>1L)nuisance <- c(nuisance,name)
  }
  min_events <- max(12L,as.integer(min_events))
  min_events <- max(min_events,4L*(length(nuisance)+2L))
  if (sum(events) < min_events)stop(
    sprintf("Too few observed survival events for the declared model: %d; require >= %d. This is a minimum software guardrail, not evidence of adequate power or stable estimation.",sum(events),min_events)
  )
  m$survival_time <- times
  m$survival_event <- events
  output <- vector("list",ncol(x))
  form <- stats::reformulate(c("assay_feature",nuisance),
    response="survival::Surv(survival_time, survival_event)")
  for (j in seq_len(ncol(x))) {
    feature <- colnames(x)[j]
    if (is.null(feature) || !nzchar(feature))feature <- paste0("F",j)
    observed <- is.finite(x[,j])
    n_observed <- sum(observed)
    n_events <- sum(events[observed])
    result <- data.frame(
      feature=feature,effect=NA_real_,standardError=NA_real_,
      pValue=NA_real_,qValue=NA_real_,hazardRatio=NA_real_,
      ciLow=NA_real_,ciHigh=NA_real_,
      proportionalHazardsP=NA_real_,globalPhP=NA_real_,
      nSubjects=n_observed,nEvents=n_events,
      missingFraction=1-n_observed/nrow(x),
      status="not_estimable",stringsAsFactors=FALSE
    )
    if (n_observed < 20L || n_events < min_events) {
      result$status <- "insufficient_observations_or_events"
      output[[j]] <- result
      next
    }
    data <- m[observed,,drop=FALSE]
    data$assay_feature <- x[observed,j]
    if(stats::sd(data$assay_feature)<1e-10) {
      result$status <- "constant_feature"
      output[[j]] <- result
      next
    }
    design <- try(stats::model.matrix(stats::delete.response(stats::terms(form)),data=data),silent=TRUE)
    if (inherits(design,"try-error") || qr(design)$rank < ncol(design)) {
      result$status <- "aliased_design"
      output[[j]] <- result
      next
    }
    fit <- try(survival::coxph(
      form,data=data,ties=ties,x=TRUE,
      control=survival::coxph.control(iter.max=60,timefix=FALSE)
    ),silent=TRUE)
    if (inherits(fit,"try-error") || any(!is.finite(stats::coef(fit))) ||
        any(abs(stats::coef(fit))>12)) {
      result$status <- "cox_nonconvergence_or_separation"
      output[[j]] <- result
      next
    }
    beta <- as.numeric(stats::coef(fit)["assay_feature"])
    v <- try(stats::vcov(fit)["assay_feature","assay_feature"],silent=TRUE)
    if(inherits(v,"try-error") || !is.finite(v) || v<=0) {
      result$status <- "invalid_cox_variance"
      output[[j]] <- result
      next
    }
    se <- sqrt(v)
    zph <- try(survival::cox.zph(fit,transform="km"),silent=TRUE)
    if(inherits(zph,"try-error") || !"assay_feature"%in%rownames(zph$table)) {
      result$status <- "ph_diagnostic_unavailable"
      output[[j]] <- result
      next
    }
    feature_ph <- as.numeric(zph$table["assay_feature","p"])
    global_ph <- if("GLOBAL"%in%rownames(zph$table))
      as.numeric(zph$table["GLOBAL","p"]) else NA_real_
    if(!is.finite(feature_ph) || !is.finite(global_ph)) {
      result$status <- "ph_diagnostic_unavailable"
      output[[j]] <- result
      next
    }
    result$effect <- beta
    result$standardError <- se
    result$pValue <- 2*stats::pnorm(-abs(beta/se))
    result$hazardRatio <- exp(beta)
    result$ciLow <- beta-1.959963984540054*se
    result$ciHigh <- beta+1.959963984540054*se
    result$proportionalHazardsP <- feature_ph
    result$globalPhP <- global_ph
    result$status <- if(feature_ph<0.05 || global_ph<0.05) "ph_assumption_warning" else "ph_test_not_rejected"
    output[[j]] <- result
  }
  results <- do.call(rbind,output)
  # The full supplied family remains in the BH denominator: unsuccessful
  # fits receive p=1, not a post-hoc reduction in the number of tests.
  family_p <- ifelse(is.finite(results$pValue),results$pValue,1)
  q <- stats::p.adjust(family_p,method="BH")
  results$qValue[is.finite(results$pValue)] <- q[is.finite(results$pValue)]
  results <- results[order(is.na(results$qValue),results$qValue),,drop=FALSE]
  rownames(results) <- NULL
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(results,file.path(output_dir,"cox_results.csv"),row.names=FALSE)
  counts <- as.list(table(results$status))
  list(
    results=results,
    summary=list(
      method=paste0("survival::coxph, ties=",ties),
      attemptedFeatures=ncol(x),estimable=sum(is.finite(results$pValue)),
      nSubjects=nrow(x),events=sum(events),tieMethod=ties,
      nuisanceTerms=nuisance,statusCounts=counts,
      importantLimit="A non-significant cox.zph test never proves PH; informative censoring, competing risks, time-dependent exposure and selective feature missingness remain untested."
    )
  )
}
