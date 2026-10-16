#!/usr/bin/env Rscript
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/validation_benchmarks.R")

# Exactly known data-generating structure, two independent omic views.
data <- simulate_mofa_known_truth(seed=20261009L,
  scenario="biological_plus_batch")
mofa_input <- mofa2_orient_reference_matrices(data$blocks)
stopifnot(
  nrow(mofa_input$transcriptomics)==32L,
  ncol(mofa_input$transcriptomics)==90L,
  identical(colnames(mofa_input$transcriptomics),data$truth$sample),
  identical(colnames(mofa_input$proteomics),data$truth$sample),
  identical(rownames(mofa_input$transcriptomics),
            colnames(data$blocks$transcriptomics)),
  identical(mofa_input$transcriptomics[1,1],
            data$blocks$transcriptomics[1,1])
)
stopifnot(identical(data$scenario,"biological_plus_batch"),
          nrow(data$truth)==90L,
          setequal(rownames(data$blocks$transcriptomics),data$truth$sample),
          setequal(rownames(data$blocks$proteomics),data$truth$sample),
          abs(stats::cor(data$truth$biological,data$truth$technical)) < 0.40)

# Factor order and sign MUST NOT affect truth recovery scores.
df <- rbind(
  data.frame(sample=data$truth$sample,factor="Factor1",
             value=-data$truth$technical),
  data.frame(sample=data$truth$sample,factor="Factor2",
             value=data$truth$biological)
)
audit <- mofa_factor_truth_diagnostics(df,data$truth)
stopifnot(
  audit$maxBiological>0.9999,
  audit$maxTechnical>0.9999,
  audit$biologicalRecovered,
  audit$batchSignalDetected,
  audit$distinctTwoFactorRecovery,
  length(unique(audit$matchedIndependentFactorNames))==2L,
  !audit$batchDominated,
  identical(audit$strongestBiologicalFactor,"Factor2"),
  identical(audit$strongestTechnicalFactor,"Factor1")
)
swapped <- df
swapped$factor <- ifelse(df$factor=="Factor1","Factor2","Factor1")
swapped$value <- -swapped$value
again <- mofa_factor_truth_diagnostics(swapped,data$truth)
stopifnot(
  isTRUE(all.equal(audit$maxBiological,again$maxBiological)),
  isTRUE(all.equal(audit$maxTechnical,again$maxTechnical))
)

# One mixed factor can correlate with BOTH truths but must not count as
# recovery of two distinct biological/technical axes.
mixed_factor <- data.frame(sample=data$truth$sample,factor="OneMixed",
  value=data$truth$biological+data$truth$technical)
mix_check <- mofa_factor_truth_diagnostics(mixed_factor,data$truth,
  biological_threshold=0.55,batch_threshold=0.55)
stopifnot(mix_check$biologicalRecovered,mix_check$batchSignalDetected,
  !mix_check$distinctTwoFactorRecovery)

# A factor that only reflects batch must be WARNED, not sold as biology.
batch_only <- simulate_mofa_known_truth(seed=20261011L,scenario="batch_only")
technical <- data.frame(sample=batch_only$truth$sample,factor="Factor1",
                        value=batch_only$truth$technical)
tech_audit <- mofa_factor_truth_diagnostics(technical,batch_only$truth)
stopifnot(tech_audit$batchDominated,
          tech_audit$batchSignalDetected,
          !tech_audit$biologicalRecovered)
bad <- df[-1,,drop=FALSE]
stopifnot(inherits(try(mofa_factor_truth_diagnostics(bad,data$truth),
                        silent=TRUE),"try-error"))
bad <- df
bad$sample[2] <- bad$sample[1]
stopifnot(inherits(try(mofa_factor_truth_diagnostics(bad,data$truth),
                        silent=TRUE),"try-error"))
bad <- df
bad$value[1] <- Inf
stopifnot(inherits(try(mofa_factor_truth_diagnostics(bad,data$truth),
                        silent=TRUE),"try-error"))

# Both data generators must produce plausible batch/non-batch partitions.
clean <- simulate_mofa_known_truth(seed=20261010L,scenario="biological")
stopifnot(!clean$expectedTechnical,clean$expectedBiological)
stopifnot(!batch_only$expectedBiological,batch_only$expectedTechnical)
cat("MOFA truth-labelled synthetic simulator, sign/order invariant factor scoring, batch-only negative control and refusal of missing/duplicated/nonfinite scores: PASS\n")
