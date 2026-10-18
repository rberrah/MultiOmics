# Reuse *only* frozen complete-pair Monte Carlo metrics, not the original
# count generator. Explicit scenario mapping is internal to aggregation;
# scientific outputs retain the original Gamma–Poisson scenario names.
source("multiomics-engine/rnaseq_deseq2_dispersion_component_metrics.R")
gp_dispersion_independent_summary <- function(dat) {
  gp<-c("gp_balanced_heterogeneous","gp_depth_heterogeneous",
        "gp_unbalanced_batch")
  original<-c("balanced_base","balanced_depth_dispersion",
              "unequal_partial_batch")
  if(!is.data.frame(dat)||!nrow(dat)||!"scenario"%in%names(dat)||
     !setequal(as.character(unique(dat$scenario)),gp) ||
     anyNA(dat$scenario) ||
     !all(dat$scenario%in%gp))
    stop("Independent null summary must use all three original Gamma–Poisson scenarios.")
  mapped<-dat
  mapped$scenario<-original[match(dat$scenario,gp)]
  if(anyNA(mapped$scenario))
    stop("A Gamma–Poisson experiment was mislabeled.")
  # Study-level uncertainty, all 12 paired cells, no gene pseudoreplication,
  # same stable statistical summary used for the original generator.
  result<-rnaseq_deseq2_dispersion_component_summary(
    mapped,per_scenario=32L,replicate_start=5001L)
  result$scenario<-gp[match(result$scenario,original)]
  if(anyNA(result$scenario)||nrow(result)!=12L)
    stop("Incomplete independent source result.")
  result$source_gamma_poisson_independently_sampled<-TRUE
  result$generator_includes_off_trend_per_gene_alpha<-TRUE
  result$generic_research_FDR_claim_permitted<-FALSE
  result
}
