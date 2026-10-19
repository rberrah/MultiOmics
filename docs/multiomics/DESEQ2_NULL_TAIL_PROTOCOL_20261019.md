# DESeq2 global-null tail diagnostic — prespecified protocol

**Protocol frozen before reading the new feature-level results: 19 October 2026.**

This experiment follows issue #44. It does **not** change DESeq2, MultiOmics production p-values, filtering or user-facing thresholds.

## Question

Previous whole-study null experiments showed more DESeq2 studies with at least one BH q≤0.05 false call than expected under the tested generators. True simulated dispersion often moved the Wald second moment closer to 1 and reduced false calls in some generators. The next question is narrower:

> Is the excess accompanied by systematically heavy native Wald/p-value tails within whole null studies, and is that tail behaviour reduced when the exact simulated NB dispersion is supplied?

## Design fixed before calculation

- Reuse the independently coded frozen NB generator in `rnaseq_deseq2_null_oracle.R`.
- Three source scenarios: `balanced_base`, `balanced_depth_dispersion`, `unequal_partial_batch`.
- **32 new independent studies per scenario**, replicate IDs **6001–6032**.
- Same 400-gene/40-subject source design, `~ batch + condition`, and the same count prefilter as the prior oracle studies.
- Two paired native DESeq2 Wald arms on every exact same count matrix:
  1. `native_default`: normal DESeq2 dispersion/size-factor estimation.
  2. `oracle_true_dispersion`: estimated DESeq2 size factors but exact generator NB dispersion supplied before `nbinomWaldTest`.
- No production inference setting is changed. For this diagnostic only, Cook filtering and independent filtering are disabled at `results()` so raw Wald/p-value tails are not conflated with screening.

## Prespecified tail summaries

For every **whole study** and arm:

- fraction of finite p-values ≤ 0.10, 0.05, 0.01, 0.001;
- fraction of finite |Wald Z| ≥ 1.644854, 1.959964, 2.575829, 3.290527;
- mean Z and mean Z²;
- 0.90, 0.95, 0.99 and 0.999 quantiles of |Z|;
- number of BH q≤0.05 calls.

Genes are **not** treated as independent biological replicates and no pooled-gene significance test will be used. Inference is summarized across the 32 independent studies in each scenario.

For each scenario/arm, report mean, median and IQR of study-level tail metrics. For the paired arm comparison, report the distribution of per-study differences `native_default - oracle_true_dispersion`.

## Interpretation boundary

A heavier default tail that shrinks under true dispersion would strengthen the hypothesis that dispersion estimation contributes to the observed simulated-null behaviour. It would not prove a universal DESeq2 defect, because the generator, finite-sample Wald approximation, dependence among genes and study design still matter.

No production method may be selected, modified or replaced based on whichever arm looks favorable. All results, including negative ones, are retained.
