# Design audit

`src/lib/multiomics/design-audit.js` — checks that run **before** a model is
chosen, and that explain why a result may not mean what it appears to mean.

Tested by `scripts/test_multiomics_design_audit.mjs`, wired into
`npm run test:multiomics`. Pure functions, no dependencies, bilingual
(`message` / `messageEn`) like the rest of `src/lib/multiomics`.

## Why

Every check here exists because it caught a real error in a real exposure
study. In that study:

* each treatment arm had **its own preparation control**, and mixing them
  reversed a published conclusion;
* one of those controls was **itself strongly bioactive** — it altered more
  features than any treatment did;
* the two preparations differed by a **2 mM cofactor** whose degradation
  products dominated the metabolome and were read as biology;
* the supplied protein table still contained **reverse-database decoys**,
  counted as results;
* zeros in a targeted panel meant **"not quantified"**, including for
  analytes present in every living cell;
* one layer had **no replication at all** yet still carried p-values and FDRs.

None of these are exotic. All of them are detectable from the metadata and the
matrix alone.

## Severity contract

Each finding is `{ code, severity, message, messageEn, detail }` with
`severity` in `error` | `warning` | `info`. `auditDesign()` sorts errors first
and returns `ok: false` with `blocking: string[]` when any error is present.
Nothing here throws on absent input: a check whose inputs are missing is
skipped silently, so a partial call is safe.

## API

### `buildContrastPlan({ arms, controlOf, untreated })`

Turns an arm → control mapping into the full contrast set:

| Kind | Meaning |
|---|---|
| `treatment_vs_own_control` | each arm against the control it was declared with |
| `control_vs_control` | every pair of controls — the contrast nobody reports |
| `control_vs_untreated` | each control against untreated cells, when present |
| `head_to_head_shared_control` | plain difference, valid only when the control is the same |
| `head_to_head_corrected` | `(A − ctrlA) − (B − ctrlB)` when the controls differ |

An arm with no declared control raises `arm_without_control` (error) rather
than silently falling back to a global reference.

```js
const plan = buildContrastPlan({
  arms: ['Untreated','Vehicle','Blank','DrugLow','DrugHigh','PoolLow','PoolHigh'],
  controlOf: { DrugLow:'Vehicle', DrugHigh:'Vehicle',
               PoolLow:'Blank',   PoolHigh:'Blank' },
  untreated: 'Untreated'
});
```

### `assessBaselineMismatch(plan, significantByContrast)`

Given a count of significant features per contrast id, raises
`baseline_mismatch` (error) when a control-versus-control contrast is larger
than the biggest treatment effect it is supposed to anchor. Counts may be
keyed in either direction (`A__vs__B` or `B__vs__A`).

This is the single highest-value check: it turns a silent interpretation error
into one sentence.

### `assessReplication(metadataRows)`

Counts **independent units** per `condition × timepoint` cell, from
`subject_id`. Repeated assays of the same subject do not count.

* `no_replication` (error) — at least one cell has a single unit. No
  within-group variance can be estimated, so any p-value reflects an assumed
  dispersion. Use `reproducibilityAcrossBlocks()` instead.
* `minimal_replication` (warning) — two units somewhere.
* `replicated` — three or more everywhere.

### `variancePartition(matrix, factors)`

Mean share of per-feature variance explained by each declared factor
(marginal, one-way). Raises `nuisance_factor_dominant` when a factor named
like a batch/replicate/series/run/plate/experiment carries variance comparable
to the largest factor — i.e. when an unblocked test would be biased and
underpowered.

Fractions can sum above 1 when factors are confounded; that is reported as
`factors_confounded` rather than hidden, because a confounded nuisance factor
is exactly what makes a design unanalysable.

### `assessPValueCalibration(pValues)`

A calibrated null puts about 5 % of tests below 0.05.

* `conservative` — under 2.5 %. Typically a hand-set dispersion (edgeR/DESeq2
  without replicates). The displayed FDR is not calibrated.
* `inflated` — over 25 % with no flat tail. Usually an unmodelled nuisance
  factor or non-independent samples.
* `too_few_tests` — fewer than 200 tests; no verdict.

### `auditFeatureTable(featureIds, matrix, options)`

* `decoy_features_present` (error) — `Reverse_`, `DECOY`, `rev|` identifiers.
* `contaminant_features_present` (warning) — `CON__`.
* `zeros_are_missing` (error) — an analyte present in every living cell (ATP,
  NAD, glutathione, glutamate, lysine, citrate, taurine, …) reported as zero
  in **every** sample. The zeros encode "not quantified".
* `all_zero_features` (warning) — otherwise, when features are zero
  everywhere. Suppress with `{ zerosAreMissing: true }` once confirmed.

### `detectNonBiologicalSignals({ featureIds, matrix, groupOf, testArm, referenceArm, options })`

Flags features whose behaviour between two arms is chemically rather than
biologically plausible:

* `implausible_fold_change` — beyond `maxAbsLog2FC` (default 5, i.e. 32×);
* `present_in_one_arm_only` / `absent_in_one_arm_only` — detected in ≥ 80 % of
  one arm and ≤ 20 % of the other.

Flagged features are then matched against **formulation families** — NADPH and
its degradation route, microsome/membrane buffer, solvents and vehicles,
medium supplements. Two or more members of one family moving together raises
`formulation_family_detected` (error): the difference between those arms
includes a difference in composition, not only a biological effect.

Present these in a *"check the formulation"* panel, not in the results.

### `reproducibilityAcrossBlocks({ featureIds, effectByBlock, options })`

When a design has no within-cell replication but repeats across blocks (time
points, sites, runs), those blocks are de facto independent experiments. A real
effect reproduces across them; noise does not.

Returns pairwise Pearson correlations of the effect-size vectors, the median,
and the features exceeding `minAbsEffect` (default 1) **with the same sign in
every block**. Raises `effects_not_reproducible` (error) when the median
correlation is below 0.1.

This is the substitute view to offer whenever `assessReplication()` returns
`no_replication`: it converts an otherwise unusable upload into a defensible
analysis.

### `auditDesign(input)`

Runs every applicable check and returns
`{ ok, blocking: string[], findings, plan }`, errors first. Every input key is
optional.

```js
const audit = auditDesign({
  armPlan, significantByContrast, metadata,
  featureIds, matrix, factors, pValues,
  nonBiological: { featureIds, matrix, groupOf, testArm, referenceArm },
  reproducibility: { featureIds, effectByBlock }
});
if (!audit.ok) render(audit.findings);
```

## Relation to the existing MS blank filter

The metabolomics QC in `deterministic-core.js` filters **analytical** blanks —
solvent, extraction and process injections, annotated through
`sample_type = blank` — by comparing biological and blank medians.

`design-audit.js` addresses a different problem: a **preparation control that
is itself an exposure arm**, such as a microsomal incubation blank. That arm
must stay `sample_type = biological`, or the MS blank filter would subtract
away the very contrast that matters. The two mechanisms are complementary and
must not be conflated.

## Relation to the confirmatory display gate

`confirmatory-result-gate.js` is a fail-closed guard on **how a result was
computed** (native R tables for every layer). `design-audit.js` is a guard on
**whether the design can answer the question at all**. A study can pass the
confirmatory gate and still fail the design audit — the real study that
motivated this module would have.

## Suggested integration

1. Run `auditDesign()` on upload, before any model is offered, and render the
   findings above the analysis controls.
2. Block the confirmatory intent while any `error` finding is unresolved.
3. When `no_replication` fires, offer the reproducibility view instead of
   p-value tables.
4. Carry `blocking` into `DOWNSTREAM_EXCHANGE_V1` so Omics2PK can see that a
   matrix came from a design with a known defect.
