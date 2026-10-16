# Raw-MS validation: real faahKO, frozen xcms, and synthetic QC stress

Status: native scientific qualification. No frontend or PMx Explain changes.

## Native R package environment (measured, not guessed)

The real publicly available `faahKO` chromatographic run passed on R **4.5.2**
and Bioconductor **3.22** using the immutable official Docker image already
used by MultiOmics.

Versions actually recorded in
[run 38913534547](https://github.com/rberrah/MultiOmics/actions/runs/38913534547),
artifact `12051317608`, include:

| Package | Observed |
|---|---|
| xcms | 4.8.0 |
| MsExperiment | 1.12.0 |
| Spectra | 1.20.1 |
| mzR | 2.44.0 |
| faahKO | 1.50.0 |
| BiocParallel | 1.44.0 |
| BiocManager | 1.30.27 |
| jsonlite | 2.0.0 |

The source-controlled lock under
`multiomics-engine/environment/raw-ms-reference-runtime.json` also records
other measured direct packages and the SHA-256 fingerprint of **all 206**
installed R packages. The verification script checks exactly the recorded
direct versions, transitive package count, R/Bioconductor release and
206-library digest, then records another inventory.

**Important limitation:** exact-version validation is NOT equivalent to a
self-contained offline installation, a frozen package-repository mirror, or
bitwise-identical chromatographic fits on different platforms. Since
CRAN/Bioconductor/apt indexes may move, mismatched versions cause a **hard
failure**, not an automatic upgrade or an undocumented alternate result.

## Real chromatography integration

The authentic public `faahKO` subset contains **two KO and two WT**
open-format CDF files. A previously qualified run yielded 33,016 individual
chromatographic peaks and 12,358 grouped features for this configuration.
The pipeline checks unaltered R native outputs, sample names and identities,
lossless feature-definition RDS and JSON-encoded list-column CSVs.

These samples do **not** contain dedicated pooled-QC or solvent-blank
injections. In the new release, the *real* faahKO output consequently
reports the impossibility of evaluating blank contamination, pooled-QC
precision, run-order drift, and replicated pooled-QC batch shift.
An honest missing-coverage alert is scientifically stronger than inventing
QC results from biological samples.

## Independent synthetic falsification fixtures

`test_ms_qc_multiscenario_stress.R` constructs reproducible, **synthetic**
intensity data with known conditions. It checks four independent scenarios:

1. Full, orthogonally controlled peak areas: a blank contaminant
   (ratio 0.4), monotonically drifting pooled QC (Spearman rho = 1),
   batch-dependent QC intensities (median ratio 1.6), group-specific
   missingness, plus an analyte that is completely missing in pooled QC.
2. No pooled QC or blank injections: the method must state which screens
   are unobservable and never certify these absent measurements.
3. Perfect confounding of biological condition with MS analytical batch:
   the summary explicitly warns the covariate effect cannot be separated.
4. Partial per-feature pooled-QC observations and invalid thresholds:
   unavailable QC measurements must remain flagged, never imputed or
   relabelled as a successful screen.

The native R backend now supplies **descriptive** flags for absolute
pooled-QC injection-order Spearman correlation and the maximum/minimum
median pooled-QC intensity across instrument batches. The predefined
heuristic defaults (rho ≥ 0.80 and median ratio ≥ 1.25) are quality
**screens**, not hypothesis tests, p-values, automatic batch correction,
or accepted instrument-agnostic clinical thresholds.

A true multi-laboratory confirmatory evaluation would require multiple
external studies with actual blanks and pooled QC injections, varied
instruments and acquisition protocols, drift correction
prespecification, independent compound identification and suitable
statistical methods. None are claimed here.

The beginner interface, public site URLs, and browser-only workflow remain
unchanged; these tests execute solely in the standalone source CI.
