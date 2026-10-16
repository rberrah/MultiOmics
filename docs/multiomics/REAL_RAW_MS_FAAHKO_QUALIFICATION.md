# Real open-format LC-MS method qualification — faahKO / xcms

Scope: independent native R software qualification. The beginner-facing web
interface and its CSV peak-area/AUC import are not modified.

The existing `multiomics-engine/tests/raw_ms_xcms_integration.R` was a
potential real-data test but was not part of the standalone qualification
workflows after the repository extraction. This document describes the new
independent `.github/workflows/native-raw-ms-faahko.yml` job.

## Real public evidence

- Source: actual four raw **netCDF (.CDF)** chromatography/mass-spectrometry
  files shipped with the Bioconductor `faahKO` research example package.
- Prespecified subsets: two KO biological files and two WT biological files.
- Genuine native R packages: `xcms`, `MsExperiment`, `Spectra`, `mzR`,
  `faahKO`, on R 4.5.2 and Bioconductor 3.22 in the immutable image already
  used for the other MultiOmics reference workflows.
- Complete processing: CentWave detection, optional Obiwarp retention
  alignment, PeakDensity correspondence, missing-peak filling and
  chromatographic peak-area export.
- The chosen algorithms and parameters are **an integration fixture, not
  a universally optimized/confirmed instrument-specific analysis method**.

## Assertions and preserved outputs

1. The genuine open raw files load and produce nonempty peaks/features.
2. Returned features×samples matrix has exactly four original study samples.
3. The saved `metabolomics_peak_area.csv` matches the numeric R matrix, not
   merely an expected filename.
4. Feature IDs are unique and agree with the separately exported definitions.
5. Sample metadata and KO/WT membership match the original input manifest.
6. Nonfinite/missing intensity accounting agrees with the model summary.
7. Vendor-native formats are refused until converted; invalid peak-picking
   parameters and unsupported formats are rejected in package-free tests.
8. All outputs, including negative findings, package versions and source SHA,
   are preserved as a downloadable Actions artifact.

The `into` values here are **integrated chromatographic peak areas**, NOT a
pharmacokinetic concentration–time AUC. A browser CSV field called "MS peak
AUC" refers to analytical peak integration, not PK systemic exposure.

## Interpretation boundaries

Four biological raw files are enough for an API/data-integrity regression
check but **not** a population FDR calibration, external omics validation,
an instrument transfer study or meaningful claim of differential abundance.
No significant finding, biological effect, or model superiority is declared
by this CI job.

Future extensions can independently test instrument-batch shifts, blank
carryover, pooled QC drift, identification certainty, missing-not-at-random
mechanisms and different chromatographic acquisition modes. They must not
be silently treated as solved by this single faahKO test.

This job does not add controls or expert settings to the beginner UI.
