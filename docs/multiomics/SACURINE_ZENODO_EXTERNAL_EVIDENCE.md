# External Sacurine raw LC-MS replay — Zenodo 1346742

The larger W4M hosted MTBLS404 36-file archive is unavailable from
GitHub Actions (provider timeouts); the attempted default EBI
ISA paths returned 404. A **separate public Zenodo research record**
provides a small, directly addressable independent raw-data fixture:

- Source: https://zenodo.org/records/1346742
- DOI: 10.5281/zenodo.1346742
- Original published Sacurine research: Thévenot et al.,
  J Proteome Res (2015), DOI 10.1021/acs.jproteome.5b00354
- Raw sources: `HU_neg_099.mzXML` and `HU_neg_178.mzXML`
  (two human urine biological acquisitions); `Blanc12.mzXML`
  and `Blanc17.mzXML` (two authentic solvent-blank acquisitions).
- Public original md5 hashes verified against Zenodo metadata *before*
  the native analysis; SHA-256 is then saved in an immutable provenance
  manifest for each genuine downloaded file.
- Frozen method execution: original Bioconductor R 4.5.2 / Bioc 3.22
  container, original 206-package R fingerprint, prespecified
  `xcms` centWave peak finding / Obiwarp alignment / density
  correspondence / gap filling / `into` chromatographic area.
- The 2 blanks are explicitly excluded from *defining* chromatographic
  feature groups but remain in the exported intensity matrix.
- Per-feature genuine blank-to-urine median intensity ratios are
  recorded only when both original blanks and both original urine
  profiles have finite numeric measurements and the urine median is
  positive; otherwise the metric is marked **not evaluable**.
- The exported CSV numerical matrix, raw `xcms` feature groups,
  original sample identifiers and observed roles are reconciled.
- Original mzXML files and research subject-level data are **not
  added to Git** or uploaded as public GitHub Actions artifacts. Only
  method/package fingerprints, source checksums, aggregate findings
  and anonymous feature-level ratios are preserved.

## Interpretation and limits

This is a genuine independently published human LC-MS acquisition;
it is not fabricated experimental data. **Two biological and two blank
samples cannot validate the population, group-comparison error rates,
biological significance, assay sensitivity, temporal drift,
or pooled-QC precision.** No pooled QC exists in this toy subset,
and instrument injection chronology is not independently verified.

The criterion of this CI is **software/external raw-data compatibility
plus measured blank-background descriptives**, NOT an assertion that a
particular feature is definitively contaminated or identified.
An analytical blank ratio above a prespecified descriptive value
is reported as a screen, not a calibrated false-positive rate.

The later full 36-file MTBLS404 sample or a complementary
`zenodo.org/records/3244991` training subset (six biological
and three pooled QC) will be needed for genuine QC precision. Do not
combine blank measurements from unrelated acquisition subsets
without explicit run metadata and technical compatibility checks.

See also `SACURINE_EXTERNAL_RAW_QUALIFICATION.md`.
