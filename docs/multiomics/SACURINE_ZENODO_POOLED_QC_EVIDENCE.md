# Real Sacurine pooled-QC LC-MS source: 9 published Galaxy tutorial mzML

The [Galaxy Training Network LC-MS tutorial](https://training.galaxyproject.org/training-material/topics/metabolomics/tutorials/lcms/tutorial.html)
cites the individual original mzML acquisitions from
[Zenodo record 3244991](https://zenodo.org/records/3244991).
We use exactly **six named urine runs** (`HU_neg_048/090/123/157/173/192`)
and **three named pooled QC runs** (`QC1_002/008/014`).

This is genuine experimental raw LC-MS, not simulated values or the
previous four-file Zenodo 1346742 toy subset.

The test first downloads original Zenodo metadata, confirms the
nine names belong to the specified immutable record, validates the
**original published per-file MD5**, and captures independent per-file
SHA-256. No patient biological labels, sample order, or diagnostic
information are inferred beyond the roles stated by the source.

The reference `xcms` environment is the previously verified R 4.5.2 /
Bioconductor 3.22 / **206-package** release. The predeclared
centWave/Obiwarp/PeakDensity/gap-filling workflow runs on original
spectra, then checks native peak definitions and per-sample exports
against the numeric R matrix. Three measured pooled-QC profiles permit
descriptive CV/RSD where at least three intensities are finite.

**Not measured:** blanks, technical batches or full QC acquisition
chronology. No contamination ratio, run-order drift,
independent biological hypothesis test, clinically validated
biomarker or controlled false discovery rate can be concluded from
this small tutorial subset.

The companion `SACURINE_ZENODO_EXTERNAL_EVIDENCE.md` uses a
different original Zenodo record to examine two *real* urine and two
*real* blank acquisitions. We do **not merge their intensity matrices**,
since instrument run compatibility and common injection order have
not been independently demonstrated. This separation avoids false
confidence from splicing biological/blank/QC observations across
distinct subsets.

Both experiments are conditional on real GitHub Actions and archived
numeric evidence; a source download failure is an availability issue,
not proof of invalid QC. No original raw human spectra are committed
to GitHub or published as CI artifacts.
