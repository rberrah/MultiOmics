# External public LC-MS qualification

## Real measured peak-area tables: independent Drosophila experiments

Upstream dataset: https://zenodo.org/records/21461823
Upstream code: https://github.com/hubertdl/Metabolomics_sampling_project_D.melanogaster
Frozen upstream commit: `7c1649ab8093e2175e1683d4ee88af6ab1aa4fea`.

CI downloads actual Experiment 1 and 2 raw peak-area CSVs and their
accompanying metadata at the exact commit. It checks each downloaded
Git blob hash, then runs source-level quantitative MS quality assurance
separately by chromatography column and ionization polarity.
These are **independent measured biological experiments but the same research
group**, not two independent laboratories and not raw mzML acquisitions.

Experiment 1 has 35 biological samples and 6 pooled-QC profiles;
Experiment 2 has 48 biological samples and 3 pooled-QC profiles. Neither
provides corresponding blank injections, and the full QC injection order
is undocumented. Biological sample population, strain, diet and pool size
must not be relabelled as analytical batches. The engine therefore records
blank contamination, QC drift and between-instrument-batch comparisons as
**not evaluable**. QC variance is calculated only from actual measured
replicate QC intensities. Feature IDs and biological metadata must match
exactly; no imputation or automatic outlier removal takes place.

## Independent laboratory follow-up: Sacurine / MTBLS404

Published source: Thévenot et al. (2015), J Proteome Res,
DOI:10.1021/acs.jproteome.5b00354.
https://workflow4metabolomics.org/lcms_dataset-sacurine
https://workflow4metabolomics.org/W4M00002

The public 36-file LC-MS mzXML subset (about 188 MB) includes 20
biological, 10 pooled-QC and 6 blank samples. The full data covers two
analytical batches, but two-batch coverage of the subset must be checked
from original injection metadata, not inferred.

Prespecified next scientific gate: obtain exact files, verify source
checksums and sample/role metadata, run frozen genuine `xcms` processing,
measure pooled-QC CV and actual blank/biological signal ratios, record
observed injection order only when documented, and compare quantitative
diagnostics with Workflow4Metabolomics. Preserve all missing or
non-evaluable cases. Do not tune peak-picking against biological effect
size or declare metabolite identification from peak areas.

These external files are NOT bundled into MultiOmics Git source, and
the Sacurine raw-data comparison must NOT be described as finished
before a complete green actual-data run.

## Scientific article claims — evidence matrix

| Claim | Evidence | Present status |
|---|---|---|
| Standalone beginner tool | Browser and R contracts, Playwright | Qualified |
| Reproducible native omics methods | Airway, Nutrimouse, frozen R/Python | Qualified within tested designs |
| Native raw MS data extraction | Authentic faahKO CDF subset | Qualified within four files |
| Planted QC failures | Synthetic truth-known blanks, drift and batches | Qualified on synthetic data |
| QC on unseen real measured profiles | Two Drosophila experiments | Under independent CI qualification |
| Real blanks, QC drift, cross-lab portability | Sacurine and additional external labs | Not yet qualified |
| Universal FDR or biomarker confirmation | Larger prospective multi-lab benchmarks needed | Not established |
| User-friendliness / UI visual parity | External novices, expert reviews, FR/EN/theme mobile | In progress |

A current article could describe architecture and qualification methodology.
A strong validation paper should report the external Sacurine and comparison
results before making general scientific-performance claims.
