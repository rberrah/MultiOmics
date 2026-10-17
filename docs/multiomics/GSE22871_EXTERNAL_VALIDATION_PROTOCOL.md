# External GEO qualification: GSE22871 (independent mice, heldout transfer NOT yet permitted)

## Registered independent cohort, not the training Nutrimouse study

The public National Center for Biotechnology Information GEO study **GSE22871** measures liver expression from **30 individual mice**, with wild-type versus PPARα-null genotype, exposed for seven days to 0, 3 or 10 mg/kg/day PFOS, **five biological animals per genotype × dose group**. The microarray is **Affymetrix GPL1261 Mouse Genome 430 2.0**. Read the actual original GEO sample title/GSM metadata; do not infer genotype/dose from sample order or chromosome patterns. Official record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22871.

The source archive is SHA-256 **`630706c5d3778068a9e08af01ce315b4886de17abefbfaae32708e7ed3d2d861`**, computed after an independent successful NCBI download in GitHub Actions [run 39016729326](https://github.com/rberrah/MultiOmics/actions/runs/39016729326). Subsequent downloads are **rejected** if any source byte changes, rather than silently interpreting a reprocessed version.

The original public **processed series matrix** is retrieved from `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22871/matrix/GSE22871_series_matrix.txt.gz` and archived with its observed SHA-256. It represents published **RMA, log2-transformed microarray probe intensities**, not raw CEL files or RNA-seq counts. The code does not masquerade as reprocessing the source CEL acquisitions. The source GitHub SHA and GEO downloaded-data SHA are separate provenance assertions; re-check if NCBI republishes data.

## Acceptance: real cohort intake (not validation of signature transfer)

1. Fetch the real GSE matrix and **reject** empty/broken/gzip-inconsistent, oversized or inaccessible sources. Preserve original compressed bytes and SHA-256 in an isolated CI artifact.
2. Match exactly 30 original GSM565106–GSM565135 IDs and matching matrix and GEO metadata order. Verify all original **30 independent animals** and the six genuine genotype × dose groups with 5 distinct within-group biological replicate identifiers.
3. Verify GPL1261 only, strict nonduplicated probe-set IDs, complete numeric finite values, expected approximate probe-row count, and sample-level observed intensity ranges. Report all values, original titles, units, groups and evidence without replacing failures with fabricated source data.
4. **Explicitly block** `Nutrimouse → GSE22871` model projection. Nutrimouse features are liver macroarray **gene symbols** among a preselected 120-feature panel; GSE22871 measurements are Affymetrix **probe IDs** on a different scale. Without a **versioned, independently checked probe-to-gene annotation**, mapping of multiple probes to one gene, and a frozen cross-platform scaling/normalization bridge, there is **no valid direct external score**.
5. No genotype-label-driven post hoc selection, no imputation claimed as raw, no assumed same physical sample across studies, no pooled pseudo-replicates, and **no patient-level PK data or concentrations**. The output declares `BLOCKED_UNTIL_ANNOTATION_AND_SCALE_BRIDGE` and zero exported per-sample signature scores by design.

## Independent official GPL1261 probe annotation audit

A dedicated pipeline fetches the original NCBI GEO platform annotation
[`GPL1261.annot.gz`](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1261/annot/GPL1261.annot.gz).
It is independently validated against the exact **45,101 original measured
probe IDs** in the SHA-256-locked GSE22871 source. The full annotation gzip
is preserved. Its qualified, immutable SHA-256 is
`394c103cf4ca99500d06369e4fb2b9b9bffc2c10fd3cbeb0b0d64c51c3fc95f8`,
independently observed from official NCBI bytes in
[Actions 39018594337](https://github.com/rberrah/MultiOmics/actions/runs/39018594337).
Future downloads **fail closed on digest drift**. The full observed digest is also
recorded in `GPL1261_annotation_mapping_audit.json`.

The genuine audit matched **45,101/45,101 source probes** against **45,101
official platform annotation rows**, yielded **22,587 distinct parsed gene-symbol
labels** and **1,196 probes with more than one gene symbol**. The prespecified
`Cyp3a11` anchor has **one candidate platform probe**, which is NOT by itself
evidence that the underlying expression measurements are on the same scale.

- Preserve every platform probe→gene-symbol pair from the *official* GEO
  annotation. Ambiguous probes mapping to multiple genes are retained with
  `multiple_genes` status; no arbitrary first-hit assignment.
- Preserve genes with multiple probe IDs and report how many source probes
  have missing annotation rows or missing/unparseable gene symbols.
  Separate machine-readable CSVs keep the ambiguous, missing and full
  mappings, rather than removing failures.
- Audit the prespecified original Nutrimouse `Cyp3a11` anchor only as
  **candidate platform symbols**. A same-named gene is neither equivalent
  to the original Nutrimouse measurement nor proof of enzyme activity.
- No gene/probe collapse is chosen by genotype/dose association, no
  average/maximum signal is computed without a declared bridge, and **zero
  external scores** are exported. The explicit status is
  `BLOCKED_SCALE_AND_PROBE_COLLAPSE_NOT_VALIDATED`.

The annotation source is an Affymetrix GPL mapping snapshot, not a clinical
genomic annotation reference; Gene Symbol aliases and updates can disagree
with modern Bioconductor `mouse4302.db`. An independent
`mouse4302.db` (version-pinned) concordance check and stable gene IDs
are an additional qualification gate before any publication-level
external projection. Annotation mapping is a necessary step, not proof of
biological transfer.

## Independent Bioconductor `mouse4302.db` second-reference audit

The latest release gate runs the **real Bioconductor** `mouse4302.db` 3.13.0
and `AnnotationDbi` with pinned R 4.5.2 / Bioconductor 3.22. It downloads
the **exact genuine GEO GSE22871 probe and sample files** from the preceding
job, rather than accepting a synthetic concordance table. The reference package
is separate from the SHA-256-pinned NCBI `GPL1261.annot.gz`.

For all 45,101 genuine platform probe IDs the independent R script extracts
the Bioconductor `SYMBOL` and `ENTREZID` associations and compares
**complete probe-level symbol sets**, not a hand-picked first hit. It reports
the exact, partly overlapping, disjoint, GEO-only, Bioconductor-only, and
both-unmapped categories. Every inconsistency and all Bioconductor Entrez IDs
are exported to machine-readable CSV without filtering on dose or genotype.
It separately tracks the prespecified mouse liver Cyp3a11 anchor.

**Important methodological limitation:** the matching criterion is literal
current gene-symbol-string agreement; a renamed alias does not imply a different
biological gene, and equal symbol strings do not establish equivalence of probe
sequence, signal scaling or platform preprocessing. Bioconductor Entrez IDs
are retained but **cannot be treated as concordant with original GEO gene IDs**
until those independently sourced IDs have also been compared. This independent
reference may legitimately disagree with the historical GPL annotation.
Consequently the CI gate requires that the **comparison actually runs and reports
all disagreements**, not that all probes agree.

The **actual independently reproduced comparison** on
[the genuine reference workflow](https://github.com/rberrah/MultiOmics/actions/runs/39019704728)
gave the following **probe counts**, from all 45,101 original GPL1261
measurements: **35,119** exact symbol-set agreements; **1,117** partial
overlaps; **2,069** disjoint symbol sets; **1,553** GEO-only assignments;
**291** Bioconductor-only assignments; and **4,952** unmapped on both sides.
These categories are mutually exclusive and sum to 45,101. Thus **9,982
probes (22.1%) do not show exact symbol-set agreement**. Their biological
meaning is not automatically disagreement at stable gene/Entrez level:
symbol renames and aliases may account for some differences. The frozen
Cyp3a11 anchor has **one candidate probe**; this is not evidence of
cross-platform numerical comparability.

No per-animal score, gene collapse, RNA expression normalization or phenotype
prediction is produced. The formal status stays
`BLOCKED_UNTIL_ANNOTATION_CONSENSUS_AND_SCALE_TRANSFER`. A future study must
predefine a representative-probe aggregation rule, train/test nonleaking
cross-platform measurement-scale transport and an external endpoint analysis,
all before using external group labels.

## Third reference audit — original NCBI Entrez Gene IDs, independently compared with Bioconductor

An additional Python verification reads the **actual `ENTREZ_GENE_ID` column
from the immutable NCBI GEO `GPL1261.annot.gz` archive**, without inferring
gene IDs from names or using gene-expression and animal group outcomes.
It compares every available numerical Gene ID against **genuine**
`AnnotationDbi::select(mouse4302.db, keytype="PROBEID", columns="ENTREZID")`
results produced by the separate pinned Bioconductor R job.

For **all 45,101 original platform probes**, the resulting report preserves
both sets of stable IDs, original raw GEO text, any unparseable entries and
the six cases *exact*, *partial*, *disjoint*, *GEO-only*, *Bioconductor-only*
and *both missing* (plus a seventh invalid-ID category). No "first ID" is
chosen when a probe maps to several genes. The independent audit also counts
cases where **Entrez IDs agree but symbol strings differ**. Such cases can
explain some name discrepancies without proving equivalent probe selectivity,
isoform measurement or cross-platform expression calibration.

Outputs: `original_NCBI_vs_Bioc_Entrez_ID/GPL1261_all_probe_entrez_concordance.csv`,
the full disagreements CSV, the CYP3A11 candidate subset, and the numerical
source-linked audit JSON. The original GEO/annotation SHA-256 and the
Bioconductor package version must match the prior validated evidence.
Failure of historical/current ID agreement is **recorded**, not "repaired"
through label-guided probe selection. Changes in retired/replaced Entrez IDs
may legitimately explain mismatches and require separate curator review.

The **real pinned-source results**, from
[GitHub Actions 39020750915](https://github.com/rberrah/MultiOmics/actions/runs/39020750915),
are **36,818** probe ID-set agreements, **1,121** partial overlaps, **367**
disjoint nonempty ID sets, **1,554** originally assigned only by NCBI,
**290** only by Bioconductor and **4,951** entirely unassigned. There were
**0 unparseable original Entrez entries**; all categories sum to the
original **45,101** measured probe IDs. Exactly **1,834 probes with
discordant gene symbol sets nevertheless have identical nonempty Entrez
sets**. This documents a real naming/reference-discrepancy subgroup;
it does **not** prove that each corresponding hybridization probe is
specific, or that the across-platform expression scales are comparable.

**Biological source article:** Rosen MB et al. (2010),
*Gene Expression Profiling in Wild-Type and PPARα-Null Mice Exposed to
Perfluorooctane Sulfonate Reveals PPARα-Independent Effects*,
*PPAR Research*, DOI
[`10.1155/2010/794739`](https://doi.org/10.1155/2010/794739).
The complete paper is openly accessible via
[PMC2948942](https://pmc.ncbi.nlm.nih.gov/articles/PMC2948942/);
no institutional credential or uploaded full text is currently required.

Even a perfect Entrez match is **not external validation of a numerical
biomarker**. The scientific status remains
`STABLE_ENTREZ_ID_COMPARISON_ONLY_EXTERNAL_PROJECTION_BLOCKED`:
no gene-level collapse, frozen cross-platform correction, PopPK prediction or
PK outcomes have been estimated. GSE22871 is a genuine new biological
cohort, but it is **not** a multiomics/PK cohort and cannot support clinical
pharmacokinetic claims.

## Remaining scientific gates for issue #20

- Freeze actual file SHA and annotation publication version and audit multi-probe ambiguity (one-to-many gene mapping). No manual selection of the probe with the strongest genotype correlation.
- Before looking at the independent outcomes, specify how to transport a train-only **cross-platform-invariant** statistic, including missing-gene policy and baseline pre-processing. The earlier Nutrimouse mean/SD are *not* directly interchangeable with GPL1261 RMA probe expression.
- Recompute external scores genuinely and independently in another implementation, preserve negative results, and add an eligible external biological endpoint (with treatment/dose factor stratification); compare matched cohorts without asserting that differing doses, genotypes, timepoints or assays are equivalent.
- No extrapolation to human CYP3A4/5, biological clinical biomarkers, tacrolimus PK/PD or MIPD. A technically qualified external dataset alone is not external validation of scores.

The other proposed GEO **GSE73298** must not be silently treated as 36 analyzed samples: its series text describes a 36-mouse design, but the accession currently lists **34 sample records**. Resolve the source discrepancy with true GSM evidence before intake or comparison.

## Reproducibility

```bash
python3 scripts/qualify_geo_gse22871.py --self-test
python3 scripts/qualify_geo_gse22871.py tmp/multiomics/real-external-GSE22871
python3 scripts/qualify_gpl1261_annotation.py --self-test
python3 scripts/qualify_gpl1261_annotation.py tmp/multiomics/real-external-GSE22871
```

All reproducibility and publication claims remain tied to the exact successful GitHub SHA and the original GEO gzip SHA. This work updates only MultiOmics, not the source/consumer repositories.
