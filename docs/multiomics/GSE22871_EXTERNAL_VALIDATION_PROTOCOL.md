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
