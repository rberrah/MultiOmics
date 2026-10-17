# External GEO qualification: GSE22871 (independent mice, heldout transfer NOT yet permitted)

## Registered independent cohort, not the training Nutrimouse study

The public National Center for Biotechnology Information GEO study **GSE22871** measures liver expression from **30 individual mice**, with wild-type versus PPARα-null genotype, exposed for seven days to 0, 3 or 10 mg/kg/day PFOS, **five biological animals per genotype × dose group**. The microarray is **Affymetrix GPL1261 Mouse Genome 430 2.0**. Read the actual original GEO sample title/GSM metadata; do not infer genotype/dose from sample order or chromosome patterns. Official record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22871.

The original public **processed series matrix** is retrieved from `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22871/matrix/GSE22871_series_matrix.txt.gz` and archived with its observed SHA-256. It represents published **RMA, log2-transformed microarray probe intensities**, not raw CEL files or RNA-seq counts. The code does not masquerade as reprocessing the source CEL acquisitions. The source GitHub SHA and GEO downloaded-data SHA are separate provenance assertions; re-check if NCBI republishes data.

## Acceptance: real cohort intake (not validation of signature transfer)

1. Fetch the real GSE matrix and **reject** empty/broken/gzip-inconsistent, oversized or inaccessible sources. Preserve original compressed bytes and SHA-256 in an isolated CI artifact.
2. Match exactly 30 original GSM565106–GSM565135 IDs and matching matrix and GEO metadata order. Verify all original **30 independent animals** and the six genuine genotype × dose groups with 5 distinct within-group biological replicate identifiers.
3. Verify GPL1261 only, strict nonduplicated probe-set IDs, complete numeric finite values, expected approximate probe-row count, and sample-level observed intensity ranges. Report all values, original titles, units, groups and evidence without replacing failures with fabricated source data.
4. **Explicitly block** `Nutrimouse → GSE22871` model projection. Nutrimouse features are liver macroarray **gene symbols** among a preselected 120-feature panel; GSE22871 measurements are Affymetrix **probe IDs** on a different scale. Without a **versioned, independently checked probe-to-gene annotation**, mapping of multiple probes to one gene, and a frozen cross-platform scaling/normalization bridge, there is **no valid direct external score**.
5. No genotype-label-driven post hoc selection, no imputation claimed as raw, no assumed same physical sample across studies, no pooled pseudo-replicates, and **no patient-level PK data or concentrations**. The output declares `BLOCKED_UNTIL_ANNOTATION_AND_SCALE_BRIDGE` and zero exported per-sample signature scores by design.

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
```

All reproducibility and publication claims remain tied to the exact successful GitHub SHA and the original GEO gzip SHA. This work updates only MultiOmics, not the source/consumer repositories.
