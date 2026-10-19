# Sacurine v3 — independent xcms vs OpenMS preprocessing comparator

**Protocol frozen on 19 October 2026 before the OpenMS result is inspected.**

This experiment follows issue #47. The authentic 12-file Sacurine v3 xcms execution is already qualified. This protocol asks a different question:

> Do independently implemented LC-MS feature-finding pipelines recover overlapping m/z–RT feature coordinates on the exact same authenticated raw injections?

It is **not** a metabolite-identification benchmark, a clinical validation, or a claim that one pipeline is superior.

## Source

Exactly the same immutable Zenodo v3 source already used by the xcms qualification:

- 6 biological urine injections;
- 3 pooled-QC injections;
- 3 analytical blanks;
- original source metadata and injection positions;
- all 12 mzML publisher MD5 values verified before either pipeline runs.

No raw mzML is committed or uploaded as an Actions artefact.

## Pipeline A — frozen MultiOmics xcms baseline

The current qualified xcms path is replayed on the 12 files with its previously frozen parameters. The comparison exports only:

- xcms feature ID;
- native median m/z `mzmed`;
- native median retention time `rtmed`.

The existing negative QC findings remain unchanged and are not filtered away.

## Pipeline B — independent pinned OpenMS build

The independent implementation uses the immutable OpenMS container digest already frozen in the protocol. Runtime inspection of the first successful source-processing artifact showed that this image identifies itself as **OpenMS 3.5.0-pre-HEAD-2025-12-12, Revision 92e4db6**, not a pure 3.5.0 release. This provenance correction was made after reading the runtime identity and does **not** change feature-finding, alignment, linking, matching tolerances or any observed result. The workflow now verifies that exact build string in addition to the pinned container digest.

The pipeline is deliberately generic and non-optimized:

1. `FeatureFinderMetabo` independently on every original mzML;
2. `MapAlignerPoseClustering` across the **9 non-blank feature maps** (6 biological + 3 pooled-QC); all 12 acquisitions are still authenticated and feature-found;
3. `FeatureLinkerUnlabeledKD` across those same 9 non-blank aligned maps to obtain the comparison consensus. Analytical blanks are excluded from consensus definition so blank-only signals cannot enlarge the biological comparison universe; the blank acquisitions remain preserved in the authenticated source and existing xcms QC evidence.

For reproducibility, each OpenMS tool writes its own complete default INI under this exact version and every subsequent command consumes that frozen INI. Parameters are therefore version-frozen defaults, not tuned after looking at Sacurine overlap.

**No silent profile→centroid conversion is allowed.** If FeatureFinderMetabo refuses the source because the mzML is not suitable centroided input, the workflow fails and that negative result is retained. A later centroiding protocol would need a separate explicit amendment.

## Predeclared coordinate comparison

The consensus OpenMS centroid coordinates are parsed from the native `consensusXML`; xcms feature IDs are never equated to OpenMS IDs.

Matching rule:

- m/z tolerance: **20 ppm**;
- RT tolerance: **30 seconds**;
- a pair is retained only if it is the **unique candidate in both directions**;
- zero candidates, multiple candidates and non-reciprocal candidates remain separate failure categories;
- no nearest-neighbour rescue;
- both denominators are reported: fraction of xcms features and fraction of OpenMS consensus features in reciprocal unique pairs.

The complete xcms-side match table is archived, together with unmatched/ambiguous counts and both pipeline version/parameter fingerprints.

There is **no predeclared minimum overlap required for PASS**. A green workflow proves that the two independent pipelines executed on the same authenticated source and that the comparison was reproduced. A low overlap is a valid scientific result.

## Interpretation limits

Even high coordinate concordance does not establish:

- chemical identity;
- isotope/adduct equivalence;
- quantitative concentration accuracy;
- validated blank subtraction, drift correction or carryover;
- interlaboratory reproducibility;
- a universally valid downstream feature-exclusion threshold.

These claims continue to require standards, annotated compounds and/or additional laboratories/instruments.
