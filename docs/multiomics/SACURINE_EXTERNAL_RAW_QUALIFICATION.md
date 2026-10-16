# External LC-MS validation gate: Sacurine / MTBLS404

## What is independently supported

The [official Workflow4Metabolomics Sacurine public source](https://workflow4metabolomics.org/lcms_dataset-sacurine)
describes negative-ion LTQ-Orbitrap LC-HRMS, with a 36-file,
approximately 188-MB mzXML subset: **20 biological samples, 10 pooled
quality-control injections, 6 blank injections**. The full two-batch
dataset is substantially larger. Reference: Thévenot et al.,
*J Proteome Res* (2015), DOI 10.1021/acs.jproteome.5b00354.

Sacurine is independent of the earlier native faahKO CDF fixture and
the measured Drosophila peak-area validation. Unlike those experiments,
the selected archive includes actual blank and pooled QC acquisitions.

## Phase A: real-source qualification (this PR)

- Read official files from the actual research provider:
  `https://download.workflow4metabolomics.org/datasets/sacurine-neg-subset_mzXML.zip`
  and `sacurine-neg_sampleMetadata.tsv`.
- Stream archive to the CI runner with a 300 MB size cap; fail on a
  missing/empty source, unsafe ZIP traversal, symlinks, duplicate
  chromatograms, unsupported file names or incorrect 20:10:6 roles.
- Record observed archive SHA-256, raw metadata SHA-256, each ZIP member
  size/CRC, and which original metadata columns match the archive file
  names. Do not infer a full injection sequence, technical batch or
  clinical subgroup from alphabetical file ordering.
- The source license, attribution and original paper remain with the
  original provider. No original raw files or individual-level clinical
  annotations are committed to this repository or uploaded as CI
  artifacts.
- ZIP integrity and filename class are **preliminary checks**, not
  successful xcms peak extraction, not measured blank contamination,
  not proof of instrument-batch balance, and not a statistical finding.
- The actual SHA-256 **must be captured from the successful first run
  and subsequently pinned fail-closed** before the raw scientific
  processing workflow can be called reproducible.

## Observed public-source availability (16 October 2026)

This is a recorded infrastructure blocker, **not** a failure of `xcms`:

- Workflow4Metabolomics download server:
  [real run 38924316667](https://github.com/rberrah/MultiOmics/actions/runs/38924316667)
  could not connect to its official `sacurine-neg_sampleMetadata.tsv`
  after three attempts (network timeouts); it never fetched the 188-MB
  archive or processed chromatograms.
- EMBL-EBI primary metadata check:
  [run 38925713284](https://github.com/rberrah/MultiOmics/actions/runs/38925713284)
  returned HTTP **404** for the two documented `s_sacurine.txt`
  paths tried; these endpoints cannot be treated as authenticated
  scientific input without new official information.
- Offline archive role/path attacks and the **local archive preflight**
  are tested separately and may pass; those results do not imply a
  successful download or external laboratory validation.
- Network jobs now require explicit `workflow_dispatch`, so temporary
  external-repository outages cannot turn routine standalone software
  CI red or encourage silently substituting incompatible data.

### Running preflight with a genuinely acquired local research archive

If the primary research provider becomes accessible through a browser,
a laboratory workstation, or an institutional data transfer route, use
the original source files **locally** rather than committing them:

```bash
python3 scripts/sacurine_source_preflight.py \
  --local-archive /path/to/sacurine-neg-subset_mzXML.zip \
  --local-metadata /path/to/sacurine-neg_sampleMetadata.tsv \
  --output /path/to/empty/sacurine-preflight-evidence
```

This path verifies complete 36-file ZIP structure, expected provisionally
classified 20/10/6 roles, archive SHA-256, source metadata checksum,
original sample identifier candidates, and ZIP traversal safety.
The command **does not delete the two source inputs**. Its compact
reports are not a replacement for checking source-provided sample
classes and injection-order fields before the full R native run.

## Phase B: scientific raw processing (next gated PR)

After Phase A passes and source roles are reconciled to original metadata:

1. Extract the exact 36 listed mzXML files into a temporary CI
   working directory with strict size and path checks.
2. Run the **already qualified R 4.5.2/Bioconductor 3.22 xcms stack**
   (xcms 4.8.0 and exact 206-package closure), with prespecified
   centWave/retention-alignment/peak-grouping parameters.
3. Compare native RDS objects, feature definitions, CSV matrices and
   sample identities; preserve version inventories and source checksums.
4. Apply the existing descriptive MS quality calculations:
   pooled-QC CV, per-feature blank/bio median ratio and explicit
   missingness/observability. Drift and batch effects only if
   verified injection-order and batch metadata permit them.
5. Archive compact quantitative reports; never upload the original
   188-MB raw files to GitHub or treat nominal discovery thresholds as
   automatically calibrated.
6. Quantitatively compare against the W4M Sacurine reference
   `W4M00002` with harmonized feature definitions and parameter
   sensitivity. Absolute feature counts alone are not equivalence.

## Scientific maturity

Preflight+authentic 36-file raw processing+independent comparator
will support claims of portable reproducible workflow and instrumental
QC detection for **one** additional external laboratory. Generalization
requires more independent instruments, cohorts and public sources.

See `docs/multiomics/EXTERNAL_PUBLIC_VALIDATION_PROTOCOL.md` for
the independent Drosophila experiments, observed measures and the
publication evidence matrix.
