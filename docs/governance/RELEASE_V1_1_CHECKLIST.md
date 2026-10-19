# MultiOmics 1.1.0 — final release checklist

This checklist begins only after PR #71 and the dedicated 1.1 release PR are green on the exact final source combination.

## Repository freeze

- [ ] Merge PR #71 only after the synchronized xcms/OpenMS comparator, standalone qualification and Omics2PK consumer contract are green.
- [ ] Merge the 1.1 release PR only after all release gates are green.
- [ ] Record the final merged `main` SHA.
- [ ] Confirm `release/research-v1.0.0` still resolves to `54dc8f0f8284cf970ed1a5be0533a298ae27639d`.
- [ ] Confirm `package.json`, `package-lock.json` and `CITATION.cff` all report `1.1.0`.
- [ ] Do not add new methods, datasets or scientific claims to 1.1.0 after this point.

## GitHub administration

The connected automation used to prepare this release cannot write repository rulesets or create GitHub Releases/tags. These two actions require repository administration in GitHub.

- [ ] Import `docs/governance/github-main-ruleset-v1.1.0.json` as an active ruleset for the default branch.
- [ ] Verify pull requests are required, force-push/deletion are blocked, and the listed checks are recognized by GitHub.
- [ ] Create immutable tag `v1.1.0` at the final merged 1.1.0 SHA.
- [ ] Create a GitHub Release from `v1.1.0`, using `docs/multiomics/RELEASE_V1_1_0_20261019.md` as the scientific release boundary.
- [ ] Do not move or recreate `v1.1.0` after publication.

## Public-site handoff

- [ ] In `rberrah/rberrah.github.io`, replace the old MultiOmics pin `b02cc48e006abf5fd0200e9e2cf376af4cfae526` with the exact final 1.1.0 SHA.
- [ ] Keep the public route `/multiomics/` and legacy aliases.
- [ ] Require the upstream qualification checks defined by the public-site integration manifest.
- [ ] Re-run mount/build/browser E2E.
- [ ] Verify the visible results page states the numerical source/evidence tier and **structural exchange verified · not PK-ready** where appropriate.
- [ ] Verify the research ZIP still passes the independent Omics2PK consumer contract.
- [ ] Close #72 only after the pinned public deployment is verified.

## Scientific limitations retained after freeze

Keep these issues open unless new evidence genuinely resolves them:

- **#44:** DESeq2 null-calibration limitation.
- **#47:** analytical LC-MS validation beyond independent preprocessing execution.
- **#20:** external biological/predictive validation of individual omics signatures for Omics2PK.

They are research limitations/perspectives, not reasons to resume general MultiOmics feature expansion.
