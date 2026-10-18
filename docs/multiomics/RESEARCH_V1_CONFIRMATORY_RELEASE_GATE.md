# MultiOmics — gate for an interpretable research v1 release

**Scope:** independent `rberrah/MultiOmics` source repository. This document
states **software acceptance and deployment handoff**, not universal
statistical or biological validation.

## Qualified starting point

At 18 October 2026, `main` commit
`c4c1c51346d4b05acd70754e7680b1f6540fdae8`
passed all eight post-merge independent GitHub Actions qualifications
(standalone UI, fixed-URL payload, R reference calibration, public Airway native R,
real Nutrimouse MOFA2/DIABLO, GSE22871 and Omics2PK exchange).

The public PMx Pages portal is a **separately pinned application**:
`rberrah/rberrah.github.io/integration/multiomics-release.json` still pointed
to source commit `b02cc48e006abf5fd0200e9e2cf376af4cfae526`
when inspected. A successful MultiOmics CI is **not evidence** the new source
is deployed at `https://rberrah.github.io/multiomics/`. Changing that pin,
running assembled Pages tests, merging the separate PMx PR and verifying the
public site must be performed by the PMx integration owner.

## Hard acceptance rules for the browser's confirmatory mode

- The user deliberately selects confirmatory rather than exploration.
- A local/authorized R backend must satisfy the exact expected health
  contract, and have a potentially suitable R package family installed.
  HTTP 200 alone or a simulated `/health` response is insufficient.
- **Groups:** an actual native R model must complete and return a *complete
  feature-wise table for each loaded omics layer*. The displayed results must
  be projected from each layer's exact returned R CSV, with matching counts
  and method key. An unrelated R model (fgsea, PCA, MOFA2 or another layer)
  cannot substitute for a missing differential result.
- **Longitudinal group-by-time:** the browser may enter confirmatory display
  only when every loaded layer has a complete native `lmerTest` table from
  the strict repeated-measures contract. The displayed effect is the explicit
  condition-by-time slope coefficient. Missing visits, singular fits,
  random-slope fallback and the linear-time assumption remain visible and
  require scientific review.
- **Survival outcome:** the browser may enter confirmatory display only when
  every loaded layer has a complete native `survival::coxph` table from the
  strict one-observation-per-subject contract. Event counts, Efron tie policy
  and `cox.zph` diagnostics are retained. A non-significant PH diagnostic is
  not proof of proportional hazards.
- **Other outcome types:** continuous, binary, count and multiclass
  association tables remain browser-exploratory. DIABLO, an HTTP 200 response,
  or the presence of an R package cannot substitute for a dedicated native
  feature-wise inference contract.
- A valid R result is **still not a certified statistical analysis**:
  design rank, biological units, data provenance, appropriate contrasts,
  missingness, batch effects, model assumptions and multiplicity remain
  research decisions. Independent review and validation of conclusions remain
  necessary.
- No fake `training_frozen_test_projection` scores may be generated for
  Omics2PK. Descriptive GSE22871 native probe profiles do not constitute
  externally validated predictors or measured concentrations.

## Regression evidence added in this release candidate

- `scripts/test_multiomics_reference_backend_preflight.mjs` checks backend
  identity and relevant R dependencies, including mixed raw-count RNA-seq +
  processed omics when limma is missing.
- `scripts/test_multiomics_confirmatory_result_gate.mjs` checks exact
  per-layer R projection and blocks partial, forged or unrelated responses
  across group, longitudinal and survival contracts.
- `scripts/test_multiomics_native_projection_designs.mjs` independently
  checks that native lmerTest/Cox CSV values replace browser numbers only for
  the supported objective and that incomplete tables fail closed.
- `tests/e2e/multiomics-confirmatory-native-gate.spec.js` checks that an
  unsupported longitudinal confirmatory workflow is blocked **before**
  transferring matrices to the reference backend.
- Existing `tests/e2e/multiomics-r-health-preflight.spec.js` verifies
  that arbitrary services returning 200 OK cannot be mistaken for the
  MultiOmics R backend.
- Mandatory `standalone-ci.yml` runs the negative controls with no
  retry-based success policy. These tests do not replace independent
  Bioconductor method qualification.

## Handoff to PMx deployment owner

1. Wait for the exact **squash-merged** source SHA to pass all required
   independent workflows; do not pin a PR branch or a stale passing SHA.
2. Propose a minimal update in the PMx-owned
   `integration/multiomics-release.json` to that immutable source SHA.
3. PMx verifies source workflows and builds the assembled Pages artefact,
   preserving `/multiomics/`, `/multiomics/tool/`, FR/EN, theme, mobile nav,
   and the old `/pharmacometrie/multiomics/` aliases.
4. Check the actual served Pages URL (hydration, raw MS CSV import,
   user-generated ZIP, local R guide, no falsely confirmed output).
5. Only then state that the latest MultiOmics is publicly deployed.

**Outstanding v1 release decisions:** make a permanent versioned research
release and frozen environment reproducible outside CI; expand R reference
qualification to new designs, MS batches with jointly observed blanks/QCs;
perform human usability testing; investigate null-simulation false calls with
sufficient Monte Carlo precision. See the statistical qualification
documentation rather than describing current results as universally
confirmatory.
