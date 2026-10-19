# MultiOmics — gate for an interpretable research v1 release

**Scope:** independent `rberrah/MultiOmics` source repository. This document
states **software acceptance and deployment handoff**, not universal
statistical or biological validation.

## Research v1 freeze candidate

Research v1 is accepted only from the final merged source revision whose
release gate, standalone qualification and required public-data evidence are
green. The immutable source identity is the final `main` SHA; historical
passing SHAs remain evidence but are not substituted for the release revision.

The frozen scientific scope is deliberately narrower than “universally
validated omics analysis”. Native-R confirmatory **preparation** is supported
only for the explicit group, repeated group×time and independent-survival
contracts below. Public real-data qualification includes Airway, Nutrimouse,
the authentic 12-file Sacurine raw-MS study, NKI/DMFS survival and GSE191328
longitudinal execution. Known limitations #44, #47 and #20 remain visible and
are part of the release contract rather than being hidden or closed by the
software freeze.

The public PMx Pages portal remains a **separately pinned application**.
A successful MultiOmics Research v1 freeze does not prove that the same SHA is
already deployed at `https://rberrah.github.io/multiomics/`. The PMx
integration owner must update its immutable source pin, run assembled-site
tests and verify the served public URL separately.

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
- `tests/e2e/multiomics-confirmatory-native-gate.spec.js` checks that a
  supported longitudinal confirmatory workflow is exposed only through the
  required native-R path and is no longer mislabeled as unsupported.
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

## Research v1 release decision

The software may be frozen as **Research v1 / 1.0.0** once the exact final
source revision passes the release CI and its required public-data evidence.
The freeze does not close the following scientific work:

- #44 remains a known limitation on universal DESeq2/FDR claims;
- #47 remains stronger analytical/chemical LC-MS validation beyond the already
  qualified authentic raw-data execution;
- #20 remains biological transferability of individual scores toward Omics2PK;
- independent human usability and study-specific statistical/biological review
  remain external tasks.

Those items block broader scientific claims, not the bounded Research v1
software/provenance baseline. A later method, dependency or scope change must
receive a new version and requalification rather than silently modifying v1.
