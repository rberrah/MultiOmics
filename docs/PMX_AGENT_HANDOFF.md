# MultiOmics → PMx Explain: source ownership and URL release contract

This contract is maintained in **rberrah/MultiOmics** only.
It documents the expected integration for the *separate* PMx Explain agent;
**no PMx source, workflow or deployment is modified here**.

## Authoritative source and provenance
- MultiOmics is the authoritative source for all omics UI, analytic JS, reference R,
  method tests and beginner demos. Do not continue editing a forked copy under PMx.
- A public PMx release must use **a full 40-digit immutable, scientifically
  qualified MultiOmics commit SHA**. Never consume a moving main branch or
  unreviewed artifact in production.
- The source emits `multiomics-fixed-urls-SHA` and
  `multiomics-root-short-url-SHA` artifacts. These are build **inputs**,
  **not** autonomous GitHub Pages deployments.
- Both payloads contain a per-file SHA-256 manifest; verify all files, reject
  extras, and allow only the declared directories before copying.

## Public URL ownership
```text
https://rberrah.github.io/multiomics/             new preferred URL
https://rberrah.github.io/multiomics/tool/
https://rberrah.github.io/multiomics/learn/

https://rberrah.github.io/pharmacometrie/multiomics/        legacy, keep live
https://rberrah.github.io/pharmacometrie/multiomics/tool/
https://rberrah.github.io/pharmacometrie/multiomics/learn/
```
- Direct root build: `BASE_PATH` unset, deploy ONLY `multiomics/` and
  `_multiomics_app/` from `short-url-overlay/site/`.
- Existing nested build: `BASE_PATH=/pharmacometrie`, deploy ONLY
  `pharmacometrie/multiomics/` and
  `pharmacometrie/_multiomics_app/` from `compat-overlay/site/`.
- Do **not** copy standalone `index.html` to the portal. Do not change PMx
  `_app/`, statistics, courses, TDM, communication or unrelated web routes.
- A static redirect, iframe or renamed folder is **not equivalent** to testing
  the live app's four demo controls, links, hydration and direct deep links.
- Keep the old nested routes valid as real apps or a tested canonical redirect;
  preserve query strings and anchors where redirecting.

## Release gates
1. Full standalone SvelteKit + JavaScript/R contracts + desktop/mobile Playwright.
2. Airway raw-count DESeq2 numerical golden replay, R Markdown and figure table checks.
3. Real Nutrimouse DIABLO and MOFA2 with negative/mispaired controls.
4. Pinned R/Bioconductor native global-null and mixed-truth stress reports.
5. Both SHA-256 overlay validation checks.
6. **On the final PMx Pages artifact**: verify both URL families, 404 free
   assets, MS-AUC CSV parsing and research ZIP integrity in a real browser.
7. Merge PMx integration **only after** all upstream checks finish green and
   the deploying site CI passes. A successful source build does not guarantee
   a successful deployment or valid biological inference.

## Operational and methodological boundary
- GitHub Pages runs **browser-only exploratory JS**; it does not execute
  DESeq2, MOFA2, mixOmics/DIABLO or other R/Bioconductor methods natively.
- Explicit, authorized local or institutional R backend is needed for R methods.
  Do not silently fall back from reference R to JS in confirmatory mode.
- Archival research artifacts must record study inputs, group contrast,
  package versions, seeds, filters, missingness decisions, SHA-256 provenance.
- Do not erase old PMx copies or tests before dual-route smoke passes and a
  recovery mechanism exists. Coordinate deletion separately with PMx agents.

Version and ownership of this document are **independent of PMx Explain**.
