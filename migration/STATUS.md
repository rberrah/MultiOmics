# Migration status — MultiOmics autonomous repository

**Date:** 2026-10-16  
**Status:** scientific source copy verified; autonomous application and hosting **not yet qualified**.

## Reproducible source snapshot

- Source: `rberrah/rberrah.github.io` at commit `dbe503f98cfde460f95a60617f3d743adb6b915b` (merged PR #112).
- Target: `rberrah/MultiOmics`, branch `migration/import-validated-multiomics-20261016`.
- Copied **244 of 244** selected source files, compared by immutable Git blob SHA against the original; **zero mismatches**.
- Complete file-by-file audit: [source-snapshot-manifest.json](source-snapshot-manifest.json).
- Sources stay in their old directory layout temporarily to preserve import paths and simplify numerical comparisons.

| Source family | Imported files |
| --- | ---: |
| JS analysis and researcher ZIP modules | 25 |
| Svelte MultiOmics routes | 3 |
| Scientific Svelte components | 3 |
| R scientific engine, tests and validators | 73 |
| Static demos, templates and explanations | 29 |
| Numeric, public-study and export test scripts | 35 |
| Browser acceptance tests | 28 |
| Scientific methods documentation | 10 |
| Prior multi-omics GitHub workflow definitions | 38 |
| **Total source-matched** | **244** |

The previous 38 GitHub Actions workflows are preserved under `migration/legacy-workflows/` **not** under `.github/workflows/`. They have not been activated because many still assume the PMx parent site, its deploy paths, environment and unrelated laboratory modules.

## Required before independent release

1. **Standalone SvelteKit shell**. Create a minimal `package.json`, lockfile, application layout, theme/language store, Vite/Svelte config and root entry. Do not copy the entire PMx portal into this repository. Ensure `/multiomics/`, `/multiomics/tool/` and `/multiomics/learn/` work.
2. **Rebuild and test** from a clean checkout, no implicit imports from missing `$lib` or parent-portal files. Make the four demo buttons appear in correct order on desktop and mobile, with hydration and no missing asset requests.
3. **Independent CI**. Port only the meaningful workflows with fresh `.github/workflows/` definitions and explicit dependencies: Node/Svelte, R/Bioconductor, Python exchange checker, Playwright, code/source hashes. Run the actual public Airway/GSE52778 paired study with 17,199 features; compare all native R results, JS/R plotted 30×8 values, ZIP replay and an HTML report knitted from the downloaded Rmd.
4. **R backend selection**. Native reference R models run locally or in an explicitly authorized server; GitHub Pages hosts static HTML/JS only. Never imply a public GitHub Pages build executes R on its own. Never upload private patient matrices to an untrusted endpoint.
5. **Hosting with no paid plan**. Enable static Pages for `rberrah/MultiOmics` only after a green standalone build. Use repository Pages base path `/MultiOmics` unless deploying via a specific custom domain. Keep PMx Explain routes intact until the standalone origin passes real public-site tests.
6. **Repository publication metadata**. Audit licensing, attribution, relevant non-patient datasets, AUTHORS/CONTRIBUTING, `CITATION.cff`, package/version pins and reproducibility instructions. Do not copy PMx/TDM citation claims automatically.
7. **Scientific limits**. Calibrated inference remains **method-, design- and data-specific**. FDR inflation in some null benchmarks, MOFA2 normalization sensitivity and uncertain DIABLO benefit must remain reported. Passing code tests is not universal confirmatory certification.
8. **Final migration audit**. Compare the old vs new sites on real demo uploads, external-data boundaries, numeric outputs, figure data, data availability, novice usability and security. Only afterwards replace source routes with compatibility redirects.

## Intentional boundaries

- **No production cutover** performed.
- **No claim that the extracted source has passed a standalone clean install/build.**
- **No new external R service** enabled.
- The old site remains the sole reference public deployment until audited cutover.
- The new repo is public to facilitate sharing and the eventual manuscript, not yet labeled a reproducible release.

See `docs/multiomics/REPOSITORY_EXTRACTION_ADR.md` for the staged extraction decision.