# MultiOmics

An open-source, browser-first environment for reproducible **single- and multi-omics research**, with an optional local R reference engine.

**Release status:** autonomous migration under review — this is **not** yet a validated standalone release or a universal confirmatory analysis platform. The upstream public tool remains available at [rberrah.github.io/pharmacometrie/multiomics/tool](https://rberrah.github.io/pharmacometrie/multiomics/tool/).

## Analyse in a browser

MultiOmics supports transcriptomics, proteomics and metabolomics (one, two or three layers). Import sample metadata and measurement matrices, or try the four built-in synthetic demos. The three individual-omics buttons sit above a wider multi-omics button.

JavaScript models are **deterministic exploratory calculations** and may differ materially from suitable reference methods. For supported scientific inference, select a valid study design and connect the real R backend. This is not a claim that a low adjusted p-value proves causality, a calibrated 5% FDR under all designs or clinical usefulness.

## Run the autonomous website locally

Requires Node.js 22 and npm.

```bash
npm ci
npm run check
npm run dev
```

Open the local address shown by Vite and choose `/multiomics/`, `/multiomics/tool/` or `/multiomics/learn/`. The site requires **no PMx Explain checkout**.

```bash
npm run test:multiomics
npm run build
npx playwright install chromium
npx playwright test --config playwright.multiomics.config.js tests/e2e/multiomics-methods-learning.spec.js
```

GitHub Pages is static: `BASE_PATH=/MultiOmics npm run build` is the deployment target when the public release has been independently qualified. A static Pages site **does not execute R on the server**.

## Optional native R reference engine

Install R and the required packages on your own machine:

```bash
Rscript multiomics-engine/install_backend_dependencies.R
Rscript multiomics-engine/run_backend.R
```

Use the **local** engine (default `http://127.0.0.1:8787`) only after checking its health in the UI. In formal/reference mode, a missing or incompatible backend must **block** native inference; it must not be silently replaced with exploratory JavaScript. Do not enter a remote R URL for confidential datasets unless institutional authorization and a privacy review exist. Consult [engine instructions](multiomics-engine/README.md) and [scientific limitations](static/multiomics/scientific-validation.md).

## Export as a reproducible research project

The research ZIP contains original inputs, complete native R result tables **when actually available**, source/version provenance and checksums, and a commented R Markdown report.

After extracting a ZIP into a local directory:

```bash
Rscript run_all.R
Rscript render_notebook.R
```

The second command installs missing **reporting** dependencies where supported, requires Pandoc, reruns supported reference R calculations and creates `analyse_reproductible.html`. Missing/mismatched scientific methods are **not** replaced with whatever package version is latest. Browser-only results are explicitly **not equivalent** to native R reference values.

## Scientific evidence and limits

- **Public Airway / GSE52778:** four biological cell lines measured twice; independent DESeq2 numerical comparison of **17,199 fitted genes**, native R vs website-plotted data (30 features × 8 samples) and archived research ZIP/Rmd replay in the original integrated CI. The autonomous run must independently pass before release.
- **MOFA2:** exploratory latent factors; batch adjustment, normalization and external replication matter.
- **DIABLO:** supervised learning needs patient-level splits, leakage controls and single-omic ablations. Its superiority is not universal.
- **False discovery rate:** stress simulations revealed possible inflation in some designs; an algorithm running successfully is not proof of 5% frequentist control. Document assumptions and uncertainty.

[Scientific release gates](docs/multiomics/RESEARCH_RELEASE_GATES.md) · [Public Airway evidence](docs/multiomics/PUBLIC_AIRWAY_SITE_R_REPLAY.md) · [Migration status](migration/STATUS.md).

## Provenance

The initial scientific snapshot contains **244 source-matched files** from `rberrah/rberrah.github.io` commit `dbe503f98cfde460f95a60617f3d743adb6b915b`, verified against immutable Git blob hashes. See `migration/source-snapshot-manifest.json`. Standalone adaptations are intentionally separate commits.

The original portal's MultiOmics links remain operational while the autonomous version is tested. No deletion, redirect cutover, scientific publication claim or clinical deployment is implied by opening this repository.
