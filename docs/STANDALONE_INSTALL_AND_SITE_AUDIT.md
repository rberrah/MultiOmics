# Independent install and PMx website compatibility audit

Status: qualification work on the independent MultiOmics repository.
This document distinguishes **verified baseline CI** from **new checks requiring PR CI**.

## 1. Source release: downloadable independently

- Public source: https://github.com/rberrah/MultiOmics
- Source ZIP: https://github.com/rberrah/MultiOmics/archive/refs/heads/main.zip
- Git: `git clone https://github.com/rberrah/MultiOmics.git`
- Development: Node.js 22+, npm; run `npm ci`, `npm run check`,
  `npm run test:multiomics`, `npm run dev`, or `npm run build`.
- No private file, parent PMx Explain checkout or cloud R server is required to
  run the browser edition. All sample CSV files and UI scripts are tracked.
- CI: independent SvelteKit build, full JS checks, real Playwright browser
  tests and package-free R contract tests already passed on
  `MultiOmics/main@9041d38c`.

## 2. Portable browser ZIP: runnable without npm

Build from a clean checkout with `npm ci && npm run build:portable`.
Python 3 is required **only for building the archive**. The generated output is:

    release/MultiOmics-portable-web.zip
    release/MultiOmics-portable-web.zip.sha256

Unzip and open a shell in the included `MultiOmics-portable-web` folder:

    node serve.mjs

Then browse http://127.0.0.1:4173/multiomics/tool/.

A fixed SHA-256 release manifest covers every included file, including the
web app, server and license files. The ZIP archive hash is published separately
for download verification. A portable Playwright suite opens the **download
runtime itself**, checks all three analysis/education entry points, verifies
Svelte hydration and computes real deterministic demonstrations, including
mobile widths. The server is local and read-only by default; it refuses remote
binds unless explicitly configured.

The standalone package is a statically built JavaScript browser edition.
It is **not** an R/Bioconductor executable and cannot run native MOFA2,
DIABLO, DESeq2 or other reference R packages by itself. Those analyses require
the complete source checkout, a local R installation and its dependencies.

## 3. Real R backend: separate installation and privacy boundary

From a full source checkout, after installing R/Bioconductor:

    Rscript multiomics-engine/install_backend_dependencies.R
    Rscript multiomics-engine/run_backend.R

The default HTTP endpoint is http://127.0.0.1:8787. Client-side analyses do
not gain native R capabilities merely by visiting GitHub Pages. The R engine
uses a browser-origin allowlist and refuses non-loopback binding without a
special explicit setting. Researchers should not expose confidential
patient-level data to unreviewed external endpoints.

The install helper resolves currently available R package versions; **it is
not a frozen environment lockfile**. For an exact manuscript reproduction,
record R versions, Bioconductor release, package versions, study matrices,
contrast, design decisions, seeds, and export checksums. Installation success
does not guarantee every R reference package fits every study design.

## 4. PMx Explain: original public URL contract

The PMx site builds from the independent MultiOmics SHA, presently pinned to
`9041d38c` via `rberrah.github.io/integration/multiomics-release.json`.
The source repo's `BASE_PATH=/pharmacometrie` build is transformed into the
manifest-verified overlay:

    dist/pharmacometrie/multiomics/
    dist/pharmacometrie/_multiomics_app/

The original public paths remain:

- https://rberrah.github.io/pharmacometrie/multiomics/
- https://rberrah.github.io/pharmacometrie/multiomics/tool/
- https://rberrah.github.io/pharmacometrie/multiomics/learn/

PMx Explain PRs #116 and #117 were merged: the first deploys only these
namespaced paths after independent scientific and PMx CI; the second verifies
live desktop/mobile use, original aliases, MS-AUC imports, source ZIP
integrity and privacy after deployment. A change in MultiOmics/main is **not
a silent auto-deploy** to PMx: a reviewed pin bump and PMx qualification are
necessary, preventing concurrent PMx agents' work from being overwritten.

The autonomous source repository runs a separate fixed-public-URL artifact
workflow on each update to its main branch. That workflow **never deploys**
to Pages itself.

## 5. Remaining scientific and product limits

- The deterministic JavaScript engine is exploratory; do not claim native
  reference package equivalence or universal type-I/FDR control.
- Airway / GSE52778 tests independently replay 17,199 genes on pinned
  reference CI, but this is validation of one design, not all possible
  experimental structures or patient datasets.
- MOFA2 Nutrimouse factor and DIABLO reference validations use authentic
  R libraries. They are not evidence of clinically validated biomarkers.
- The browser option to contact Reactome or another external service is
  opt-in; avoid moving researcher data unless independently authorized.
- The content/code licenses and contributions remain subject to their
  actual tracked license files, third-party terms, and a separate legal
  review of the licensing proposal in PR #4.

## 6. Release acceptance and rollback

1. Checkout the independent repository with **no PMx Explain source**.
2. Verify `npm ci`, `npm run check`, JS tests, Svelte build and all browser
   demonstrations; download and extract the **new portable ZIP**.
3. Verify the SHA-256 manifest, archive hash, URL routing, 404/read-only
   behavior, desktop/mobile hydration and demo results.
4. Run the fixed-URL overlay workflow, verify exact namespaced paths and
   SHA-256 manifest; do not upload the standalone root `index.html` to PMx.
5. Only then update the immutable source pin in a **separate PMx PR**,
   preserving all concurrently edited PMx features and test workflows.
6. Rollback is a reversion to the last proven source SHA, without changing the
   public origin or URL paths.
