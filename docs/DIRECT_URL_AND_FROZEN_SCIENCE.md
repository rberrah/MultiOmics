# Direct /multiomics/ URL and reproducible reference environment

## 1. Direct public URL without breaking former links

**Target canonical path:** `https://rberrah.github.io/multiomics/`
with `/multiomics/tool/` and `/multiomics/learn/`.

**Current production:** the `/multiomics/` address is an existing portal
alias pointing to `/pharmacometrie/multiomics/`. GitHub Pages cannot
transparently reverse-proxy this without a new site artifact. Simply
changing the MultiOmics navigation would break deep links.

This source repository prepares **two separate, SHA-256 verified** builds:

- **Legacy PMx root:** `BASE_PATH=/pharmacometrie`, exporting only
  `pharmacometrie/multiomics/` and `pharmacometrie/_multiomics_app/`.
- **Direct domain root:** `BASE_PATH` unset, exporting only
  `multiomics/` and `_multiomics_app/`. This is a real hydrated static
  application, not a meta refresh or an iframe.

Run `.github/workflows/pmx-fixed-url-artifact.yml` to obtain both
**non-deploying** artifacts and manifests, or locally:

```bash
npm ci
BASE_PATH=/pharmacometrie npm run build
node scripts/build_pmx_compat_overlay.mjs
node scripts/test_pmx_compat_overlay.mjs
env -u BASE_PATH npm run build
node scripts/build_short_url_overlay.mjs
node scripts/test_short_url_overlay.mjs
```

For the production switch, PMx Explain's deploying repository must, in an
isolated integration PR, overlay the root `multiomics/` and root
`_multiomics_app/` directories **after** assembling the existing portal.
Do not copy the root `index.html` or unrelated assets. Retain a legacy
`/pharmacometrie/multiomics/` functional application or safe alias to the
new canonical route; preserve deep links, query strings and fragments.
The old URL should not disappear, because it may be referenced by
scientific communications and saved research projects.

**This repository does not change PMx Explain, its concurrent agents, or
its live deployment.** The root URL is only ready for release, not deployed,
until that reviewed PMx PR passes the final assembled-site browser
acceptance and public smoke.

## 2. Fixed reference runtime (independent of browser UI)

Canonical R validation baseline:

- R **4.5.2**
- Bioconductor **3.22** (official pairing for R 4.5)
- Immutable OCI image by full digest:
  `bioconductor/bioconductor_docker:RELEASE_3_22-r-4.5.2@sha256:e34177896c1579f66d48b756e4229b58237499183a089da084af72d7aba1d40a`

The machine-readable contract lives under
`multiomics-engine/environment/reference-runtime.json`. The verifier
checks R/Bioconductor versions, required packages and, once populated
from actual CI package installation, **exact installed package
versions**. The entire installed package inventory is saved with each
scientific job to detect otherwise invisible dependency drift. It is
not appropriate to invent package versions: the contract is populated
from an actual reproducible CI environment and then enforced fail-closed.

An immutable image locks its original OS/R contents. **Installing
previously missing R packages from moving CRAN/Bioconductor indexes is
not reproducible unless those packages are independently pinned as well.**
Do not describe this as a complete `renv.lock` or an offline
Bioconductor distribution. A full transitive lock / binary mirror can
be added independently if exact offline reinstallation is required.

## 3. Extended methods qualification — interface unchanged

Separate from exploratory browser tests:

- Genuine **DESeq2**, **edgeR quasi-likelihood** and **limma-voom** under
  prespecified all-null negative-binomial synthetic studies. Report the
  number of false discoveries and exact Monte Carlo uncertainty. A
  regression threshold is not a proof of 5% FDR.
- The same three native R models under planted nonzero effects (1%, 5%,
  20% of features), with FDR, power, effect direction and CI coverage.
  No method is selected post hoc based on apparent benchmark performance.
- Real **Airway** native R vs full tables and R Markdown replay.
- Real **Nutrimouse** MOFA2 and DIABLO with matched animals, label
  permutations and single-omic ablation tests, after pinning the
  compatible R/Bioconductor release.
- Fast package-free numerical/contract cases reject false confidence
  from confounding, incomplete RNA-seq sample manifests, heldout
  leakage, unsupported FDR interpretations, low-quality MS and
  potential missing-not-at-random proteomics data.

Source scientific results, including unflattering null calibration or
negative MOFA/DIABLO signals, must remain published as evidence. This
does **not** establish global Type-I/FDR calibration or biological/clinical
validity for arbitrary input datasets. The beginner interface and its
single-omic/three-omic demonstrations are not modified.
