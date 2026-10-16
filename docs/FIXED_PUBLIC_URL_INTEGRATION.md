**Current deployment (16 October 2026):** This integration has been implemented by `rberrah/rberrah.github.io` PR #116 (commit `f81e2748`) and monitored by PR #117 (commit `1320880a`). The public deployment of PR #116 was successful. PMx Explain owns Pages, pins `MultiOmics@9041d38c` and builds only the allowlisted `/pharmacometrie/multiomics/` and `/_multiomics_app/` overlay. A new MultiOmics/main commit does NOT automatically change the public app: coordinate an immutable-pin update through PMx CI. The below cutover section is the **historical design and rollback recipe**, not an outstanding first-deployment instruction.

---

# Same public URLs, independent MultiOmics repository

Decision (2026-10-16): the repository rberrah/MultiOmics is the future
authoritative application source, but **the public origin and paths do not
move**:

- https://rberrah.github.io/pharmacometrie/multiomics/
- https://rberrah.github.io/pharmacometrie/multiomics/tool/
- https://rberrah.github.io/pharmacometrie/multiomics/learn/

Other MultiOmics assets and legacy routes must continue to work.
No redirect to a separate GitHub Pages project URL is part of this migration.

## Isolation from the PMx Explain team

This integration is implemented exclusively in the MultiOmics repository.
It does NOT modify rberrah/rberrah.github.io/main, PMx CI, concurrent
PMx Explain PRs, or GitHub Pages settings. A later, small and explicitly
reviewed PMx integration PR must be based on the then-current PMx head.

The standalone SvelteKit app uses _multiomics_app instead of _app,
so its compiled asset bundle cannot overwrite the PMx Explain asset bundle
even when both are published on the same rberrah.github.io origin.

## Produce a non-deploying integration payload

    npm ci
    BASE_PATH=/pharmacometrie npm run build
    node scripts/build_pmx_compat_overlay.mjs
    node scripts/test_pmx_compat_overlay.mjs

Output:

    compat-overlay/
      manifest.json
      site/
        pharmacometrie/
          multiomics/         # prerendered pages, examples, template files
          _multiomics_app/    # isolated hashed SvelteKit assets

The manifest gives the relative path, bytes and SHA-256 of each file.
The packaging script rejects files outside the two approved directories
and checks that HTML uses /pharmacometrie/_multiomics_app/ instead of
the shared PMx /pharmacometrie/_app/ namespace.

The independent GitHub Actions workflow pmx-fixed-url-artifact.yml builds
and tests this payload without uploading it to GitHub Pages. Only its CI
artifact is available for inspection.

## Controlled future cutover (NOT performed here)

After scientific, privacy, desktop and mobile browser qualification, create
a **separate PMx deployment PR**. Pin an immutable MultiOmics release SHA
and build it as above, or verify an equivalent authentic artifact. Overlay
only dist/pharmacometrie/multiomics/ and
dist/pharmacometrie/_multiomics_app/ onto the existing PMx deployment.

The current PMx deploy-pages.yml copies portal/ to dist/, then its full
build/ to dist/pharmacometrie/. The two namespaced directories belong
AFTER these existing copy steps and BEFORE the public-artifact audit and
Pages upload. Never copy the MultiOmics root index.html or its full build/
over the PMx site. Keep the legacy portal aliases and rollback procedure.

Public acceptance: the existing three small mono-omics demonstration
buttons and one large multi-omics button, MS CSV import stages,
FR/EN desktop/mobile, keyboard use, valid navigation and examples, full
ZIP/Rmd replay and data tables, privacy/no surprise network transfer,
and the real deployed URLs after GitHub Pages success.

## Scientific and rollout constraints

GitHub Pages serves static JS/HTML only. Reference R analyses need local or
explicitly approved compute. CI success does not certify universal FDR,
biological validity, clinical use or a hosted R backend.

Do not delete the original PMx route source in the cutover PR. Source removal
is separate and gated by stable public URLs, live smoke tests and independent
MultiOmics releases. Until then both source trees may coexist temporarily.
