// @ts-nocheck
import { test, expect } from '@playwright/test';

test('Scientific summary displays imported diagnostics literally, never as executable HTML', async ({ page }) => {
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 30000 });
  await page.getByTestId('multiomics-result-technical-detail').locator('summary').first().click();
  const panel = page.getByTestId('multiomics-scientific-summary');
  await expect(panel).toBeVisible({ timeout: 30000 });

  await page.evaluate(() => {
    window.__PMX_SCIENCE_INJECTED__ = 0;
    const marker = '<img src=x onerror="window.__PMX_SCIENCE_INJECTED__=1">';
    const diagnostics = {
      status: 'review_required',
      warnings: ['Un en-tête importé contient ' + marker],
      notes: ['Résultats à lire littéralement : < 5 observations'],
      blockers: [],
      methodEligibility: {
        mofa2: { status: 'not_recommended', sharedSubjects: 4, reasons: ['Missing ID ' + marker] },
        diablo: { status: 'not_applicable', smallestClass: 1, reasons: [] }
      }
    };
    window.__PMX_MULTIOMICS_DIAGNOSTICS__ = diagnostics;
    window.dispatchEvent(new CustomEvent('pmx-multiomics-diagnostics', { detail: diagnostics }));
  });

  const details = panel.locator('details.scientific-extra');
  await expect(details).toContainText('Un en-tête importé contient');
  await details.locator('summary').first().click();
  await expect(details).toContainText('<img src=x onerror=');
  await expect(details).toContainText('< 5 observations');
  expect(await page.evaluate(() => ({
    injected: window.__PMX_SCIENCE_INJECTED__,
    injectedHtml: document.querySelectorAll('[data-testid="multiomics-scientific-summary"] img').length
  }))).toEqual({ injected: 0, injectedHtml: 0 });
  await expect(panel).toContainText('Reproduire et rapporter');
});

test('Novice scientific summary uses full precision BH values and never counts Reactome pathways', async ({ page }) => {
  // Serve a controlled DOM fixture, not a mutable Svelte-managed subtree.
  // Exercise the production static browser script against its actual markup.
  await page.route(/\/multiomics\/tool\/?(?:\?.*)?$/, async (route) => {
    await route.fulfill({ status: 200, contentType: 'text/html', body: `<!doctype html>
      <html lang="fr"><head><meta charset="utf-8">
      <script defer src="/multiomics/scientific-summary-ui.js"></script></head>
      <body><main><div data-testid="multiomics-results"><article>
        <div class="layer-title"><span class="method-tag">RNA precision audit</span></div>
        <div class="feature-table" data-testid="multiomics-feature-inference">
          <div class="feature-head"><b>Variable</b><b>Fold ratio</b><b>p</b><b>q BH</b></div>
          <div><code>G1</code><span>2.00×</span><span>0.010</span><span data-scientific-q="0.05004">0.0500</span></div>
          <div><code>G2</code><span>0.50×</span><span>0.010</span><span data-scientific-q="0.04996">0.0500</span></div>
        </div></article>
        <table><thead><tr><th>Reactome pathway</th><th>Assay-universe FDR</th></tr></thead>
          <tbody><tr><td>Unrelated test</td><td>0.00001</td></tr></tbody></table>
      </div></main></body></html>`});
  });
  await page.goto('/multiomics/tool/');
  const panel = page.getByTestId('multiomics-scientific-summary');
  await expect(panel).toBeVisible();
  const evidence = panel.locator('.scientific-summary-grid');
  await expect(evidence).toContainText('RNA precision audit');
  await expect(evidence).toContainText('1/2 lignes affichées');
  await expect(evidence).not.toContainText('Unrelated test');
  await expect(panel).toContainText('q BH');
});
