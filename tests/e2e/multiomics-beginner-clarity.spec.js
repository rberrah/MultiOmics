// @ts-nocheck
import { test, expect } from '@playwright/test';

test('FDR explanation counts real feature-grid q-values, not unrelated table q-values', async ({ page }) => {
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 30_000 });

  await page.evaluate(() => {
    const results = document.querySelector('[data-testid="multiomics-results"]');
    const fixture = document.createElement('div');
    fixture.className = 'feature-table';
    fixture.dataset.testid = 'fdr-grid-fixture';
    fixture.innerHTML = '<span data-scientific-q="0.02">0,020</span>' +
      '<span data-scientific-q="0.08">0,080</span>' +
      '<span data-scientific-q="">—</span>' +
      '<span data-scientific-q="NaN">invalid</span>' +
      '<span data-scientific-q="1.5">invalid</span>';
    results.append(fixture);

    const unrelatedTable = document.createElement('table');
    unrelatedTable.dataset.testid = 'unrelated-fdr-table';
    unrelatedTable.innerHTML = '<thead><tr><th>Pathway</th><th>q BH</th></tr></thead>' +
      '<tbody><tr><td>Not a feature</td><td>0.001</td></tr></tbody>';
    results.append(unrelatedTable);
    window.dispatchEvent(new Event('pmx-multiomics-analysis'));
  });

  const policy = page.getByTestId('multiomics-fdr-policy');
  await expect(policy).toBeVisible();
  const actualCountsMatchScreen = async (lang) => page.evaluate((locale) => {
    const results = document.querySelector('[data-testid="multiomics-results"]');
    const raw = [...(results?.querySelectorAll('.feature-table [data-scientific-q]') || [])]
      .map(node => node.getAttribute('data-scientific-q'))
      .filter(value => value !== null && value.trim() !== '');
    const qs = raw.map(Number).filter(q => Number.isFinite(q) && q >= 0 && q <= 1);
    const primary = qs.filter(q => q <= 0.05).length;
    const exploratory = qs.filter(q => q > 0.05 && q <= 0.10).length;
    const text = document.querySelector('[data-testid="multiomics-fdr-policy"] .fdr-policy-count')?.textContent || '';
    const first = locale === 'en'
      ? `${primary} displayed row(s) meet q ≤ 0.05`
      : `${primary} ligne(s) affichée(s) satisfont q ≤ 0,05`;
    const second = locale === 'en'
      ? `${exploratory} additional row(s)`
      : `${exploratory} ligne(s) supplémentaire(s)`;
    return qs.length >= 2 && text.includes(first) && text.includes(second);
  }, lang);

  // The current Svelte view may update between the initial demo and the
  // injected fixture. Always compare against q-values still on screen.
  await expect.poll(() => actualCountsMatchScreen('fr'), { timeout: 10_000 }).toBe(true);
  await expect(policy).toContainText('dix premières lignes visibles par omique');

  await page.locator('.language-toggle').getByRole('button', { name: 'EN', exact: true }).click();
  await expect.poll(() => actualCountsMatchScreen('en'), { timeout: 10_000 }).toBe(true);
  await expect(policy).toContainText('first ten rows shown per omic');
});

test('beginner workflow is browser-first, exploratory and keeps R advanced', async ({ page }) => {
  await page.goto('/multiomics/tool?lang=fr');

  const policy = page.getByTestId('multiomics-exploratory-policy');
  await expect(policy).toBeVisible();
  await expect(policy).toContainText('Données exploratoires');
  await expect(policy).toContainText('code R');

  const advanced = page.getByTestId('multiomics-advanced-analysis-options');
  await expect(advanced).toBeVisible();
  await expect(advanced).not.toHaveAttribute('open', '');

  const run = page.getByTestId('multiomics-run');
  await expect(run).toContainText('Explorer mes données');

  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 30_000 });
  await expect(page.getByTestId('multiomics-research-bundle-download'))
    .toContainText('dossier + code R');
});

test('Readiness text derived from visible metadata is escaped, never injected as HTML', async ({ page }) => {
  await page.goto('/multiomics/tool?lang=fr');
  const readiness = page.getByTestId('multiomics-run-readiness');
  await expect(readiness).toBeVisible();
  await page.evaluate(() => {
    document.querySelector('.status strong').textContent = '<img src=x onerror=alert(1)>';
  });
  await expect(readiness).toContainText('<img src=x onerror=alert(1)>');
  await expect(readiness.locator('img')).toHaveCount(0);
});

test('expert signature upload refuses invalid untraceable scores without affecting default demo', async ({page}) => {
  await page.goto('/multiomics/tool/?lang=en');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30_000});
  await page.getByText('Other downloads', {exact:true}).click();
  const upload=page.getByTestId('multiomics-signature-upload');
  await expect(upload).toBeVisible();
  await upload.setInputFiles({
    name:'unverified_group_effect.csv',mimeType:'text/csv',
    buffer:Buffer.from('subject_id,sample_id,assay_id,omic,signature_id,score\\nP01,S01,RNA1,transcriptomics,group_logfc,2.4\\n')
  });
  await expect(page.getByTestId('multiomics-signature-error')).toContainText('missing required');
  await expect(page.getByTestId('multiomics-omics2pk-readiness')).toContainText('Exchange unavailable');
  await page.getByTestId('multiomics-research-bundle-download').click();
  await expect(page.getByTestId('multiomics-research-bundle-error')).toContainText('Individual signature import is invalid');
});


test('results explain evidence tier and verify Omics2PK structural exchange without claiming PK readiness', async ({page}) => {
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30_000});

  const interpretation=page.getByTestId('multiomics-result-interpretation');
  await expect(interpretation).toBeVisible();
  await expect(interpretation).toContainText('Source numérique');
  await expect(interpretation).toContainText('Exploratoire navigateur');
  await expect(interpretation).toContainText(/méthode de référence|reference method/);

  const reuse=page.getByTestId('multiomics-omics2pk-readiness');
  await expect(reuse).toContainText('Échange individuel disponible');
  await expect(reuse).toContainText('8 sujets');
  await expect(reuse).toContainText('49 assays');
  await expect(reuse).toContainText('PK non prête');
  await expect(reuse).toContainText(/vraies doses|genuine PK doses/);

  const panel=page.getByTestId('multiomics-omics2pk-export-panel');
  await expect(panel).toContainText(/matrices individuelles|Individual matrices/);
  await expect(panel).toContainText(/doses, temps et concentrations|doses, times and concentrations/);

  const downloadPromise=page.waitForEvent('download');
  await page.getByTestId('multiomics-omics2pk-export').click();
  await downloadPromise;
  await expect(page.getByTestId('multiomics-omics2pk-export-status'))
    .toContainText('non prêt pour une analyse PK');
});
