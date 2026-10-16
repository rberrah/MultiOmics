// @ts-nocheck
import { test, expect } from '@playwright/test';

// Contrast is checked against each element's actual painted parent surface,
// in both explicit light/dark themes. This catches regressions where a CSS
// variable is undefined and falls back to white inside the dark theme.
async function assertContrast(page, textSelector, surfaceSelector, label) {
  const ratio = await page.evaluate(({ textSelector, surfaceSelector }) => {
    const el = document.querySelector(textSelector);
    const surface = document.querySelector(surfaceSelector);
    if (!el || !surface) return null;
    const color = (value) => {
      const numbers = value.match(/[\d.]+/g)?.map(Number);
      if (!numbers || numbers.length < 3) return null;
      return numbers.slice(0, 3).map(v => {
        const c = v / 255;
        return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
      });
    };
    const lum = (rgb) => rgb ? rgb[0] * 0.2126 + rgb[1] * 0.7152 + rgb[2] * 0.0722 : null;
    const fg = lum(color(getComputedStyle(el).color));
    const bg = lum(color(getComputedStyle(surface).backgroundColor));
    if (fg == null || bg == null) return null;
    return (Math.max(fg, bg) + 0.05) / (Math.min(fg, bg) + 0.05);
  }, { textSelector, surfaceSelector });
  expect(ratio, label + ': CSS colors must be computable').not.toBeNull();
  expect(ratio, label + ': text contrast should satisfy WCAG AA 4.5:1').toBeGreaterThanOrEqual(4.5);
}

for (const theme of ['light', 'dark']) {
  test('Multiomics onboarding and visual headlines pass contrast in ' + theme + ' theme', async ({ page }) => {
    await page.emulateMedia({ colorScheme: theme });
    await page.goto('/multiomics/tool');
    await page.evaluate(mode => document.documentElement.setAttribute('data-theme', mode), theme);
    await expect(page.getByTestId('multiomics-quick-start')).toBeVisible();
    await assertContrast(page, '.quick-start strong', '.quick-start', theme + ' quick-start heading');
    await assertContrast(page, '.quick-start span', '.quick-start', theme + ' quick-start description');
    await assertContrast(page, '.quick-start small', '.quick-start', theme + ' quick-start synthetic-data notice');
    for (const key of ['rna', 'protein', 'metabolite'])
      await expect(page.getByTestId('multiomics-quick-demo-' + key)).toBeVisible();

    await page.getByTestId('multiomics-quick-demo').click();
    await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 25000 });
    await expect(page.getByTestId('multiomics-visualizations')).toBeVisible();
    await assertContrast(page, '.visual-head h2', '.visual-panel', theme + ' visualization headline');
    await assertContrast(page, '.visual-head p:not(.eyebrow)', '.visual-panel', theme + ' visualization description');
    await assertContrast(page, '.figure-title h3', '.figure-card', theme + ' figure headline');
    await assertContrast(page, '.figure-title p', '.figure-card', theme + ' figure description');
  });
}

test('The beginner RNA shortcut runs a genuine one-omic analysis and labels it', async ({ page }) => {
  await page.goto('/multiomics/tool');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 25000 });
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText('Omique unique');
  await expect(page.getByTestId('multiomics-results')).toContainText('Résultats de la démonstration');
  await expect(page.getByTestId('multiomics-single-omic-notice')).toBeVisible();
});


test('Scientific evidence status is explicit and external data transfer is opt-in', async ({ page }) => {
  await page.goto('/multiomics/tool');
  await page.getByTestId('multiomics-run-options').locator('summary').first().click();
  await expect(page.getByLabel(/Interroger Reactome|Query Reactome/)).not.toBeChecked();
  await expect(page.getByLabel(/Résoudre les métabolites non canoniques|Resolve selected non-canonical/)).not.toBeChecked();
  await expect(page.getByLabel(/Mode backend R|R backend mode/i)).toHaveValue('browser');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-scientific-assurance')).toBeVisible({ timeout: 25000 });
  await expect(page.getByTestId('multiomics-scientific-assurance')).toContainText(/validation scientifique|scientific validation/i);
  const evidence = await page.getByTestId('multiomics-scientific-assurance').innerText();
  expect(evidence).toMatch(/Exécution réussie|Successful execution/);
  const reportPromise = page.waitForEvent('download');
  await page.getByRole('button', { name: /Télécharger le rapport complet|Download full report/ }).click();
  const report = await reportPromise;
  expect(report.suggestedFilename()).toBe('multiomics_reproducible_report.html');
});
