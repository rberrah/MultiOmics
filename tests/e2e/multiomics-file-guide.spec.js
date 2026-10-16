// @ts-nocheck
import { test, expect } from '@playwright/test';

test('multi-omics upload area explains the minimum file structure to beginners', async ({ page }) => {
  await page.goto('/multiomics/tool');

  const guide = page.getByTestId('multiomics-file-guide');
  await expect(guide).toBeVisible();
  await expect(guide).toContainText('Que faut-il charger ?');
  await expect(guide).toContainText('1 à 3 matrices omiques + tableau des échantillons (ou assistant)');
  await expect(guide).toContainText(/Deux omiques différentes permettent une analyse multi-omique/i);
  await expect(guide).toContainText('subject_id · sample_id · assay_id · omic');
  await expect(guide).toContainText('feature_id | RNA001 | RNA002');
  await expect(guide).toContainText('Comment les fichiers se relient');
  await expect(guide).toContainText('Exemple mono-omique (ARN uniquement)');
  await expect(guide).toContainText('RNA001  transcriptomics');
  await expect(guide).toContainText('l’assistant ne doit pas les inventer');
  await expect(guide).toContainText('échantillon biologique indépendant');

  await page.locator('.language-toggle').getByRole('button', { name: 'EN', exact: true }).click();
  await expect(guide).toContainText('What do I need to upload?');
  await expect(guide).toContainText('1–3 omics matrices + a sample sheet (or guided creation)');
  await expect(guide).toContainText('Single-omic example (RNA only)');
  await expect(guide).toContainText('provide actual sample metadata');
});
