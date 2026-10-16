// @ts-nocheck
import { test, expect } from '@playwright/test';

// The statistical warning must survive language changes without asserting
// that a finite Monte Carlo reference benchmark certifies a user study.
test('novice RNA-seq calibration explanation is precise in French and English', async ({page}) => {
  await page.goto('/multiomics/tool?lang=fr');
  const validation = page.getByTestId('multiomics-validation-evidence');
  await expect(validation).toBeVisible();
  await validation.locator('details > summary').first().click();
  let text = page.getByTestId('rnaseq-global-null-calibration-explanation');
  await expect(text).toBeVisible();
  await expect(text).toContainText('12 études standard et 8 études');
  await expect(text).toContainText('26,5 %');
  await expect(text).toContainText('ne prouve pas un contrôle du FDR à 5 %');
  await page.locator('.language-toggle').getByRole('button', {name:'EN',exact:true}).click();
  await validation.locator('details > summary').first().click();
  text = page.getByTestId('rnaseq-global-null-calibration-explanation');
  await expect(text).toBeVisible();
  await expect(text).toContainText('12 standard and 8 studies');
  await expect(text).toContainText('26.5%');
  await expect(text).toContainText('not evidence that FDR is controlled at 5%');
});
