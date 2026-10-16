// @ts-nocheck
import { test, expect } from '@playwright/test';

test('novice bilingual RNA-seq mixed truth explanation distinguishes FDR from power', async ({page}) => {
  await page.goto('/multiomics/tool?lang=fr');
  const validation = page.getByTestId('multiomics-validation-evidence');
  await expect(validation).toBeVisible();
  await validation.locator('details > summary').first().click();
  let warning = page.getByTestId('rnaseq-mixed-truth-calibration-explanation');
  await expect(warning).toBeVisible();
  await expect(warning).toContainText('1 %, 5 % ou 20 %');
  await expect(warning).toContainText('faux positifs');
  await expect(warning).toContainText('puissance');
  await expect(warning).toContainText('11,1 %');
  await expect(warning).toContainText('sans permettre de classer fiablement');
  await expect(warning).toContainText("Quelques simulations artificielles ne valident pas");
  await page.locator('.language-toggle').getByRole('button', {name:'EN',exact:true}).click();
  await validation.locator('details > summary').first().click();
  warning = page.getByTestId('rnaseq-mixed-truth-calibration-explanation');
  await expect(warning).toBeVisible();
  await expect(warning).toContainText('1%, 5% or 20%');
  await expect(warning).toContainText('false positives');
  await expect(warning).toContainText('power');
  await expect(warning).toContainText('11.1%');
  await expect(warning).toContainText('not a reliable comparison');
  await expect(warning).toContainText('cannot validate your biological study');
});
