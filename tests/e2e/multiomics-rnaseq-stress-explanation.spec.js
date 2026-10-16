// @ts-nocheck
import { test, expect } from '@playwright/test';

test('beginners get accurate bilingual caveats on RNA-seq low power and small cohorts', async ({page}) => {
  await page.goto('/multiomics/tool?lang=fr');
  const section = page.getByTestId('multiomics-validation-evidence');
  await expect(section).toBeVisible();
  await section.locator('details > summary').first().click();
  let warning = page.getByTestId('rnaseq-stress-grid-explanation');
  await expect(warning).toBeVisible();
  await expect(warning).toContainText('six sujets par groupe');
  await expect(warning).toContainText('log2 fold change de 0,35');
  await expect(warning).toContainText('six patients par groupe suffisent');
  await expect(warning).toContainText('la puissance mesure');
  await expect(warning).toContainText('2,5 %');
  await expect(warning).toContainText('0 %');
  await expect(warning).toContainText('trop incertain pour comparer');
  await page.locator('.language-toggle').getByRole('button', {name:'EN',exact:true}).click();
  await section.locator('details > summary').first().click();
  warning = page.getByTestId('rnaseq-stress-grid-explanation');
  await expect(warning).toBeVisible();
  await expect(warning).toContainText('six subjects per group');
  await expect(warning).toContainText('log2 fold change 0.35');
  await expect(warning).toContainText('six patients per group are enough');
  await expect(warning).toContainText('power describes');
  await expect(warning).toContainText('2.5%');
  await expect(warning).toContainText('0%');
  await expect(warning).toContainText('too small and uncertain to compare');
});
