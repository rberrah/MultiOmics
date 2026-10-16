// @ts-nocheck
import {test,expect} from '@playwright/test';

test('novice RNA-seq filtering warning preserves scientific context FR/EN', async ({page}) => {
  await page.goto('/multiomics/tool?lang=fr');
  const validation=page.getByTestId('multiomics-validation-evidence');
  await expect(validation).toBeVisible();
  await validation.locator('details > summary').first().click();
  let info=page.getByTestId('rnaseq-prefilter-sensitivity-guide');
  await expect(info).toBeVisible();
  await expect(info).toContainText('20 comptages dans 40 % des échantillons');
  await expect(info).toContainText('recalculent la normalisation');
  await expect(info).toContainText("la chaîne d'analyse complète");
  await expect(info).toContainText('ne faut pas choisir après coup');
  await expect(info).toContainText('190 gènes sur 400');
  await expect(info).toContainText('9,6 des 20 vrais signaux');
  await expect(info).toContainText('ne valide pas automatiquement votre étude');
  await page.locator('.language-toggle').getByRole('button',
    {name:'EN',exact:true}).click();
  await validation.locator('details > summary').first().click();
  info=page.getByTestId('rnaseq-prefilter-sensitivity-guide');
  await expect(info).toBeVisible();
  await expect(info).toContainText('20 counts in 40% of samples');
  await expect(info).toContainText('re-estimate normalization');
  await expect(info).toContainText('full analysis pipeline');
  await expect(info).toContainText('Do not select whichever filter');
  await expect(info).toContainText('190 of 400 genes');
  await expect(info).toContainText('9.6 of the 20 planted true signals');
  await expect(info).toContainText('not validation of your own study');
});
