// @ts-nocheck
import {test,expect} from '@playwright/test';

test('beginner documentation rejects zero-discovery equals certified FDR, in both languages',async ({page})=>{
  await page.goto('/multiomics/tool?lang=fr');
  const evidence=page.getByTestId('multiomics-validation-evidence');
  await expect(evidence).toBeVisible();
  await evidence.locator('details > summary').first().click();
  let note=page.getByTestId('rnaseq-fdr-small-sample-guard');
  await expect(note).toBeVisible();
  await expect(note).toContainText('huit études simulées');
  await expect(note).toContainText('0 % à 0 %');
  await expect(note).toContainText('nombre d\'études indépendantes est insuffisant');
  await expect(note).toContainText('Zéro erreur observée ne signifie pas validation scientifique');
  await page.locator('.language-toggle').getByRole('button',{name:'EN',exact:true}).click();
  await evidence.locator('details > summary').first().click();
  note=page.getByTestId('rnaseq-fdr-small-sample-guard');
  await expect(note).toBeVisible();
  await expect(note).toContainText('eight simulated studies');
  await expect(note).toContainText('0% to 0%');
  await expect(note).toContainText('across all scenario-and-method comparisons');
  await expect(note).toContainText('Zero observed errors is not scientific certification');
});
