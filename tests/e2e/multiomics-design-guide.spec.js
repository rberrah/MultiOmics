// @ts-nocheck
import { test, expect } from '@playwright/test';

test('multi-omics design choices are translated into licence-level explanations', async ({ page }) => {
  await page.goto('/multiomics/tool');

  const guide = page.getByTestId('multiomics-design-guide');
  await expect(guide).toBeVisible();
  await expect(guide).toContainText('Traduire le plan d’étude en langage simple');
  await expect(guide).toContainText('Unités indépendantes');
  await expect(guide).toContainText('un réplicat technique n’est pas un nouveau sujet biologique');

  const scope=page.getByTestId('multiomics-design-evidence-scope');
  await expect(scope).toContainText('Deux groupes de sujets distincts');
  await page.getByTestId('multiomics-simple-design').selectOption('paired');
  await expect(scope).toContainText('un test indépendant compterait deux fois');
  await page.getByTestId('multiomics-simple-design').selectOption('crossover');
  await expect(scope).toContainText('Modèle non disponible');
  await page.getByTestId('multiomics-simple-design').selectOption('repeated');
  await expect(guide).toContainText('Mesures répétées');
  await expect(scope).toContainText('Les visites d’un même sujet ne sont pas indépendantes');
  await expect(guide).toContainText('Le même sujet est mesuré plusieurs fois');

  await page.getByLabel(/Objectif principal|Main objective/).selectOption('outcome');
  await page.getByTestId('multiomics-simple-outcome').selectOption('survival');
  await expect(guide).toContainText('il faut une durée de suivi ET');

  await page.locator('.language-toggle').getByRole('button', { name: 'EN', exact: true }).click();
  await expect(guide).toContainText('Translate the study design into plain language');
  await expect(guide).toContainText('a technical replicate is not a new biological subject');
  await expect(scope).toContainText('Repeated visits from one subject are not independent');
});
