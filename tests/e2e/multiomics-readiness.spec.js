// @ts-nocheck
import { test, expect } from '@playwright/test';

test('multi-omics pre-run gate explains readiness in licence-level language', async ({ page }) => {
  await page.goto('/multiomics/tool');

  const readiness = page.getByTestId('multiomics-run-readiness');
  await expect(readiness).toBeVisible();
  await expect(readiness).toContainText('L’analyse est-elle prête ?');
  await expect(readiness).toContainText('À compléter');
  await expect(readiness).toContainText('une seule pour l’analyse mono-omique');
  await expect(readiness).toContainText('tableau des échantillons');

  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 20_000 });
  await expect(readiness).toContainText('Prêt pour explorer');
  await expect(readiness).toContainText('Les contrôles d’entrée permettent de lancer le calcul');
  await expect(readiness).toContainText('pas de valider un résultat scientifique');
  await page.getByTestId('multiomics-analysis-intent').locator('select').selectOption('confirmatory');
  await expect(readiness).toContainText('Prêt pour les contrôles R');
  await expect(readiness).not.toContainText('Prêt à analyser');
});
