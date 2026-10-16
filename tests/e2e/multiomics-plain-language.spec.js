// @ts-nocheck
import { test, expect } from '@playwright/test';
import fs from 'node:fs/promises';

test('multi-omics French import screen avoids unexplained metadata jargon', async ({ page }) => {
  await page.goto('/multiomics/tool');

  await expect(page.getByTestId('multiomics-sample-sheet-guide')).toContainText('tableau des échantillons');
  const glossary = page.getByTestId('multiomics-quick-glossary');
  await expect(glossary).toBeVisible();
  await expect(glossary.getByText('Lexique express · six mots utilisés dans tout l’outil')).toBeVisible();
  await glossary.locator('summary').click();
  await expect(glossary).toContainText('Série technique');
  await expect(glossary).toContainText('Facteur d’ajustement');
  await expect(glossary).toContainText('Résultat corrigé (q/FDR)');
  await expect(glossary).toContainText('Validation croisée');

  const validation = page.getByTestId('multiomics-validation-evidence');
  await expect(validation).toBeVisible();
  await expect(validation).toContainText('Validation du logiciel');
  await expect(validation).toContainText('Signaux publics connus');
  await expect(validation).toContainText('Contrôles négatifs');
  await expect(validation).toContainText('Provenance reproductible');
  await expect(validation).toContainText('Ce que cela ne prouve pas');
  await expect(validation).toContainText('93 %');
  await expect(validation).toContainText('58,6 %');
  await expect(validation).toContainText('MNAR');
  await expect(validation).toContainText('200 patientes MAINZ');
  await expect(validation).toContainText('198 patientes TRANSBIG');
  await expect(validation).toContainText('0,8954');
  await expect(validation).toContainText('337 patientes Agilent NKI');
  await expect(validation).toContainText('0,9855');
  await expect(validation).toContainText('0,2793');
  await expect(validation).toContainText('0,2793 → 0,3021');
  await expect(validation).toContainText('nouvelle validation confirmatoire');
  await expect(validation).toContainText('40 souris Nutrimouse');
  await expect(validation).toContainText('(BER) de 0,05');
  await expect(validation).toContainText('Aucun gain prédictif lié à la combinaison des omiques n’est démontré.');
  await expect(validation).toContainText('Les 120 gènes avaient été présélectionnés');
  await expect(validation).toContainText('BER 0,05 pour l’ARN seul');
  await expect(validation).toContainText('Une AUC élevée ne garantit pas des probabilités justes.');

  const metaboliteConfidence = page.getByTestId('metabolomics-identification-confidence');
  await expect(metaboliteConfidence).toBeVisible();
  await expect(metaboliteConfidence).toContainText('Confiance d’identification des métabolites');
  await expect(metaboliteConfidence).toContainText('MSI 1');
  await expect(metaboliteConfidence).toContainText('MSI 4');
  await expect(metaboliteConfidence).toContainText('séparé du mapping ChEBI/HMDB/KEGG');

  await page.locator('.language-toggle').getByRole('button', { name: 'EN', exact: true }).click();
  await expect(glossary).toContainText('Quick glossary · six terms used throughout the tool');
  await expect(glossary).toContainText('Technical series');
  await expect(glossary).toContainText('Adjustment factor');
  await expect(validation).toContainText('Software validation');
  await expect(validation).toContainText('58.6%');
  await expect(validation).toContainText('200 MAINZ participants');
  await expect(validation).toContainText('198 different-study TRANSBIG participants');
  await expect(validation).toContainText('0.8954');
  await expect(validation).toContainText('337 NKI Agilent participants');
  await expect(validation).toContainText('0.9855');
  await expect(validation).toContainText('0.2793');
  await expect(validation).toContainText('0.2793 → 0.3021');
  await expect(validation).toContainText('not fresh confirmatory validation');
  await expect(validation).toContainText('40 Nutrimouse animals');
  await expect(validation).toContainText('BER) of 0.05');
  await expect(validation).toContainText('No added predictive value from combining omics was demonstrated.');
  await expect(validation).toContainText('120 genes were previously selected');
  await expect(validation).toContainText('BER 0.05 for RNA alone');
  await expect(validation).toContainText('A high AUC does not mean trustworthy individual probabilities.');
  await expect(metaboliteConfidence).toContainText('Metabolite identification confidence');
  await page.locator('.language-toggle').getByRole('button', { name: 'FR', exact: true }).click();
  await expect(glossary).toContainText('Lexique express · six mots utilisés dans tout l’outil');

  await expect(page.getByText('Aucun fichier sélectionné').first()).toBeVisible();
  await expect(page.getByText('Détail des colonnes du tableau des échantillons')).toBeVisible();
  await expect(page.getByText('Tableau des échantillons').first()).toBeVisible();
  await expect(page.getByText('Métadonnées échantillons')).toHaveCount(0);
});

test('guided mode lets a licence-level user choose a biological question instead of a method name', async ({ page }) => {
  await page.goto('/multiomics/tool');

  const guide = page.getByTestId('multiomics-beginner-guide');
  await expect(guide).toBeVisible();
  await expect(guide).toContainText('Mode guidé · niveau licence');
  await expect(guide).toContainText('Commencez par la question biologique');
  await expect(guide).toContainText('Qu’est-ce qui varie ensemble entre mes omiques ?');
  await expect(guide).toContainText('Qu’est-ce qui diffère entre mes groupes ?');
  await expect(guide).toContainText('Qu’est-ce qui est associé à mon critère clinique ou expérimental ?');
  await expect(guide).toContainText('Qu’est-ce qui change au cours du temps ?');

  await guide.locator('summary').first().click();
  await guide.getByRole('button', { name: /Qu’est-ce qui diffère entre mes groupes/ }).click();
  await expect(page.getByLabel(/Objectif principal|Main objective/)).toHaveValue('groups');
  await expect(page.getByTestId('multiomics-guided-plan')).toContainText('analyse d’abord chaque omique');
  await expect(guide.getByText('Détails scientifiques de la route choisie')).toBeVisible();

  await guide.getByRole('button', { name: /Qu’est-ce qui est associé à mon critère clinique ou expérimental/ }).click();
  await expect(page.getByLabel(/Objectif principal|Main objective/)).toHaveValue('outcome');
  await expect(page.getByTestId('multiomics-guided-plan')).toContainText('choisit le modèle adapté au type de critère');
});

test('multi-omics result summary uses plain-language study checks and reproducible scientific synthesis', async ({ page }) => {
  await page.route('https://reactome.org/AnalysisService/**', async (route) => {
    await route.fulfill({
      status: 200,
      contentType: 'application/json',
      body: JSON.stringify({
        summary: { token: 'PLAIN_LANGUAGE_TEST', type: 'OVERREPRESENTATION' },
        pathwaysFound: 0,
        identifiersNotFound: 0,
        pathways: []
      })
    });
  });

  await page.goto('/multiomics/tool');
  await page.getByTestId('metabolomics-identification-confidence').locator('summary').first().click();
  await page.locator('#pmx-metabolomics-msi-level').selectOption('msi3');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 20_000 });
  await page.getByTestId('multiomics-run-options').locator('summary').first().click();
  await page.getByLabel(/Interroger Reactome|Query Reactome/i).check();
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:20000});

  const readiness = page.getByTestId('multiomics-plain-readiness');
  await expect(readiness).toBeVisible({ timeout: 10_000 });
  await expect(readiness).toContainText('Vérification avant d’interpréter les résultats');
  await expect(readiness).toContainText('Chevauchement entre omiques');
  await expect(readiness).toContainText('Séries techniques');
  await expect(readiness).toContainText('Données manquantes et mesures à vérifier');
  await expect(readiness.getByText('À propos de la puissance statistique')).toBeVisible();

  const readingGuide = page.getByTestId('multiomics-result-reading-guide');
  await expect(readingGuide).toBeVisible();
  await expect(readingGuide).toContainText('Comprendre le résultat en quatre passages');
  await expect(readingGuide).toContainText('Qualité d’abord');
  await expect(readingGuide).toContainText('Dans chaque omique');
  await expect(readingGuide).toContainText('Entre les omiques');
  await expect(readingGuide).toContainText('Interprétation biologique');
  await expect(readingGuide).toContainText('Ce que vous pouvez conclure');
  await expect(readingGuide).toContainText('Ce que vous ne devez pas conclure');

  const synthesis = page.getByTestId('multiomics-scientific-summary');
  await expect(synthesis).toBeVisible({ timeout: 10_000 });
  await expect(synthesis).toContainText('Synthèse scientifique');
  await expect(synthesis).toContainText('Que permettent réellement de conclure ces résultats ?');
  await expect(synthesis).toContainText('Association ≠ prédiction ≠ causalité');
  await expect(synthesis).toContainText(/q\/FDR/i);

  const fdrPolicy = page.getByTestId('multiomics-fdr-policy');
  await expect(fdrPolicy).toBeVisible({ timeout: 10_000 });
  await expect(fdrPolicy).toContainText('Deux niveaux de FDR');
  await expect(fdrPolicy).toContainText('Preuve principale · q ≤ 0,05');
  await expect(fdrPolicy).toContainText('Signal exploratoire · 0,05 < q ≤ 0,10');
  await expect(fdrPolicy).toContainText('p-values brutes');
  await expect(fdrPolicy).toContainText('q-value n’est pas la probabilité');

  const analysisDownloadPromise = page.waitForEvent('download');
  await page.locator('.export-detail > summary').click();
  await page.getByRole('button', { name: /Données détaillées \(JSON\)/ }).click();
  const analysisDownload = await analysisDownloadPromise;
  const analysisPath = await analysisDownload.path();
  expect(analysisPath).toBeTruthy();
  const analysisJson = JSON.parse(await fs.readFile(analysisPath, 'utf8'));
  expect(analysisJson.protocol.metabolomicsIdentificationConfidence).toBe('msi3');
  expect(analysisJson.metabolomicsAnnotationConfidence.scheme).toContain('Metabolomics Standards Initiative');
  expect(analysisJson.metabolomicsAnnotationConfidence.exactPathwayMappingSuppressed).toBe(true);
  expect(analysisJson.reactome.metabolomicsConfidenceSuppressed).toBe(true);
  expect(analysisJson.reactome.combined).toMatchObject({
    token: null,
    pathwaysFound: 0,
    pathways: [],
    suppressed: true
  });
  expect(analysisJson.reactome.consensus).toEqual([]);

  await page.getByTestId('multiomics-scientific-summary').locator('details.scientific-extra > summary').click();
  const manifest = page.getByTestId('multiomics-methods-manifest');
  await expect(manifest).toBeVisible();
  await expect(manifest).toContainText('Reproduire et rapporter l’analyse');
  await expect(manifest.getByRole('button', { name: 'Rapport méthodes (.md)' })).toBeVisible();
  await expect(manifest.getByRole('button', { name: 'Manifeste (.json)' })).toBeVisible();

  const jsonDownload = page.waitForEvent('download');
  await manifest.getByRole('button', { name: 'Manifeste (.json)' }).click();
  const downloadedManifest = await jsonDownload;
  expect(downloadedManifest.suggestedFilename()).toBe('multiomics-reproducibility-manifest.json');

  await page.getByTestId('multiomics-quality-details').locator('summary').first().click();
  await expect(page.getByText('Vérification des séries techniques')).toBeVisible();
  await expect(page.getByText('Audit des batches techniques')).toHaveCount(0);

  const interpretation = page.getByTestId('multiomics-interpretation');
  await expect(interpretation).toContainText('Différence observée');
  await expect(interpretation).toContainText('Résultat après correction statistique');
  await expect(interpretation).not.toContainText('Effet biologique');
  await expect(interpretation).not.toContainText('Significativité');
});
