// @ts-nocheck
import { test, expect } from '@playwright/test';

test('RNA, proteomics and metabolomics can each be analysed without fake second layer', async ({ page }) => {
  await page.goto('/multiomics/tool');
  for (const [layer,button] of [
    ['transcriptomics','multiomics-demo-rna'],
    ['proteomics','multiomics-demo-protein'],
    ['metabolomics','multiomics-demo-metabolite']
  ]) {
    await page.getByTestId(button).click();
    const results=page.getByTestId('multiomics-results');
    await expect(results).toBeVisible({timeout:35_000});
    await expect(page.getByTestId('multiomics-analysis-mode')).toContainText(/Omique unique|Single omic/);
    await expect(page.getByTestId('multiomics-single-omic-notice')).toBeVisible();
    await expect(results).toContainText(layer==='transcriptomics' ? /Transcriptomique|Transcriptomics/
      : layer==='proteomics' ? /Protéomique|Proteomics/ : /Métabolomique|Metabolomics/);
    await expect(page.getByTestId('multiomics-quality-summary')).toContainText(/Méthodes réellement exécutées|Methods actually run/);
    await page.getByLabel(/Objectif principal|Main objective/).selectOption('explore');
    await page.getByTestId('multiomics-run').click();
    await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:35_000});
    await expect(page.getByTestId('multiomics-analysis-mode')).toContainText(/Omique unique|Single omic/);
    // Re-enable a time-course demo for the next iteration.
    await page.getByLabel(/Objectif principal|Main objective/).selectOption('time');
  }
});

test('MS import preserves vendor annotation fields and distinguishes features from identified compounds', async ({page})=>{
  await page.goto('/multiomics/tool');
  await page.waitForLoadState('networkidle');
  const rows=['feature_id,assay_id,auc,mz,rt,adduct,condition,subject_id,sample_id'];
  for(let i=0;i<8;i++){
    for(const [feature,mz,rt,adduct,ctrl,treated] of [
      ['MS_001','123.45','2.35','[M+H]+',100,450],
      ['MS_002','567.89','4.56','[M-H]-',350,80]
    ]) rows.push([feature,'INJ'+(i+1),i<4?ctrl+i:treated+i,mz,rt,adduct,
      i<4?'control':'treated','S'+(i+1),'S'+(i+1)].join(','));
  }
  await page.getByTestId('multiomics-ms-auc-upload').setInputFiles({
    name:'vendor_original.csv',mimeType:'text/csv',buffer:Buffer.from(rows.join('\n'))
  });
  await expect(page.getByTestId('multiomics-ms-auc-import-status')).toContainText('2 signaux MS');
  await expect(page.getByTestId('multiomics-ms-annotations-preserved')).toContainText('mz');
  await expect(page.getByTestId('multiomics-ms-annotations-preserved')).toContainText('rt');
  await expect(page.getByTestId('multiomics-ms-annotations-preserved')).toContainText('adduct');
  await expect(page.getByTestId('multiomics-single-omic-notice')).toBeVisible();
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30_000});
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText('Omique unique');
});

test('A single omics matrix can be explored after explicit independent-sample confirmation', async ({page}) => {
  await page.goto('/multiomics/tool');
  await page.getByRole('heading', {name:/Analyser une ou plusieurs omiques|Analyze one or several omics/}).waitFor();
  await page.locator('.uploads input[type="file"]').nth(1).setInputFiles({
    name:'rna_feature_matrix.csv',mimeType:'text/csv',
    buffer:Buffer.from([
      'feature_id,S01,S02,S03,S04,S05',
      'G1,10,12,16,14,11',
      'G2,90,55,70,40,65',
      'G3,5,8,12,10,9'
    ].join('\n'))
  });
  const helper=page.getByTestId('multiomics-single-sheet-helper');
  await expect(helper).toBeVisible();
  const create=page.getByTestId('multiomics-generate-explore-metadata');
  await expect(create).toBeDisabled();
  await helper.getByRole('checkbox').check();
  await create.click();
  await expect(page.getByTestId('multiomics-generated-sheet-notice')).toContainText(/aucun groupe inféré|no groups inferred/);
  await expect(page.getByTestId('multiomics-run')).toBeEnabled();
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30_000});
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText(/Omique unique|Single omic/);
});
