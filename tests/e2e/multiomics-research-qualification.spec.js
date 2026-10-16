// @ts-nocheck
import {test,expect} from '@playwright/test';
import fs from 'node:fs/promises';

test('novice can export reproducible, private methods audit without raw patient rows',async({page})=>{
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated
  ),{timeout:12000}).toBe('true');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await page.getByTestId('multiomics-results').locator('.export-detail summary').click();
  const [download]=await Promise.all([
    page.waitForEvent('download'),
    page.getByTestId('multiomics-methods-audit-export').click()
  ]);
  expect(download.suggestedFilename()).toBe('multiomics_methods_scientific_audit.json');
  const audit=JSON.parse(await fs.readFile(await download.path(),'utf8'));
  expect(audit.format).toBe('pmx-omics-methods-audit/v1');
  expect(audit.scientific_certification).toBe(false);
  expect(audit.clinical_validation).toBe(false);
  expect(audit.omics).toHaveProperty('transcriptomics');
  expect(audit.omics.transcriptomics).toHaveProperty('reference_executed');
  expect(audit.reproducibility_limits.length).toBeGreaterThanOrEqual(4);
  expect(JSON.stringify(audit)).not.toContain('archivedResult');
  expect(JSON.stringify(audit)).not.toContain('rawMatrix');
});

test('novice can find MS technical quality warning without extra mandatory inputs',async({page})=>{
  await page.goto('/multiomics/tool?lang=fr');
  const panel=page.getByTestId('metabolomics-identification-confidence');
  await panel.locator('details.msi-expert > summary').click();
  let boundary=page.getByTestId('multiomics-ms-qc-evidence-boundary');
  await expect(boundary).toContainText('Sans injections QC ou blancs');
  await expect(boundary).toContainText('une variable peut être absente de ces injections');
  await page.locator('.language-toggle').getByRole('button',{name:'EN',exact:true}).click();
  await panel.locator('details.msi-expert > summary').click();
  boundary=page.getByTestId('multiomics-ms-qc-evidence-boundary');
  await expect(boundary).toContainText('If QC/blank injections are absent');
  await expect(boundary).toContainText('a given molecular feature can be absent from those controls');
});
