// @ts-nocheck
import {test,expect} from '@playwright/test';
import fs from 'node:fs/promises';

test('Portable project exports exact input files and reopens locally without silent external requests',async({page})=>{
  await page.goto('/multiomics/tool');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await page.getByTestId('multiomics-results').locator('.export-detail summary').click();
  const [download]=await Promise.all([
    page.waitForEvent('download'),
    page.getByTestId('multiomics-project-export').click()
  ]);
  expect(download.suggestedFilename()).toBe('multiomics_project_reproducible_v1.json');
  const downloadedPath=await download.path();
  const raw=JSON.parse(await fs.readFile(downloadedPath,'utf8'));
  expect(raw.format).toBe('pmx-multiomics-portable-project');
  expect(raw.inputs.metadata.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(raw.inputs.transcriptomics.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(raw.inputs.metadata.content).toContain('subject_id');
  expect(raw.archivedResult).toHaveProperty('layers');
  expect(raw.settings.demoLoaded).toBe(true);
  expect(JSON.stringify(raw.settings)).not.toContain('referenceBackendUrl');

  await page.reload();
  await page.getByTestId('multiomics-project-import-panel').locator('summary').click();
  let externalConnections=0;
  page.on('request',request=>{
    if(/reactome\.org|ebi\.ac\.uk|\/run(\?|$)/.test(request.url()))externalConnections++;
  });
  await page.getByTestId('multiomics-project-import-file').setInputFiles(downloadedPath);
  await expect(page.getByRole('status').filter({hasText:'SHA-256'})).toBeVisible({timeout:15000});
  await expect(page.getByTestId('multiomics-results')).toHaveCount(0);
  await expect(page.getByTestId('multiomics-run')).toBeEnabled();
  expect(externalConnections).toBe(0);
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText('Omique unique');
  expect(externalConnections).toBe(0);
});

test('Tampered project input is rejected instead of silently restoring modified patient data',async({page})=>{
  await page.goto('/multiomics/tool');
  await page.getByTestId('multiomics-project-import-panel').locator('summary').click();
  const payload={
    format:'pmx-multiomics-portable-project',schemaVersion:1,
    inputs:{metadata:{
      filename:'modified.csv',mime:'text/csv',
      sizeBytes:8,sha256:'0'.repeat(64),
      content:'a,b\n1,2\n'
    }},
    settings:{analysisIntent:'exploratory'}
  };
  await page.getByTestId('multiomics-project-import-file').setInputFiles({
    name:'corrupted.json',mimeType:'application/json',
    buffer:Buffer.from(JSON.stringify(payload))
  });
  await expect(page.getByRole('alert')).toContainText('SHA-256');
  await expect(page.getByTestId('multiomics-results')).toHaveCount(0);
});
