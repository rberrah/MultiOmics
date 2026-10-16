// @ts-nocheck
// Professional/reviewer acceptance: user can complete an actual analysis,
// EXPORT all needed scientific metadata, then restore original inputs and rerun.
// This is reproducibility of a synthetic tutorial, NOT clinical qualification.
import {test,expect} from '@playwright/test';
import fs from 'node:fs/promises';

function validateResult(data,layerCount){
  expect(data).toHaveProperty('engine');
  expect(data).toHaveProperty('protocol');
  const entries=Object.entries(data.layers||{});
  expect(entries.length).toBe(layerCount);
  for(const [name,layer] of entries){
    expect(['transcriptomics','proteomics','metabolomics']).toContain(name);
    expect(layer).toHaveProperty('rows');
    expect(layer.rows.length).toBeGreaterThan(0);
    for(const row of layer.rows){
      expect(typeof row.feature).toBe('string');
      expect(row.feature.length).toBeGreaterThan(0);
      for(const key of ['pValue','qValue']){
        if(row[key]!==null&&row[key]!==undefined){
          expect(typeof row[key]).toBe('number');
          expect(Number.isFinite(row[key])).toBe(true);
          expect(row[key]).toBeGreaterThanOrEqual(0);
          expect(row[key]).toBeLessThanOrEqual(1);
        }
      }
    }
  }
  const marker=JSON.stringify(data).toLowerCase();
  expect(marker).not.toMatch(/"scientific_certification"\s*:\s*true/);
  expect(marker).not.toMatch(/"clinical_validation"\s*:\s*true/);
}
/** @param {import('@playwright/test').Page} page */
async function getJsonExport(page,button){
  if(!(await page.locator('.export-detail').evaluate(el=>el.open)))
    await page.locator('.export-detail > summary').click();
  const [download]=await Promise.all([
    page.waitForEvent('download',{timeout:15000}),
    button(page).click()
  ]);
  return JSON.parse(await fs.readFile(await download.path(),'utf8'));
}
const detailButton=page=>page.getByRole('button',{name:/Données détaillées|Detailed data/i});
const methodsButton=page=>page.getByTestId('multiomics-methods-audit-export');
const projectButton=page=>page.getByTestId('multiomics-project-export');

test('Scientific research export: a complete one-omic study reopens and recomputes the same numeric results',async({page})=>{
  test.setTimeout(90000);
  const errs=[];
  page.on('pageerror',e=>errs.push(e.message));
  let unapproved=0;
  page.on('request',req=>{
    if(/reactome\.org|ebi\.ac\.uk|127\.0\.0\.1:8787\/run/i.test(req.url()))unapproved++;
  });
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const initial=await getJsonExport(page,detailButton);
  validateResult(initial,1);
  const audit=await getJsonExport(page,methodsButton);
  expect(audit.format).toBe('pmx-omics-methods-audit/v1');
  expect(audit.scientific_certification).toBe(false);
  expect(audit.clinical_validation).toBe(false);
  expect(audit.omics.transcriptomics.reference_executed).toBe(false);
  expect(audit.reference_runtime_summary).toBeNull();
  const portable=await getJsonExport(page,projectButton);
  expect(portable.format).toBe('pmx-multiomics-portable-project');
  expect(portable.inputs.metadata.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(portable.inputs.transcriptomics.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(portable.privacy).toMatch(/SENSITIVE/);
  await page.reload();
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  await page.getByTestId('multiomics-project-import-panel').locator('summary').click();
  await page.getByTestId('multiomics-project-import-file').setInputFiles({
    name:'scientific-project.json',mimeType:'application/json',
    buffer:Buffer.from(JSON.stringify(portable),'utf8')
  });
  await expect(page.getByRole('status').filter({hasText:'SHA-256'}))
    .toBeVisible({timeout:20000});
  await expect(page.getByTestId('multiomics-results')).toHaveCount(0);
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const rerun=await getJsonExport(page,detailButton);
  validateResult(rerun,1);
  const numericRows=result=>result.layers.transcriptomics.rows.map(
    x=>({feature:x.feature,p:x.pValue,q:x.qValue,effect:x.effect})
  );
  expect(numericRows(rerun)).toEqual(numericRows(initial));
  expect(unapproved).toBe(0);
  expect(errs).toEqual([]);
});

test('Scientific researcher workflows: each single-omic example produces valid, downloadable data without a hidden reference engine',async({page})=>{
  test.setTimeout(100000);
  await page.goto('/multiomics/tool?lang=en');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  for(const [id,layer] of [
    ['multiomics-quick-demo-rna','transcriptomics'],
    ['multiomics-quick-demo-protein','proteomics'],
    ['multiomics-quick-demo-metabolite','metabolomics']
  ]){
    await page.getByTestId(id).click();
    await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
    const details=await getJsonExport(page,detailButton);
    validateResult(details,1);
    expect(Object.keys(details.layers)).toEqual([layer]);
    const audit=await getJsonExport(page,methodsButton);
    expect(audit.scientific_certification).toBe(false);
    expect(audit.omics).toHaveProperty(layer);
    expect(audit.omics[layer].reference_executed).toBe(false);
    expect(audit.reproducibility_limits.length).toBeGreaterThanOrEqual(4);
  }
});

test('Scientific researcher workflows: three-omics demo produces genuine individual layer tables and documents integration limits',async({page})=>{
  test.setTimeout(70000);
  await page.goto('/multiomics/tool?lang=en');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  await page.getByTestId('multiomics-quick-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const details=await getJsonExport(page,detailButton);
  validateResult(details,3);
  const audit=await getJsonExport(page,methodsButton);
  expect(Object.keys(audit.omics).sort()).toEqual(
    ['metabolomics','proteomics','transcriptomics']);
  expect(audit.clinical_validation).toBe(false);
  expect(audit.scientific_certification).toBe(false);
  for(const layer of Object.values(audit.omics))
    expect(layer.reference_executed).toBe(false);
});


test('Researcher supplies an original custom RNA matrix; source identities are retained and statistics are exported without invented groups',async({page})=>{
  test.setTimeout(65000);
  await page.goto('/multiomics/tool?lang=fr');
  // On public/static builds, SSR inputs are visible BEFORE Svelte handlers attach.
  // A true researcher must upload only after hydration, otherwise the file is
  // selected in the native input but the analysis state is not updated.
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  await page.locator('.uploads input[type="file"]').nth(1).setInputFiles({
    name:'external-lab-test-matrix.csv',mimeType:'text/csv',
    buffer:Buffer.from([
      'feature_id,S01,S02,S03,S04,S05,S06,S07,S08',
      'GENE1,10,12,15,20,18,17,14,13',
      'GENE2,30,31,35,39,45,49,41,37',
      'GENE3,55,58,61,63,60,57,56,53',
      'GENE4,8,11,10,9,14,16,17,13'
    ].join('\n'))
  });
  await expect(page.getByTestId('multiomics-single-sheet-helper')).toBeVisible();
  const helper=page.getByTestId('multiomics-single-sheet-helper');
  await helper.getByRole('checkbox').check();
  await page.getByTestId('multiomics-generate-explore-metadata').click();
  await expect(page.getByTestId('multiomics-generated-sheet-notice'))
    .toContainText(/aucun groupe inféré|no groups inferred/);
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const data=await getJsonExport(page,detailButton);
  validateResult(data,1);
  const manifest=await getJsonExport(page,methodsButton);
  expect(manifest.scientific_certification).toBe(false);
  expect(manifest.experimental_design).not.toBe('externally_validated');
  const project=await getJsonExport(page,projectButton);
  expect(project.inputs.transcriptomics.filename).toBe('external-lab-test-matrix.csv');
  expect(project.inputs.transcriptomics.content).toContain('GENE4');
  expect(project.inputs.transcriptomics.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(project.inputs.metadata.content).toContain('subject_id');
});
