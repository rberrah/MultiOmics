// @ts-nocheck
import {test,expect} from '@playwright/test';
import fs from 'node:fs/promises';
function fixture(){
  const samples=Array.from({length:8},(_,i)=>'ASSAY_'+String(i+1).padStart(2,'0'));
  const metadata=['subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate,sample_type'];
  for(let i=0;i<8;i++){
    const sid='SUBJECT_'+String(Math.floor(i/2)+1).padStart(2,'0');
    metadata.push([sid,samples[i],samples[i],'transcriptomics',
      i%2?'treated':'control','T0','RUN1','1','biological'].join(','));
  }
  const raw=['feature_id,'+samples.join(',')];
  for(let gene=1;gene<=80;gene++){
    raw.push('ENSG_DEMO_'+String(gene).padStart(4,'0')+','+
      samples.map((_,i)=>Math.round(65+gene%21+i*4+(i%2&&gene<=10?25:0))).join(','));
  }
  return {metadata:metadata.join('\n')+'\n',matrix:raw.join('\n')+'\n'};
}
test('Novice researcher sees four biological pairs, downloads template, and cannot run paired RNA without actual R',async({page})=>{
  test.setTimeout(60_000);
  const errs=[];
  page.on('pageerror',e=>errs.push(e.message));
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:14000}).toBe('true');
  await page.getByLabel(/Objectif principal|Main objective/).selectOption('groups');
  await page.getByTestId('multiomics-simple-design').selectOption('paired');
  await expect(page.getByTestId('multiomics-paired-help')).toContainText('Quatre paires');
  const [download]=await Promise.all([
    page.waitForEvent('download',{timeout:12000}),
    page.getByTestId('multiomics-paired-template-download').click()
  ]);
  expect(download.suggestedFilename()).toBe('EXEMPLE_FICTIF_4_PAIRES_metadata.csv');
  const template=(await fs.readFile(await download.path())).toString('utf8');
  expect(template).toContain('DEMO_SUBJECT_04_B');
  expect(template.split('\n').filter(Boolean).length).toBe(9);
  const files=fixture();
  await page.locator('.uploads input[type="file"]').nth(0).setInputFiles({
    name:'researcher_paired_metadata.csv',mimeType:'text/csv',
    buffer:Buffer.from(files.metadata,'utf8')
  });
  await page.locator('.uploads input[type="file"]').nth(1).setInputFiles({
    name:'researcher_paired_raw_counts.csv',mimeType:'text/csv',
    buffer:Buffer.from(files.matrix,'utf8')
  });
  await expect(page.getByTestId('multiomics-paired-count')).toContainText('4',{timeout:15000});
  await expect(page.getByTestId('multiomics-paired-count')).toContainText('8');
  await page.getByTestId('multiomics-run').click();
  await expect(page.locator('p.error,[role="alert"]').filter({hasText:/R requis|Requires R|moteur R local/})
    .first()).toBeVisible({timeout:10000});
  expect(errs).toEqual([]);
});
