// @ts-nocheck
import {test,expect} from '@playwright/test';

test('confirmatory R check rejects an unrelated 200-OK service before uploading study data',async({page})=>{
  let uploaded=0;
  await page.route('http://127.0.0.1:8787/health',async route=>{
    await route.fulfill({status:200,contentType:'application/json',
      body:JSON.stringify({status:'ok',version:'1.4.0',packages:{DESeq2:true},
        capabilities:{reference_analysis:true}})});
  });
  await page.route('http://127.0.0.1:8787/run',async route=>{
    uploaded++;
    await route.fulfill({status:200,contentType:'application/json',
      body:JSON.stringify({status:'ok',methods:{}})});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:25000});
  await page.getByTestId('multiomics-analysis-intent').locator('select').selectOption('confirmatory');
  await page.getByTestId('multiomics-run-options').locator('summary').click();
  await expect(page.getByLabel(/Mode backend R|R backend mode/)).toHaveValue('required');
  await page.getByTestId('multiomics-run-options').getByRole('button',{name:'Tester'}).click();
  await expect(page.getByTestId('multiomics-run-options')).toContainText('ne présente pas le contrat attendu');
  await page.getByTestId('multiomics-run').click();
  await expect(page.locator('.error').filter({hasText:'contrat attendu'})).toBeVisible({timeout:25000});
  expect(uploaded).toBe(0);
});

test('recognized R health remains only method readiness, never scientific certification',async({page})=>{
  await page.route('http://127.0.0.1:8787/health',async route=>{
    await route.fulfill({status:200,contentType:'application/json',body:JSON.stringify({
      status:'ok',engine:'PMx Explain reference R backend',version:'1.4.0',
      runtime:{RVersion:'4.5.2',RPlatform:'x86_64'},
      capabilities:{reference_analysis:true,mofa2_orientation_contract:'mofa2-feature-rows-sample-columns-v2'},
      packages:{DESeq2:true,edgeR:true,limma:true,lmerTest:true,
        survival:true,MOFA2:true,mixOmics:true,fgsea:true}
    })});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-run-options').locator('summary').click();
  await page.getByLabel(/Mode backend R|R backend mode/).selectOption('auto');
  await page.getByTestId('multiomics-run-options').getByRole('button',{name:'Tester'}).click();
  await expect(page.getByTestId('multiomics-run-options')).toContainText('les méthodes réellement exécutées seront contrôlées');
  await expect(page.getByTestId('multiomics-run-options')).not.toContainText('scientifiquement validé');
});
