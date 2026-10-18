// @ts-nocheck
import {test,expect} from '@playwright/test';

test('Outdated MOFA2 localhost backend is refused before uploading research matrices',async({page})=>{
  let uploads=0;
  await page.route('http://127.0.0.1:8787/health',async route=>{
    await route.fulfill({status:200,headers:{'Access-Control-Allow-Origin':'*'},contentType:'application/json',
      body:JSON.stringify({status:'ok',engine:'PMx Explain reference R backend',
        version:'1.3.0',runtime:{RVersion:'4.5.2'},
        packages:{DESeq2:true,edgeR:true,limma:true,lmerTest:true,
          survival:true,MOFA2:true,mixOmics:true,fgsea:true},
        capabilities:{reference_analysis:true}})});
  });
  await page.route('http://127.0.0.1:8787/run',async route=>{
    uploads++;
    await route.fulfill({status:200,headers:{'Access-Control-Allow-Origin':'*'},contentType:'application/json',
      body:JSON.stringify({status:'ok',methods:{mofa2:{status:'ok'}}})});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:25000});
  // The general synthetic demo defaults to a longitudinal model, not MOFA2.
  // Switch to the exploratory 3-omic plan before checking MOFA2 compatibility.
  await page.getByLabel(/Objectif principal|Main objective/).selectOption('explore');
  await page.getByTestId('multiomics-run-options').locator('summary').click();
  await page.getByLabel(/Mode backend R|R backend mode/).selectOption('auto');
  await page.getByTestId('multiomics-run').click();
  const warning=page.getByTestId('multiomics-mofa-version-warning');
  await expect(warning).toBeVisible({timeout:25000});
  await expect(warning).toContainText('Mettez à jour');
  await expect(page.getByTestId('multiomics-results')).toBeVisible();
  expect(uploads).toBe(0);
});
