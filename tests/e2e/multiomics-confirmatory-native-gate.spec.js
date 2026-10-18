// @ts-nocheck
import {test,expect} from '@playwright/test';

// This uses the ordinary synthetic longitudinal demo and deliberately enters
// confirmatory mode. The browser does not yet project complete R longitudinal
// coefficients into its displayed figures; no /run upload is justified.
test('unsupported longitudinal confirmatory intent blocks before any R data upload',async({page})=>{
  let uploads=0;
  await page.route('http://127.0.0.1:8787/run',async route=>{
    uploads++;
    await route.fulfill({status:200,headers:{'Access-Control-Allow-Origin':'*'},
      contentType:'application/json',body:'{"status":"ok","methods":{}}'});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:25000});
  await page.getByTestId('multiomics-analysis-intent').locator('select').selectOption('confirmatory');
  await expect(page.getByTestId('multiomics-confirmatory-projection-unsupported')).toContainText('limité aux comparaisons de groupes');
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toHaveCount(0);
  await expect(page.locator('.error').filter({hasText:'mode confirmatoire'})).toBeVisible();
  expect(uploads).toBe(0);
});
