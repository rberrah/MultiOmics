import { test, expect } from '@playwright/test';

test('Downloaded portable package serves real isolated routes, HTML and scripts', async ({ page }) => {
  for (const route of ['/','/multiomics/','/multiomics/learn/','/multiomics/tool/']) {
    const response=await page.goto(route,{waitUntil:'domcontentloaded'});
    expect(response?.status(),'HTTP '+route).toBe(200);
  }
  /** @type {string[]} */
  const pageErrors=[];
  page.on('pageerror',err=>pageErrors.push(err.message));
  await page.goto('/multiomics/tool/');
  await page.waitForFunction(()=>document.documentElement.dataset.pmxMultiomicsHydrated==='true');
  await expect(page.getByTestId('multiomics-quick-demo-rna')).toBeVisible();
  await expect(page.getByTestId('multiomics-quick-demo-protein')).toBeVisible();
  await expect(page.getByTestId('multiomics-quick-demo-metabolite')).toBeVisible();
  await expect(page.getByTestId('multiomics-quick-demo')).toBeVisible();
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  expect(pageErrors).toEqual([]);
});

test('Portable standalone tool works at mobile width without PMx assets', async ({page})=>{
  await page.setViewportSize({width:390,height:844});
  await page.goto('/multiomics/tool/');
  await page.waitForFunction(()=>document.documentElement.dataset.pmxMultiomicsHydrated==='true');
  const overflow=await page.evaluate(()=>document.documentElement.scrollWidth-document.documentElement.clientWidth);
  expect(overflow).toBeLessThanOrEqual(12);
  await page.getByTestId('multiomics-quick-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const assets=await page.evaluate(()=>[...document.querySelectorAll('script[src]')].map(s=>s.getAttribute('src')));
  expect(assets.some(x=>x?.includes('/pharmacometrie/'))).toBe(false);
});
