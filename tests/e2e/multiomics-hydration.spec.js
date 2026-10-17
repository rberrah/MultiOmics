import {test,expect} from '@playwright/test';

test('static analyzer boots without hydration errors and remains interactive', async ({page}) => {
  /** @type {string[]} */
  const javascriptErrors = [];
  page.on('pageerror',error=>javascriptErrors.push(error.message));
  page.on('console',message=>{
    if(message.type()==='error' && /hydrate|HierarchyRequestError|appendChild/i.test(message.text()))
      javascriptErrors.push(message.text());
  });
  await page.goto('/multiomics/tool/');
  await expect(page.getByTestId('multiomics-quick-demo-rna')).toBeVisible();
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  expect(javascriptErrors).toEqual([]);
});