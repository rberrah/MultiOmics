// @ts-nocheck
import {test,expect} from '@playwright/test';

// The ordinary synthetic demo is longitudinal. Once confirmatory intent is
// selected, longitudinal group-by-time is a supported native-R projection
// family and must no longer be rejected by the UI before R preflight.
test('longitudinal confirmatory intent is offered when native R projection is supported',async({page})=>{
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:25000});
  await page.getByTestId('multiomics-advanced-analysis-options').locator('summary').click();
  await page.getByTestId('multiomics-analysis-intent').locator('select').selectOption('confirmatory');
  await expect(page.getByLabel(/Mode backend R|R backend mode/)).toHaveValue('required');
  await expect(page.getByTestId('multiomics-confirmatory-projection-unsupported')).toHaveCount(0);
});
