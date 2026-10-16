import { test, expect } from '@playwright/test';

test('four real demo choices, responsive order, hydration and result computation', async ({ page }) => {
  const errors = [];
  page.on('pageerror', e => errors.push(e.message));
  await page.goto('/multiomics/tool/');
  const ids = ['multiomics-quick-demo-rna','multiomics-quick-demo-protein','multiomics-quick-demo-metabolite','multiomics-quick-demo'];
  const buttons = ids.map(id => page.getByTestId(id));
  for (const button of buttons) await expect(button).toBeVisible();
  for (let i = 1; i < buttons.length; i++) {
    const a = await buttons[i-1].boundingBox(); const b = await buttons[i].boundingBox();
    expect(a && b).toBeTruthy();
    expect(b.y >= a.y - 3 || (Math.abs(b.y-a.y)<3 && b.x>a.x)).toBeTruthy();
  }
  const first = await buttons[0].boundingBox(); const last = await buttons[3].boundingBox();
  expect(last.y).toBeGreaterThan(first.y);
  await expect.poll(() => page.locator('html').getAttribute('data-pmx-multiomics-hydrated')).toBe('true');
  await buttons[3].click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 30000 });
  expect(errors).toEqual([]);
});

test('independent learning and overview routes; language switch', async ({ page }) => {
  await page.goto('/multiomics/');
  await expect(page.getByRole('link', {name:/Ouvrir l’outil|Open the tool/})).toBeVisible();
  await page.goto('/multiomics/learn/');
  await expect(page.getByRole('heading', {name:/Choisir une méthode|Choose a method/})).toBeVisible();
  await page.getByRole('button', {name:'Switch language'}).click();
  await expect(page.locator('html')).toHaveAttribute('lang','en');
});
