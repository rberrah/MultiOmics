import { test, expect } from '@playwright/test';

test('standalone shell keeps navigation, margins and ecosystem links coherent', async ({ page }) => {
  await page.goto('/multiomics/');
  await expect(page.getByRole('link', { name: 'MultiOmics', exact: true })).toBeVisible();
  await expect(page.getByRole('link', { name: /Site principal|Main site/ }).first()).toHaveAttribute('href', 'https://rberrah.github.io/');
  await expect(page.getByRole('link', { name: /Autres disciplines|Other disciplines/ })).toHaveAttribute('href', /\/multiomics\/disciplines\/$/);
  await expect(page.getByRole('link', { name: /Ressources|Resources/ })).toHaveAttribute('href', /\/multiomics\/resources\/$/);

  const shell = page.locator('main#content');
  const box = await shell.boundingBox();
  const viewport = page.viewportSize();
  if (!box || !viewport) throw new Error('Missing shell geometry');
  const edges = await shell.evaluate(node => {
    const s = getComputedStyle(node);
    return { left: parseFloat(s.paddingLeft), right: parseFloat(s.paddingRight) };
  });
  expect(box.x + edges.left).toBeGreaterThanOrEqual(15);
  expect(viewport.width - (box.x + box.width) + edges.right).toBeGreaterThanOrEqual(15);
});

test('language and theme controls share PMx conventions and persist', async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.setItem('pk-theme', 'dark');
    localStorage.setItem('pharmacometrie-language', 'en');
  });
  await page.goto('/multiomics/');
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'dark');
  await expect(page.getByRole('link', { name: 'Other disciplines' })).toBeVisible();
  await expect(page.getByRole('link', { name: 'Resources' })).toBeVisible();

  await page.getByTestId('language-toggle').getByRole('button', { name: 'FR' }).click();
  await expect(page.getByRole('link', { name: 'Autres disciplines' })).toBeVisible();
  expect(await page.evaluate(() => localStorage.getItem('pharmacometrie-language'))).toBe('fr');

  await page.getByTestId('theme-toggle').click();
  expect(await page.evaluate(() => localStorage.getItem('pk-theme'))).toBe('system');
});

test('Other disciplines presents PMx and Resources remains a first-class page', async ({ page }) => {
  await page.goto('/multiomics/disciplines/');
  await expect(page.getByRole('heading', { name: /MultiOmics ne s’arrête pas|MultiOmics does not stop/ })).toBeVisible();
  await expect(page.getByRole('link', { name: /Ouvrir PMx Explain|Open PMx Explain/ })).toHaveAttribute('href', 'https://rberrah.github.io/pharmacometrie/');
  await expect(page.getByText(/Omics2PK/).first()).toBeVisible();

  await page.goto('/multiomics/resources/');
  await expect(page.getByRole('heading', { name: /Retrouver les méthodes|Find methods/ })).toBeVisible();
  await expect(page.getByRole('main').getByRole('link', { name: /GitHub/ })).toHaveAttribute('href', 'https://github.com/rberrah/MultiOmics');
  await expect(page.getByRole('link', { name: /Retour au site|Back to main site/ })).toHaveAttribute('href', 'https://rberrah.github.io/');
});

test('shell remains usable on mobile without losing navigation', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('/multiomics/tool/');
  const menu = page.getByTestId('nav-toggle');
  await expect(menu).toBeVisible();
  await expect(menu).toHaveAttribute('aria-expanded', 'false');
  await menu.click();
  await expect(menu).toHaveAttribute('aria-expanded', 'true');
  await expect(page.getByRole('link', { name: /Ressources|Resources/ })).toBeVisible();
  await expect(page.getByTestId('language-toggle')).toBeVisible();
  await expect(page.getByTestId('theme-toggle')).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(menu).toHaveAttribute('aria-expanded', 'false');
  await expect(menu).toBeFocused();
  await menu.click();
  await page.getByTestId('nav-resources').click();
  await expect(page).toHaveURL(/\/multiomics\/resources\/$/);
  await expect(menu).toHaveAttribute('aria-expanded', 'false');
  const overflow = await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
  expect(overflow).toBeLessThanOrEqual(12);
});


test('PMx compact editorial header preserves ecosystem navigation and theme', async ({ page }) => {
  await page.goto('/multiomics/');
  await expect(page).toHaveTitle(/MultiOmics/);
  await expect(page).not.toHaveTitle(/— PMx Explain/);
  const header = page.getByTestId('site-header');
  const mark = header.locator('.mark');
  const brand = page.getByTestId('logo-link');
  await expect(brand).toHaveAttribute('href', /\/multiomics\/$/);
  await expect(mark).toHaveText('Om');
  await expect(header.getByTestId('nav-pmx')).toHaveAttribute('href', 'https://rberrah.github.io/pharmacometrie/');
  await expect(header.getByTestId('nav-main-site')).toHaveAttribute('href', 'https://rberrah.github.io/');
  for (const key of ['analyze', 'learn', 'disciplines', 'resources']) {
    await expect(header.getByTestId('nav-' + key)).toBeVisible();
  }
  const metrics = await header.evaluate((node) => {
    const style = getComputedStyle(node);
    const badge = node.querySelector('.mark');
    if (!(badge instanceof Element)) throw new Error('Missing PMx-style monogram');
    const markStyle = getComputedStyle(badge);
    return { sticky: style.position, font: markStyle.fontFamily,
      background: markStyle.backgroundColor, minHeight: parseFloat(style.minHeight) };
  });
  expect(metrics.sticky).toBe('sticky');
  expect(metrics.minHeight).toBeGreaterThanOrEqual(58);
  expect(metrics.font.toLowerCase()).toContain('source serif');
  expect(metrics.background).not.toBe('rgba(0, 0, 0, 0)');
  await page.addInitScript(() => localStorage.setItem('pk-theme', 'dark'));
  await page.reload();
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'dark');
  await expect(header.getByTestId('nav-resources')).toBeVisible();
});

test('mobile navigation links and click-away work at 320px without overflow', async ({ page }) => {
  await page.setViewportSize({ width: 320, height: 720 });
  await page.goto('/multiomics/');
  const toggle = page.getByTestId('nav-toggle');
  await toggle.click();
  await expect(page.getByTestId('nav-disciplines')).toBeVisible();
  await expect(page.getByTestId('nav-pmx')).toBeVisible();
  await expect(page.getByTestId('nav-main-site')).toBeVisible();
  await page.getByRole('main').click({ position: { x: 6, y: 6 } });
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  const geometry = await page.evaluate(() => ({
    excess: document.documentElement.scrollWidth - document.documentElement.clientWidth,
    width: document.documentElement.clientWidth
  }));
  expect(geometry.excess).toBeLessThanOrEqual(12);
  expect(geometry.width).toBe(320);
});
