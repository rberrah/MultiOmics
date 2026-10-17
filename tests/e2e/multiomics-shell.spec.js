import {test,expect} from '@playwright/test';

test('PMx identical editorial dropdown shell with mO monogram', async ({page}) => {
  await page.goto('/multiomics/?lang=en');
  const header=page.getByTestId('site-header');
  await expect(header.locator('.mark')).toHaveText('mO');
  await expect(header.getByTestId('logo-link')).toHaveAttribute('href',/\/multiomics\/$/);
  for(const id of ['learn','analyze','disciplines','resources']){
    const details=header.getByTestId('goal-'+id);
    await expect(details).toBeVisible();
    await expect(details.locator('summary')).toHaveCount(1);
    await expect(details).not.toHaveAttribute('open','');
  }
  await header.getByTestId('goal-learn').locator('summary').click();
  await expect(header.getByTestId('nav-learn')).toBeVisible();
  await expect(header.getByTestId('nav-design')).toBeVisible();
  await header.getByTestId('goal-resources').locator('summary').click();
  await expect(header.getByTestId('goal-learn')).not.toHaveAttribute('open','');
  await expect(header.getByTestId('nav-glossary')).toHaveAttribute('href',/\/multiomics\/glossary\/$/);
  await expect(header.getByTestId('nav-references')).toHaveAttribute('href',/\/multiomics\/references\/$/);
  await expect(header.getByTestId('nav-main-site')).toHaveAttribute('href','https://rberrah.github.io/');
  await page.keyboard.press('Escape');
  await expect(header.getByTestId('goal-resources')).not.toHaveAttribute('open','');
  await expect(header.getByTestId('goal-resources').locator('summary')).toBeFocused();
  const css=await header.evaluate(e=>({position:getComputedStyle(e).position,height:e.getBoundingClientRect().height}));
  expect(css.position).toBe('sticky');
  expect(css.height).toBeGreaterThan(50);
});

test('Other disciplines dropdown directs to PMx Explain and Omics2PK', async ({page}) => {
  await page.goto('/multiomics/?lang=en');
  const group=page.getByTestId('goal-disciplines');
  await group.locator('summary').click();
  await expect(group.getByTestId('nav-disciplines')).toHaveAttribute('href',/\/multiomics\/disciplines\/$/);
  await expect(group.getByTestId('nav-pmx')).toHaveAttribute('href','https://rberrah.github.io/pharmacometrie/');
  await expect(group.getByTestId('nav-omics2pk')).toHaveAttribute('href','https://rberrah.github.io/Omics2PK/');
  await page.goto('/multiomics/glossary/?lang=en');
  await expect(page.getByRole('heading',{name:'Omics glossary'})).toBeVisible();
  await page.goto('/multiomics/references/?lang=en');
  await expect(page.getByRole('heading',{name:'Scientific references'})).toBeVisible();
});

test('mobile native dropdown menu, a11y escape and click-away', async ({page}) => {
  await page.setViewportSize({width:320,height:720});
  await page.goto('/multiomics/?lang=fr');
  const burger=page.getByTestId('nav-toggle');
  await expect(burger).toBeVisible();
  await burger.click();
  await expect(burger).toHaveAttribute('aria-expanded','true');
  await page.getByTestId('goal-resources').locator('summary').click();
  await expect(page.getByTestId('nav-main-site')).toBeVisible();
  await expect(page.getByTestId('language-toggle')).toBeVisible();
  await expect(page.getByTestId('theme-toggle')).toBeVisible();
  await page.getByTestId('nav-glossary').click();
  await expect(page).toHaveURL(/\/multiomics\/glossary\//);
  await expect(burger).toHaveAttribute('aria-expanded','false');
  await burger.click();
  await page.keyboard.press('Escape');
  await expect(burger).toHaveAttribute('aria-expanded','false');
  const overflow=await page.evaluate(()=>document.documentElement.scrollWidth-document.documentElement.clientWidth);
  expect(overflow).toBeLessThanOrEqual(12);
});

test('FR/EN and shared light/dark theme conventions survive navigation', async ({page})=>{
  await page.addInitScript(()=>{
    localStorage.setItem('pk-theme','dark');
    localStorage.setItem('pharmacometrie-language','en');
  });
  await page.goto('/multiomics/');
  await expect(page.locator('html')).toHaveAttribute('data-theme','dark');
  await expect(page.getByTestId('goal-resources').locator('summary')).toContainText('Resources');
  await page.getByTestId('language-toggle').getByRole('button',{name:'FR'}).click();
  await expect(page.getByTestId('goal-resources').locator('summary')).toContainText('Ressources');
  expect(await page.evaluate(()=>localStorage.getItem('pharmacometrie-language'))).toBe('fr');
  await page.getByTestId('theme-toggle').click();
  expect(await page.evaluate(()=>localStorage.getItem('pk-theme'))).toBe('system');
});
