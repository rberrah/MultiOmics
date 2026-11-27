import {test, expect} from '@playwright/test';

test('Mo Game keeps progress and answer feedback during language changes', async ({page}) => {
  await page.goto('/multiomics/game/?lang=en');
  await expect(page.getByRole('heading',{name:'Mo Game',exact:true})).toBeVisible();
  await expect(page.getByTestId('mo-game-study')).toContainText('17,199');
  await page.getByRole('button',{name:'Donor-paired comparison'}).click();
  await expect(page.getByRole('status')).toContainText('Correct.');
  await page.getByTestId('language-toggle').getByRole('button',{name:'FR',exact:true}).click();
  await expect(page.getByRole('status')).toContainText('Correct.');
  await expect(page.getByRole('button',{name:'Comparaison appariée par donneur'})).toHaveAttribute('aria-pressed','true');
  await page.getByRole('button',{name:'Continuer',exact:true}).click();
  await expect(page.getByRole('heading',{name:'Airway · Étape 2 / 4'})).toBeVisible();
  await page.getByRole('button',{name:/Nutrimouse/}).click();
  await expect(page.getByTestId('mo-game-study')).toContainText('40 souris');
  await page.getByRole('button',{name:/Airway/}).click();
  await expect(page.getByRole('heading',{name:'Airway · Étape 2 / 4'})).toBeVisible();
  await page.getByRole('link',{name:'Retour à MultiOmics'}).click();
  await expect(page.getByTestId('mo-game-cta')).toBeVisible();
});

test('Mo Game mobile layout and incomplete answer guard', async ({page}) => {
  await page.setViewportSize({width:320,height:720});
  await page.goto('/multiomics/game/?lang=en');
  await page.getByRole('button',{name:'Independent groups'}).click();
  await expect(page.getByRole('button',{name:'Continue',exact:true})).toBeDisabled();
  await expect(page.getByRole('status')).toContainText('Try again.');
  expect(await page.evaluate(()=>document.documentElement.scrollWidth-document.documentElement.clientWidth)).toBeLessThanOrEqual(2);
});
