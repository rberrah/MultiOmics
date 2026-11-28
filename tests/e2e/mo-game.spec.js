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

test('real-data illustrations, thresholds and split selection', async ({page})=>{
  await page.goto('/multiomics/game/?lang=en');
  const airway=page.getByTestId('public-study-airway');
  await expect(airway).toBeVisible();
  await expect(page.getByTestId('public-discovery-count')).toContainText('4,100');
  await page.getByTestId('public-q-threshold').selectOption({label:'0.01'});
  await expect(page.getByTestId('public-discovery-count')).toContainText('2,940');
  await page.getByRole('button',{name:/Nutrimouse/}).click();
  await expect(page.getByTestId('public-study-nutrimouse')).toBeVisible();
  await page.getByTestId('public-mouse-subset').selectOption('heldout');
  await expect(page.getByTestId('public-mouse-count')).toHaveText('10');
  await expect(page.locator('[data-mouse]')).toHaveCount(10);
  await page.getByTestId('language-toggle').getByRole('button',{name:'FR',exact:true}).click();
  await expect(page.getByTestId('public-mouse-count')).toHaveText('10');
  await page.goto('/multiomics/demo/?lang=en');
  await expect(page.getByTestId('public-study-airway')).toBeVisible();
  await expect(page.getByTestId('public-study-nutrimouse')).toBeVisible();
});

test('interpretation mission uses real results and a separate interactive batch simulation',async({page})=>{
  await page.goto('/multiomics/game/?lang=en');
  await page.getByRole('button',{name:/^Interpretation/}).click();
  await expect(page.getByTestId('mission-brief')).toContainText('Read a threshold');
  await page.getByTestId('public-q-threshold').selectOption({label:'0.01'});
  await page.getByRole('button',{name:'Fewer discoveries, without refitting'}).click();
  await page.getByRole('button',{name:'A hint',exact:true}).click();
  await expect(page.locator('.hint-box')).toContainText('Adjusted p-values stay fixed');
  await page.getByRole('button',{name:'Continue',exact:true}).click();
  await page.getByRole('button',{name:'In this contrast, estimated expression is higher under dexamethasone'}).click();
  await page.getByRole('button',{name:'Continue',exact:true}).click();
  await expect(page.getByTestId('game-batch-lab')).toBeVisible();
  await expect(page.getByTestId('batch-difference')).toHaveText('2');
  await page.getByTestId('cross-batches').check();
  await expect(page.getByTestId('batch-difference')).toHaveText('0');
  await page.getByTestId('language-toggle').getByRole('button',{name:'FR',exact:true}).click();
  await expect(page.getByTestId('cross-batches')).toBeChecked();
  await expect(page.getByTestId('batch-identifiability')).toContainText('Comparaisons intra-lot possibles');
  await page.getByRole('button',{name:'Le plan ne permet pas de séparer les deux effets'}).click();
  await page.getByRole('button',{name:'Continuer',exact:true}).click();
  await expect(page.getByTestId('public-study-nutrimouse')).toBeVisible();
});
