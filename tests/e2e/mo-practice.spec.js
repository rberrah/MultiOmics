import {test, expect} from '@playwright/test';

test('independent omics practice, worked answers and in-place FR/EN', async ({page}) => {
  await page.goto('/multiomics/practice/?lang=en');
  await expect(page.getByRole('heading',{name:'Practice',exact:true})).toBeVisible();
  const answers=page.locator('main details');
  await expect(answers).toHaveCount(3);
  await answers.first().locator('summary').click();
  await expect(answers.first()).toContainText('Four, not eight.');
  await expect(answers.first()).toHaveAttribute('open','');
  await page.getByTestId('language-toggle').getByRole('button',{name:'FR',exact:true}).click();
  await expect(page.getByRole('heading',{name:'Exercices',exact:true})).toBeVisible();
  await expect(answers.first()).toHaveAttribute('open','');
  await expect(answers.first()).toContainText('Quatre, et non huit.');
  await expect(answers.nth(1)).not.toHaveAttribute('open','');
  await page.getByRole('link',{name:'Retour à Omics Explain'}).click();
  await expect(page.getByRole('link',{name:'Faire les exercices'})).toBeVisible();
});
