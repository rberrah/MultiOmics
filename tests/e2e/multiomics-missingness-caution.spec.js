// @ts-nocheck
import {test,expect} from '@playwright/test';

test('Any missing measurement causes a visible plain-language caveat, even with no group-rate imbalance',async({page})=>{
  // This is an adversarial presentation test, not a proof that observed
  // missingness reveals the unobservable MCAR/MAR/MNAR mechanism.
  await page.route('**/multiomics/demo_proteomics.csv',async route=>{
    const original=await route.fetch();
    const lines=(await original.text()).trim().split(/\r?\n/);
    const firstProtein=lines[1].split(',');
    firstProtein[1]=''; // just one missing value: median missingness can be zero
    lines[1]=firstProtein.join(',');
    await route.fulfill({status:200,contentType:'text/csv',body:lines.join('\n')+'\n'});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-quick-demo-protein').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const caveat=page.getByTestId('multiomics-missingness-caution');
  await expect(caveat).toBeVisible();
  await expect(caveat).toContainText('même si leur proportion est identique entre groupes');
  await expect(caveat).toContainText('analyse de sensibilité');
  const assurance=page.getByTestId('multiomics-scientific-assurance');
  await assurance.getByText('Voir les limites spécifiques à cette analyse').click();
  await expect(assurance).toContainText('Certaines mesures sont manquantes');
});
