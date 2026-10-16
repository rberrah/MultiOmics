// @ts-nocheck
import {test, expect} from '@playwright/test';

for(const width of [1280,390])for(const lang of ['fr','en']){
 test(`Omics-PKPD cross-domain learning is usable in ${lang} at ${width}px`,
   async({page})=>{
   await page.setViewportSize({width,height:950});
   const pageErrors=[];
   page.on('pageerror',error=>pageErrors.push(error.message));
   await page.goto(`/domaines/?lang=${lang}`);
   const content=page.getByTestId('cross-domain-courses');
   await expect(content).toBeVisible();
   await expect(page.getByTestId('cross-domain-open-tool'))
     .toHaveAttribute('href','https://rberrah.github.io/multiomics/tool/');
   for(const id of ['match','covariate','pd','mechanistic']){
     const lesson=page.getByTestId(`cross-domain-lesson-${id}`);
     await expect(lesson).toBeVisible();
     await expect(lesson.locator('.lesson-sources summary')).toContainText(
       lang==='en'?'Sources & further reading':'Sources et lectures complémentaires');
     await expect(lesson.locator('.lesson-sources a')).toHaveCount(2);
   }
   await expect(page.getByTestId('cross-domain-auc')).toContainText(
     lang==='en'?'200.0 ng·h/mL':'200,0 ng·h/mL');
   await expect(page.getByTestId('cross-domain-cl')).toContainText(
     lang==='en'?'25.0 L/h':'25,0 L/h');
   await page.getByTestId('cross-domain-signature').fill('2');
   await expect(page.getByTestId('cross-domain-ratio'))
     .toContainText(lang==='en'?'0.61':'0,61');
   await page.getByTestId('cross-domain-beta').fill('0');
   await expect(page.getByTestId('cross-domain-ratio'))
     .toContainText(lang==='en'?'1.00':'1,00');
   for(const [id,answer] of [['match',1],['covariate',1],['pd',1],['mechanistic',1]]){
     await page.getByTestId(`cross-domain-answer-${id}-${answer}`).click();
     await expect(page.getByTestId(`cross-domain-feedback-${id}`))
       .toContainText(lang==='en'?'Correct.':'Correct.');
   }
   await expect(page.getByTestId('cross-domain-progress'))
     .toContainText('4 / 4');
   expect(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1)).toBe(true);
   expect(pageErrors).toEqual([]);
 });
}

test('Cross-domain navigation is discoverable and accessible from main menu',async({page})=>{
 await page.goto('/?lang=fr');
 const menu=page.getByTestId('goal-connections');
 await expect(menu).toBeVisible();
 await menu.locator('summary').click();
 await expect(menu.getByTestId('nav-cross-domains'))
   .toHaveAttribute('href',/\/domaines\/?$/);
 await menu.getByTestId('nav-cross-domains').click();
 await expect(page.getByTestId('cross-domain-courses')).toBeVisible();
});

test('Misconceptions about causality and time are corrected, not claimed as clinical proof',async({page})=>{
 await page.goto('/domaines/?lang=fr');
 await page.getByTestId('cross-domain-answer-pd-0').click();
 await expect(page.getByTestId('cross-domain-feedback-pd')).toContainText('À revoir.');
 await expect(page.getByTestId('cross-domain-lesson-mechanistic'))
   .toContainText('Expression ARN et activité enzymatique ne sont pas équivalentes.');
 await expect(page.getByTestId('cross-domain-courses'))
   .toContainText('sans données patients ni recommandation de dose');
});
