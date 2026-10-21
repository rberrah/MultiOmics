import {test,expect} from '@playwright/test';

test('Beginner omics lessons: method selection, warnings and return to analyzer (FR)',async({page})=>{
  /** @type {string[]} */
  const errors=[];
  page.on('pageerror',e=>errors.push(e.message));
  await page.goto('/multiomics/learn?lang=fr');
  await expect(page.getByRole('heading',{name:/Choisir une méthode et comprendre ses résultats/i})).toBeVisible();
  await expect(page.getByRole('heading',{name:/Comment lire les figures et statistiques/i})).toBeVisible();
  const deseq=page.locator('details.lesson').filter({has:page.locator('summary', {hasText:'DESeq2 / edgeR-QL'})});
  await deseq.locator('summary').click();
  await expect(deseq).toContainText(/comptages RNA-seq bruts entiers/i);
  await expect(page.locator('.caution')).toContainText(/JavaScript.*simplifiés/i);
  await expect(page.getByText(/ANALYSE_REPRODUCTIBLE.Rmd/)).toBeVisible();
  await expect(page.getByText(/Rscript render_notebook.R/)).toBeVisible();
  await page.getByRole('link',{name:/Ouvrir MultiOmics/}).click();
  await expect(page).toHaveURL(/\/multiomics\/tool/);
  expect(errors).toEqual([]);
});

test('Omics homepage links to a readable course',async({page})=>{
  await page.goto('/multiomics');
  const course=page.getByRole('link',{name:/Cours : méthodes et résultats|Learn methods and interpretation/});
  await expect(course).toBeVisible();
  await course.click();
  await expect(page).toHaveURL(/\/multiomics\/learn/);
});

test('Quick-start demo has three single-omic choices above a full-width multi-omic action',async({page})=>{
  await page.setViewportSize({width:1280,height:900});
  await page.goto('/multiomics/tool?lang=fr');
  const action=page.getByTestId('multiomics-quick-start');
  await expect(action).toBeVisible();
  const selectors=[
    'multiomics-quick-demo-rna',
    'multiomics-quick-demo-protein',
    'multiomics-quick-demo-metabolite'
  ];
  const singles=selectors.map(id=>page.getByTestId(id));
  const multi=page.getByTestId('multiomics-quick-demo');
  await expect(singles[0]).toHaveText('Transcriptomique');
  await expect(singles[1]).toHaveText('Protéomique');
  await expect(singles[2]).toHaveText('Métabolomique');
  await expect(multi).toContainText('Multi-omiques');
  const boxes=await Promise.all([...singles,multi].map(x=>x.boundingBox()));
  expect(boxes.every(x=>x!==null)).toBe(true);
  const [rna,protein,metabolite,all]=boxes;
  if(!rna || !protein || !metabolite || !all) throw new Error('Missing one or more quick-demo button positions');
  expect(Math.abs(rna.y-protein.y)).toBeLessThan(3);
  expect(Math.abs(protein.y-metabolite.y)).toBeLessThan(3);
  expect(all.y).toBeGreaterThan(rna.y+rna.height-3);
  expect(all.width).toBeGreaterThan(rna.width*2.5);
});

test('Mobile 390px: three single-omic shortcuts and a full-width multi-omics button',async({page})=>{
  await page.setViewportSize({width:390,height:844});
  await page.goto('/multiomics/tool?lang=fr');
  const ids=['multiomics-quick-demo-rna','multiomics-quick-demo-protein','multiomics-quick-demo-metabolite','multiomics-quick-demo'];
  const boxes=[];
  for(const id of ids){
    const el=page.getByTestId(id);
    await expect(el).toBeVisible();
    const area=await el.boundingBox();
    if(!area)throw new Error('Mobile shortcut cannot be measured: '+id);
    boxes.push(area);
  }
  const [a,b,c,all]=boxes;
  expect(Math.abs(a.y-b.y)).toBeLessThan(3);
  expect(Math.abs(b.y-c.y)).toBeLessThan(3);
  expect(all.y).toBeGreaterThan(a.y+a.height-3);
  expect(all.width).toBeGreaterThan(a.width*2.5);
  const overflow=await page.evaluate(()=>document.documentElement.scrollWidth-document.documentElement.clientWidth);
  expect(overflow).toBeLessThanOrEqual(12);
  await page.getByTestId('multiomics-quick-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(page.getByTestId('multiomics-visualizations')).toBeVisible();
});

test('Guided demo explains the existing synthetic longitudinal case and can launch it directly',async({page})=>{
  await page.goto('/multiomics/demo?lang=fr');
  await expect(page.getByRole('heading',{name:/Cas pratique.*traitement/i})).toBeVisible();
  await expect(page.getByTestId('guided-demo-design')).toContainText('8 sujets');
  await expect(page.getByTestId('guided-demo-design')).toContainText('T0');
  await expect(page.getByTestId('guided-demo-interpretation')).toContainText('tryptophane');
  await expect(page.getByTestId('guided-demo-limit')).toContainText(/synthétique/i);
  const launch=page.getByTestId('guided-demo-launch');
  await expect(launch).toHaveAttribute('href',/\/multiomics\/tool\/\?demo=multi$/);
  await launch.click();
  await expect(page).toHaveURL(/\/multiomics\/tool\/\?demo=multi$/);
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText(/Multi-omique/);
});
