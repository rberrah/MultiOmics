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

test('Guided demo is a PMx-style worked example driven by the exact Tool output',async({page})=>{
  const snapshot=()=>page.evaluate(()=>{
    const r=window.__PMX_MULTIOMICS_ANALYSIS__;
    if(!r)throw new Error('No published MultiOmics result');
    return {
      protocol:r.protocol,
      metadata:{
        subjects:r.metadataSummary?.subjects,
        samples:r.metadataSummary?.samples,
        assays:r.metadataSummary?.assays,
        conditions:r.metadataSummary?.conditions,
        timepoints:r.metadataSummary?.timepoints
      },
      layers:Object.fromEntries(Object.entries(r.layers||{}).map(([name,layer])=>[name,{
        contrast:layer.contrast,
        mode:layer.mode,
        inferenceMethod:layer.inferenceMethod,
        inferencePolicy:layer.inferencePolicy,
        effectScale:layer.effectScale,
        groupSizes:layer.groupSizes,
        selectionRule:layer.selectionRule,
        steps:layer.steps,
        rows:(layer.rows||[]).map(row=>({
          feature:row.feature,effect:row.effect,foldRatio:row.foldRatio,
          pValue:row.pValue,qValue:row.qValue,nReference:row.nReference,nComparison:row.nComparison
        }))
      }])),
      crossOmics:{
        method:r.crossOmics?.method,
        testedPairs:r.crossOmics?.testedPairs,
        significantPairs:r.crossOmics?.significantPairs,
        pairs:(r.crossOmics?.pairs||[]).map(pair=>({
          layerA:pair.layerA,featureA:pair.featureA,layerB:pair.layerB,featureB:pair.featureB,
          rReference:pair.rReference,rComparison:pair.rComparison,deltaR:pair.deltaR,
          pValue:pair.pValue,qValue:pair.qValue,pattern:pair.pattern
        }))
      },
      visualizations:r.visualizations
    };
  });

  await page.goto('/multiomics/demo?lang=fr');
  await expect(page.getByRole('heading',{name:/MultiOmics.*questionnaire.*sorties exactes/i})).toBeVisible();
  await expect(page.getByTestId('guided-demo-exact-results')).toContainText('STAT1',{timeout:30000});
  await expect(page.getByTestId('guided-demo-questionnaire')).toContainText('Mesures répétées');
  await expect(page.getByTestId('guided-demo-methods')).toContainText('deterministic permutation');
  await expect(page.getByTestId('guided-demo-methods')).toContainText('difference-in-differences');
  await expect(page.getByTestId('guided-demo-live-visualizations')).toBeVisible();
  await expect(page.getByTestId('guided-demo-limit')).toContainText(/synthétiques/i);
  const guided=await snapshot();

  await page.goto('/multiomics/tool/?demo=multi&lang=fr');
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(page.getByTestId('multiomics-visualizations')).toBeVisible();
  const tool=await snapshot();

  expect(guided).toEqual(tool);
  expect(guided.metadata).toEqual({
    subjects:8,samples:16,assays:49,
    conditions:['control','treatment'],timepoints:['T0','T12']
  });
  for(const layer of Object.values(guided.layers)){
    expect(layer.mode).toBe('difference-in-differences');
    expect(layer.inferenceMethod).toBe('deterministic permutation');
    expect(layer.groupSizes).toEqual([4,4]);
  }
});
