// @ts-nocheck
import {test,expect} from '@playwright/test';

test('A new researcher can finish a one-omic demo on a phone without R or jargon-heavy setup',async({page})=>{
  await page.setViewportSize({width:390,height:844});
  let backendCalls=0;
  page.on('request',request=>{
    if(/127\.0\.0\.1:8787\/run|reactome\.org\/AnalysisService/.test(request.url()))
      backendCalls++;
  });
  await page.goto('/multiomics/tool?lang=fr');
  const novice=page.getByTestId('multiomics-quick-start');
  await expect(novice).toBeVisible();
  await expect(novice).toContainText('Découvrir en 1 clic');
  await expect(novice).toContainText('Aucun fichier à préparer');
  await page.getByTestId('multiomics-quick-demo-rna').focus();
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText('Omique unique');
  expect(backendCalls).toBe(0);
  const initialOverflow=await page.evaluate(()=>
    document.documentElement.scrollWidth-document.documentElement.clientWidth);
  expect(initialOverflow).toBeLessThanOrEqual(12);
});


test('Three-omics beginner demo never requires R or Reactome service',async({page})=>{
  let external=0;
  page.on('request',request=>{
    if(/127\.0\.0\.1:8787\/run|reactome\.org\/AnalysisService/.test(request.url()))external++;
  });
  await page.goto('/multiomics/tool');
  await page.getByTestId('multiomics-quick-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText(/multi|omiques/i);
  expect(external).toBe(0);
});

test('Beginner one-click demo instructions meet WCAG AA 4.5:1 text contrast',async({page})=>{
  await page.goto('/multiomics/tool?lang=fr');
  const ratios=await page.evaluate(()=>{
    const surface=document.querySelector('.quick-start');
    if(!surface)throw new Error('Beginner intro panel missing');
    const background=getComputedStyle(surface).backgroundColor;
    const parse=(value)=>{
      const parts=value.match(/[\d.]+/g)?.slice(0,3).map(Number);
      if(!parts||parts.length<3)throw new Error('Cannot evaluate intro contrast: '+value);
      return parts;
    };
    const luminance=(value)=>{
      const components=parse(value).map(n=>{
        const srgb=n/255;
        return srgb<=0.04045?srgb/12.92:((srgb+0.055)/1.055)**2.4;
      });
      return components[0]*0.2126+components[1]*0.7152+components[2]*0.0722;
    };
    const bg=luminance(background);
    return ['strong','span','small'].map(selector=>{
      const node=surface.querySelector('.quick-start-copy '+selector);
      if(!node)throw new Error('Missing beginner instructions: '+selector);
      const fg=luminance(getComputedStyle(node).color);
      return {selector,ratio:(Math.max(bg,fg)+0.05)/(Math.min(bg,fg)+0.05)};
    });
  });
  for(const item of ratios)
    expect(item.ratio,'Low-contrast beginner text: '+item.selector).toBeGreaterThanOrEqual(4.5);
});

test('Beginner-facing DIABLO/MOFA explanations define terms and avoid biomarker certification',async({page})=>{
  await page.route('http://127.0.0.1:8787/health',async route=>{
    await route.fulfill({status:200,contentType:'application/json',
      body:JSON.stringify({status:'ok',engine:'PMx reference',
        packages:{mixOmics:true,MOFA2:true},
        capabilities:{mofa2_orientation_contract:'mofa2-feature-rows-sample-columns-v2'}})});
  });
  await page.route('http://127.0.0.1:8787/run',async route=>{
    await route.fulfill({status:200,contentType:'application/json',
      body:JSON.stringify({status:'ok',engine:{name:'PMx R reference',version:'1.3.0'},
        methods:{
          diablo_heldout:{status:'ok',method:'Internal heldout DIABLO',
            summary:{meanBER:0.25,folds:[{},{},{}],
              signatureStability:{
                transcriptomics:{meanPairwiseJaccard:0.50,repeatedlySelected:['G1','G2']},
                proteomics:{meanPairwiseJaccard:0.34,repeatedlySelected:['P1']}
              }}},
          mofa2:{status:'ok',method:'MOFA2',
            summary:{inputAudit:{nSharedSubjects:48},matrixOrientation:'feature_rows_subject_columns_for_MOFA2',subjectIdentityPreserved:true,implementationContract:'mofa2-feature-rows-sample-columns-v2'}}
        }})});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:20000});
  await page.getByTestId('multiomics-run-options').locator('summary').first().click();
  await page.getByLabel(/Mode backend R|R backend mode/).selectOption('auto');
  await page.getByTestId('multiomics-run').click();
  await page.getByTestId('multiomics-quality-details').locator('summary').first().click();
  await expect(page.getByTestId('reference-backend-results')).toBeVisible({timeout:20000});
  const stability=page.getByTestId('diablo-signature-stability');
  await expect(stability).toBeVisible();
  await stability.locator('summary').click();
  await expect(stability).toContainText('Jaccard : de 0');
  await expect(stability).toContainText('ne valide pas un biomarqueur');
  const mofa=page.getByTestId('mofa2-input-audit');
  await expect(mofa).toContainText('facteur latent');
  await expect(mofa).toContainText('effet de lot');
});


test('A failed demo file download shows a clear message and a retry actually works',async({page})=>{
  const jsErrors=[];
  page.on('pageerror',error=>jsErrors.push(error.message));
  await page.route('**/multiomics/demo_transcriptomics.csv',route=>
    route.fulfill({status:503,contentType:'text/plain',body:'Temporary unavailable'}));
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  const warning=page.getByTestId('multiomics-demo-error');
  await expect(warning).toBeVisible({timeout:15000});
  await expect(warning).toContainText('Réessayez');
  await expect(page.getByTestId('multiomics-results')).toHaveCount(0);
  await expect(page.getByTestId('multiomics-quick-demo-rna')).toBeEnabled();
  expect(jsErrors).toEqual([]);
  await page.unroute('**/multiomics/demo_transcriptomics.csv');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(warning).toHaveCount(0);
  await expect(page.getByTestId('multiomics-analysis-mode')).toContainText('Omique unique');
  expect(jsErrors).toEqual([]);
});


test('Real MOFA2 negative pairing result is understandable and disclosed in French and English',async({page})=>{
  for(const [locale,claim] of [
    ['fr', 'ne prouve donc pas que les deux omiques sont correctement appariées'],
    ['en', 'does not prove that two omics were correctly paired']
  ]){
    await page.goto('/multiomics/tool?lang='+locale);
    const panel=page.getByTestId('multiomics-validation-evidence');
    await expect(panel).toBeVisible();
    await panel.locator('details').first().locator('summary').first().click();
    const negative=page.getByTestId('mofa-real-negative-evidence');
    await expect(negative).toBeVisible();
    await expect(negative).toContainText(claim);
    await expect(negative).toContainText('0.5313'.replace('.',locale==='fr'?',':'.'));
    await expect(negative).toContainText(locale==='fr'?'normalisation':'normalization');
    await expect(negative).toContainText(locale==='fr'?'0,9359':'0.9359');
    await expect(negative).toContainText(locale==='fr'?'Le résultat dépend donc fortement du prétraitement':'Changing the preprocessing substantially changed');
  }
});


test('R1 confirmatory user-selected RNA direction is sent verbatim to the reference R backend',async({page})=>{
  let submitted=null;
  await page.route('http://127.0.0.1:8787/health',async route=>
    route.fulfill({status:200,contentType:'application/json',
      body:JSON.stringify({status:'ok',version:'1.4.0',
        capabilities:{reference_analysis:true,
          mofa2_orientation_contract:'mofa2-feature-rows-sample-columns-v2'}})}));
  await page.route('http://127.0.0.1:8787/run',async route=>{
    submitted=route.request().postDataJSON();
    await route.fulfill({status:200,contentType:'application/json',
      body:JSON.stringify({status:'ok',methods:{
        transcriptomics_differential_deseq2:{
          status:'blocked',message:'Synthetic demo: no R1 scientific confirmation'
        }
      }})});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await page.getByLabel(/Objectif principal|Main objective/).selectOption('groups');
  // The RNA demo starts with a repeated-measures design. Configure the
  // independent design explicitly through the genuine, visible advanced UI.
  await page.getByTestId('multiomics-advanced-settings').locator('summary').click();
  await page.getByLabel(/Structure du design|Design structure/).selectOption('independent');
  await page.getByLabel(/Design longitudinal|Longitudinal design/).selectOption('no');
  await page.getByTestId('multiomics-analysis-intent').locator('select').selectOption('confirmatory');
  const chooser=page.getByTestId('multiomics-rnaseq-contrast');
  await expect(chooser).toBeVisible();
  const fields=chooser.locator('select');
  await expect(fields).toHaveCount(2);
  const levels=await fields.first().locator('option').evaluateAll(nodes=>
    nodes.map(n=>n.value).filter(Boolean));
  expect(levels.length).toBe(2);
  await expect(fields.nth(0)).toHaveValue('');
  await expect(fields.nth(1)).toHaveValue('');
  await fields.nth(0).selectOption(levels[0]);
  await fields.nth(1).selectOption(levels[1]);
  await page.getByTestId('multiomics-run').click();
  await expect.poll(()=>submitted?.protocol?.analysisIntent,{timeout:30000}).toBe('confirmatory');
  expect(submitted.protocol.rnaPrimaryContrast).toEqual(['condition',levels[1],levels[0]]);
  expect(submitted.protocol.rnaPrimaryContrast[1]).not.toBe(
    submitted.protocol.rnaPrimaryContrast[2]);
  // A blocked synthetic response is intentionally not a pass/certificate.
  expect(submitted.dataTypes.transcriptomics).toBe('raw_counts');
});


test('Beginner R backend setup is discoverable, precise and mobile readable without installing R',async({page})=>{
  await page.setViewportSize({width:390,height:844});
  let submissions=0;
  page.on('request',request=>{
    if(request.url().endsWith('127.0.0.1:8787/run'))submissions++;
  });
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:18000}).toBe('true');
  const help=page.getByTestId('multiomics-local-r-setup');
  await expect(help.locator('summary')).toContainText('Comment activer le moteur R');
  await help.locator('summary').click();
  await expect(help).toContainText('GitHub Pages ne peut pas exécuter R');
  await expect(help).toContainText('Installer R et Git');
  await expect(help).toContainText('Auto · utiliser si disponible');
  await expect(help.locator('pre')).toContainText('Rscript multiomics-engine/run_backend.R');
  await expect(help.locator('pre')).toContainText('Rscript multiomics-engine/install_backend_dependencies.R');
  await expect(help).toContainText('127.0.0.1');
  await expect(help.locator('a')).toHaveAttribute('href',
     'https://github.com/rberrah/rberrah.github.io/blob/main/multiomics-engine/README.md');
  const overflow=await page.evaluate(()=>
    document.documentElement.scrollWidth-document.documentElement.clientWidth);
  expect(overflow).toBeLessThanOrEqual(12);
  expect(submissions).toBe(0);
  await page.goto('/multiomics/tool?lang=en');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:18000}).toBe('true');
  const english=page.getByTestId('multiomics-local-r-setup');
  await expect(english.locator('summary')).toContainText('How to enable R');
  await english.locator('summary').click();
  await expect(english).toContainText('No installation is needed to explore data');
  await expect(english).toContainText('Never send patient data');
  await expect(english).toContainText('Auto · use when available');
  expect(submissions).toBe(0);
});


test('R beginner setup: Browser only refuses health checks; Auto can verify local R without uploading a study',async({page})=>{
  let health=0,uploaded=0;
  await page.route('http://127.0.0.1:8787/health',async route=>{
    health++;
    await route.fulfill({status:200,contentType:'application/json',
      body:JSON.stringify({status:'ok',version:'1.8.0',capabilities:{reference_analysis:true}})});
  });
  await page.route('http://127.0.0.1:8787/run',async route=>{
    uploaded++;
    await route.fulfill({status:500,body:'Unexpected research upload'});
  });
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  const advanced=page.getByTestId('multiomics-run-options');
  await advanced.locator(':scope > summary').click();
  const mode=page.getByLabel('Mode backend R');
  await expect(mode).toHaveValue('browser');
  await advanced.getByRole('button',{name:'Tester'}).click();
  await expect(advanced).toContainText('Backend R désactivé');
  expect(health).toBe(0);
  expect(uploaded).toBe(0);
  await mode.selectOption('auto');
  await advanced.getByRole('button',{name:'Tester'}).click();
  await expect(advanced).toContainText('Backend R disponible');
  expect(health).toBe(1);
  expect(uploaded).toBe(0);
});
