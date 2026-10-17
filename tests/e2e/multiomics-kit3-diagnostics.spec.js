import {test,expect} from '@playwright/test';

test('Kit3 browser route and hydration stability diagnostic', async ({page})=>{
  const events=[];
  page.on('console', msg=>{
    if(events.length<120)events.push({type:'console-'+msg.type(),message:msg.text().slice(0,600)});
  });
  page.on('pageerror', error=>{
    if(events.length<120)events.push({type:'pageerror',message:error.message.slice(0,1000),stack:error.stack?.slice(0,1500)});
  });
  page.on('framenavigated', frame=>{
    if(frame===page.mainFrame()&&events.length<120)events.push({type:'navigation',url:frame.url()});
  });
  page.on('requestfailed',request=>{
    if(events.length<120)events.push({type:'requestfailure',url:request.url(),error:request.failure()?.errorText});
  });
  await page.addInitScript(()=>{
    window.__moRuns={documentLoads:(Number(sessionStorage.getItem('mo_loads')||'0')+1)};
    sessionStorage.setItem('mo_loads',String(window.__moRuns.documentLoads));
    let seen=0;
    const observer=new MutationObserver(()=>{
      const el=document.querySelector('[data-testid="multiomics-quick-demo-rna"]');
      if(el && !window.__moRuns.first)window.__moRuns.first=performance.now();
      if(el && window.__moRuns.last!==el){window.__moRuns.replacements=(window.__moRuns.replacements||0)+1;window.__moRuns.last=el;}
      if(++seen<200)window.__moRuns.observerNotifications=seen;
    });
    observer.observe(document,{subtree:true,childList:true});
  });
  const response=await page.goto('/multiomics/tool',{waitUntil:'domcontentloaded'});
  await page.waitForTimeout(4000);
  const state=await page.evaluate(()=>({
    url:location.href,readyState:document.readyState,title:document.title,
    counts:{rna:document.querySelectorAll('[data-testid="multiomics-quick-demo-rna"]').length,
      quick:document.querySelectorAll('[data-testid="multiomics-quick-start"]').length},
    renderText:document.body?.textContent?.slice(0,160),
    loads:window.__moRuns?.documentLoads,
    domReplacements:window.__moRuns?.replacements,
    mutationNotifications:window.__moRuns?.observerNotifications,
    activeTag:document.activeElement?.tagName
  }));
  console.log('KIT3_RUNTIME_DIAGNOSTIC',JSON.stringify({httpStatus:response?.status(),state,events}));
  expect(response?.status()).toBe(200);
  expect(state.counts.rna).toBe(1);
  expect(state.domReplacements).toBeLessThan(6);
  expect(events.filter(x=>x.type==='pageerror')).toHaveLength(0);
});