import {browser} from '$app/environment';
import {writable} from 'svelte/store';
const sharedKey='pharmacometrie-language';
const legacyKey='multiomics-language';
const supported=new Set(['fr','en']);
function initial() {
  if(!browser)return 'fr';
  const fromUrl=new URLSearchParams(window.location.search).get('lang');
  if(fromUrl!==null && supported.has(fromUrl))return fromUrl;
  try {
    const shared=localStorage.getItem(sharedKey);
    if(shared!==null && supported.has(shared))return shared;
    const legacy=localStorage.getItem(legacyKey);
    if(legacy!==null && supported.has(legacy))return legacy;
  } catch {}
  return 'fr';
}
const source=writable(initial());
if(browser)source.subscribe(value=>{
  const normalized=supported.has(value)?value:'fr';
  try {
    localStorage.setItem(sharedKey,normalized);
    localStorage.setItem(legacyKey,normalized);
  } catch {}
  document.documentElement.lang=normalized;
});
export const language={
  subscribe:source.subscribe,
  set(value){
    const normalized=supported.has(value)?value:'fr';
    source.set(normalized);
    if(browser){
      const u=new URL(location.href);
      if(normalized==='fr')u.searchParams.delete('lang');
      else u.searchParams.set('lang',normalized);
      if(u.href!==location.href)history.replaceState(history.state,'',u);
    }
  }
};
