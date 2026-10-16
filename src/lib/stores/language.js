import {browser} from '$app/environment';
import {writable} from 'svelte/store';
const key='multiomics-language';
const supported=new Set(['fr','en']);
function initial() {
  if(!browser)return 'fr';
  const fromUrl=new URLSearchParams(window.location.search).get('lang');
  if(supported.has(fromUrl))return fromUrl;
  try {const local=localStorage.getItem(key);if(supported.has(local))return local;}catch{}
  return 'fr';
}
const source=writable(initial());
if(browser)source.subscribe(value=>{
  const normalized=supported.has(value)?value:'fr';
  try{localStorage.setItem(key,normalized);}catch{}
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
