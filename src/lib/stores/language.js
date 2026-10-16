import { browser } from '$app/environment';
import { writable } from 'svelte/store';
const key = 'multiomics-language';
function initial() {
  if (!browser) return 'fr';
  const param = new URLSearchParams(location.search).get('lang');
  if (param === 'en' || param === 'fr') return param;
  try { return localStorage.getItem(key) === 'en' ? 'en' : 'fr'; } catch { return 'fr'; }
}
const store = writable(initial());
if (browser) store.subscribe(lang => {
  document.documentElement.lang = lang;
  try { localStorage.setItem(key, lang); } catch { /* private mode */ }
});
export const language = {
  subscribe: store.subscribe,
  set(lang) {
    const next = lang === 'en' ? 'en' : 'fr';
    store.set(next);
    if (browser) {
      const url = new URL(location.href);
      if (next === 'en') url.searchParams.set('lang', 'en');
      else url.searchParams.delete('lang');
      history.replaceState(history.state, '', url);
    }
  }
};
