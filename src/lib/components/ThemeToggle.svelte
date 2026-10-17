<script>
  import { onMount } from 'svelte';
  import { language } from '#lib/stores/language.js';

  const KEY = 'pk-theme';
  /** @type {'system'|'light'|'dark'} */
  let mode = $state('system');

  /** @param {'system'|'light'|'dark'} value */
  function apply(value) {
    const root = document.documentElement;
    if (value === 'system') root.removeAttribute('data-theme');
    else root.setAttribute('data-theme', value);
  }

  onMount(() => {
    try {
      const saved = localStorage.getItem(KEY);
      if (saved === 'light' || saved === 'dark' || saved === 'system') mode = saved;
    } catch {}
    apply(mode);
  });

  function cycle() {
    mode = mode === 'system' ? 'light' : mode === 'light' ? 'dark' : 'system';
    apply(mode);
    try { localStorage.setItem(KEY, mode); } catch {}
  }

  const label = () => {
    if ($language === 'en') {
      if (mode === 'light') return 'Theme: light';
      if (mode === 'dark') return 'Theme: dark';
      return 'Theme: system';
    }
    if (mode === 'light') return 'Thème : clair';
    if (mode === 'dark') return 'Thème : sombre';
    return 'Thème : système';
  };
</script>

<button class="theme-toggle" type="button" onclick={cycle} title={label()} aria-label={label()} data-testid="theme-toggle">
  {#if mode === 'system'}
    <svg viewBox="0 0 24 24" aria-hidden="true"><rect x="3" y="4" width="18" height="13" rx="2"/><path d="M8 20h8"/></svg>
  {:else if mode === 'light'}
    <svg viewBox="0 0 24 24" aria-hidden="true"><circle cx="12" cy="12" r="4.2"/><path d="M12 2.5v2M12 19.5v2M2.5 12h2M19.5 12h2M5.2 5.2l1.4 1.4M17.4 17.4l1.4 1.4M18.8 5.2l-1.4 1.4M6.6 17.4l-1.4 1.4"/></svg>
  {:else}
    <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M20 14.5A8.5 8.5 0 1 1 9.5 4a6.8 6.8 0 0 0 10.5 10.5Z"/></svg>
  {/if}
</button>

<style>
  .theme-toggle{display:inline-flex;align-items:center;justify-content:center;width:34px;height:34px;padding:0;background:transparent;border:1px solid var(--border-subtle);border-radius:var(--radius);color:var(--text-secondary);cursor:pointer;transition:color .18s ease,border-color .18s ease,background-color .18s ease}
  .theme-toggle:hover{color:var(--accent-pk);border-color:var(--accent-pk);background:var(--bg-secondary)}
  svg{width:17px;height:17px;fill:none;stroke:currentColor;stroke-width:1.7;stroke-linecap:round;stroke-linejoin:round}
</style>
