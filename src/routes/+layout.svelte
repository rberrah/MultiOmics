<script>
  import '$lib/styles/theme.css';
  import '@fontsource-variable/inter';
  import '@fontsource-variable/source-serif-4';
  import '@fontsource-variable/jetbrains-mono';
  import { base } from '$app/paths';
  import { page } from '$app/stores';
  import { afterNavigate } from '$app/navigation';
  import { language } from '$lib/stores/language';
  import ThemeToggle from '$lib/components/ThemeToggle.svelte';
  import LanguageToggle from '$lib/components/LanguageToggle.svelte';

  let { children } = $props();

  /** @param {string} fr @param {string} en */
  const t = (fr, en) => $language === 'en' ? en : fr;
  const mainSite = 'https://rberrah.github.io/';
  const pmxExplain = 'https://rberrah.github.io/pharmacometrie/';
  /** @type {HTMLElement | undefined} */
  let navigation = $state();
  let menuOpen = $state(false);
  let links = $derived([
    { key: 'analyze', path: '/multiomics/tool/', label: t('Analyser', 'Analyze') },
    { key: 'learn', path: '/multiomics/learn/', label: t('Cours', 'Learn') },
    { key: 'disciplines', path: '/multiomics/disciplines/', label: t('Autres disciplines', 'Other disciplines') },
    { key: 'resources', path: '/multiomics/resources/', label: t('Ressources', 'Resources') }
  ]);

  /** @param {string} path */
  const isActive = (path) => $page.url.pathname.startsWith(base + path);
  function closeNavigation() { menuOpen = false; }
  /** @param {MouseEvent} event */
  function handleOutsideClick(event) {
    if (!menuOpen || !(event.target instanceof Node)) return;
    if (navigation?.contains(event.target)) return;
    if (event.target instanceof Element && event.target.closest('.burger')) return;
    closeNavigation();
  }
  /** @param {KeyboardEvent} event */
  function handleEscape(event) {
    if (event.key === 'Escape' && menuOpen) {
      closeNavigation();
      const toggle = document.querySelector('.burger');
      if (toggle instanceof HTMLElement) toggle.focus();
    }
  }

  afterNavigate(closeNavigation);
  $effect(() => { document.documentElement.lang = $language === 'en' ? 'en' : 'fr'; });
</script>

<svelte:window onclick={handleOutsideClick} onkeydown={handleEscape} />
<a class="skip-link" href="#content">{t('Aller au contenu', 'Skip to content')}</a>

<div class="app">
  <header class="site-header" data-testid="site-header">
    <div class="header-shell">
      <a class="logo" href={base + '/multiomics/'} aria-label="MultiOmics" data-testid="logo-link">
        <span class="mark" aria-hidden="true">Om</span>
        <span class="word" aria-hidden="true">Multi<em>Omics</em></span>
        <span class="tagline" aria-hidden="true">Omics → PK/PD</span>
      </a>

      <button
        class="burger"
        type="button"
        aria-label={t('Ouvrir le menu', 'Open menu')}
        aria-expanded={menuOpen}
        aria-controls="site-nav"
        onclick={() => (menuOpen = !menuOpen)}
        data-testid="nav-toggle"
      ><span></span><span></span><span></span></button>

      <nav
        id="site-nav"
        class:open={menuOpen}
        bind:this={navigation}
        aria-label={t('Navigation principale', 'Main navigation')}
        data-testid="primary-nav"
      >
        {#each links as link}
          <a
            class:active={isActive(link.path)}
            aria-current={isActive(link.path) ? 'page' : undefined}
            href={base + link.path}
            onclick={closeNavigation}
            data-testid={'nav-' + link.key}
          >{link.label}</a>
        {/each}
        <a
          class="ecosystem-link"
          href={pmxExplain}
          target="_blank"
          rel="noopener noreferrer"
          data-testid="nav-pmx"
        >PMx Explain ↗</a>
        <a
          class="ecosystem-link main-site-link"
          href={mainSite}
          target="_blank"
          rel="noopener noreferrer"
          data-testid="nav-main-site"
        >{t('Site principal', 'Main site')} ↗</a>
        <div class="tools" data-testid="header-controls">
          <LanguageToggle />
          <ThemeToggle />
        </div>
      </nav>
    </div>
  </header>

  <main id="content" class="page-shell">
    {@render children()}
  </main>

  <footer class="site-footer">
    <div class="footer-shell">
      <div>
        <strong>MultiOmics</strong>
        <p>{t('Analyse mono- et multi-omique reproductible, reliée à l’écosystème PMx.', 'Reproducible single- and multi-omics analysis connected to the PMx ecosystem.')}</p>
      </div>
      <nav aria-label={t('Liens de pied de page', 'Footer links')}>
        <a href={pmxExplain}>PMx Explain</a>
        <a href={mainSite}>{t('Site principal', 'Main site')} ↗</a>
        <a href="https://github.com/rberrah/MultiOmics">GitHub ↗</a>
      </nav>
    </div>
  </footer>

  <div class="disclaimer" role="note">
    {t(
      'Outil de recherche exploratoire — les calculs ne constituent pas une certification de validité biologique ou clinique.',
      'Exploratory research tool — calculations do not certify biological or clinical validity.'
    )}
  </div>
</div>

<style>
  .app { min-height: 100vh; display: flex; flex-direction: column; }
  .site-header {
    position: sticky; top: 0; z-index: 50;
    min-height: var(--header-h);
    background: color-mix(in srgb, var(--bg-primary) 85%, transparent);
    backdrop-filter: blur(14px);
    border-bottom: 1px solid var(--border-subtle);
  }
  .header-shell {
    width: 100%; max-width: var(--maxw); margin: 0 auto;
    padding: var(--space-3) var(--space-6);
    display: flex; justify-content: space-between; align-items: center;
    gap: var(--space-4);
  }
  .logo {
    display: flex; align-items: center; gap: var(--space-3);
    text-decoration: none; color: var(--text-primary); white-space: nowrap;
    flex-shrink: 0;
  }
  .mark {
    font-family: var(--font-heading); font-weight: 900; color: var(--bg-primary);
    background: var(--accent-pk); width: 34px; height: 34px;
    border-radius: 8px; display: grid; place-items: center;
    font-size: .95rem; letter-spacing: -.04em;
  }
  .word {
    font-family: var(--font-heading); font-size: 1.05rem; font-weight: 700;
    letter-spacing: -.02em;
  }
  .word em { color: var(--accent-pk); font-style: normal; }
  .tagline {
    font-family: var(--font-mono); font-size: .64rem; line-height: 1.2;
    color: var(--text-muted); border-left: 1px solid var(--border-strong);
    padding-left: var(--space-3); letter-spacing: -.02em;
  }
  nav#site-nav { display: flex; align-items: center; gap: var(--space-1); min-width: 0; }
  nav#site-nav > a {
    font-size: var(--text-sm); text-decoration: none; color: var(--text-secondary);
    padding: var(--space-2) var(--space-3); border-radius: 6px; font-weight: 600;
    white-space: nowrap; transition: color .18s ease, background-color .18s ease;
  }
  nav#site-nav > a:hover { color: var(--text-primary); background: var(--bg-secondary); }
  nav#site-nav > a.active {
    color: var(--accent-pk);
    background: color-mix(in srgb, var(--accent-pk) 10%, transparent);
  }
  nav#site-nav > a.ecosystem-link { font-weight: 500; }
  .tools { display: flex; gap: var(--space-2); align-items: center; margin-left: var(--space-2); }
  .burger {
    display: none; flex-direction: column; gap: 4px;
    background: none; border: 1px solid var(--border-subtle);
    border-radius: var(--radius); cursor: pointer; padding: 8px;
    color: var(--text-primary);
  }
  .burger span { width: 22px; height: 2px; background: currentColor; display: block; }
  .page-shell {
    flex: 1; width: 100%; max-width: var(--maxw);
    margin: 0 auto; padding: var(--space-8) var(--space-6) var(--space-24);
  }
  .site-footer {
    border-top: 1px solid var(--border-subtle);
    background: var(--bg-secondary); color: var(--text-secondary);
  }
  .footer-shell {
    max-width: var(--maxw); margin: 0 auto;
    padding: var(--space-8) var(--space-6);
    display: flex; justify-content: space-between; align-items: start; gap: var(--space-8);
  }
  .footer-shell strong { color: var(--text-primary); }
  .footer-shell p { font-size: var(--text-sm); max-width: 62ch; margin: var(--space-2) 0 0; }
  .footer-shell nav { display: flex; flex-wrap: wrap; justify-content: flex-end; gap: var(--space-4); }
  .footer-shell a { font-size: var(--text-sm); }
  .disclaimer {
    text-align: center; padding: 6px var(--space-4);
    background: var(--disclaimer-bg); color: var(--disclaimer-text);
    font-family: var(--font-mono); font-size: var(--text-xs);
  }
  .logo:focus-visible, .burger:focus-visible, nav#site-nav a:focus-visible {
    outline: 2px solid var(--focus-ring); outline-offset: 3px;
  }
  @media(max-width: 1180px) {
    .tagline { display: none; }
  }
  @media(max-width: 1050px) {
    .header-shell { flex-wrap: wrap; }
    .burger { display: flex; margin-left: auto; }
    nav#site-nav {
      display: none; flex-direction: column; align-items: stretch;
      flex-basis: 100%; width: 100%; gap: var(--space-1);
      max-height: calc(100dvh - var(--header-h)); overflow-y: auto;
      padding: var(--space-2) 0 var(--space-3);
      border-top: 1px solid var(--border-subtle);
    }
    nav#site-nav.open { display: flex; }
    nav#site-nav > a { padding: var(--space-3); white-space: normal; }
    .tools { margin: var(--space-2) 0 0; }
  }
  @media(max-width: 640px) {
    .header-shell { padding: var(--space-3) var(--space-4); }
    .page-shell { padding: var(--space-6) var(--space-4) var(--space-24); }
    .footer-shell { flex-direction: column; gap: var(--space-4); padding: var(--space-8) var(--space-4); }
    .footer-shell nav { justify-content: flex-start; }
  }
</style>
