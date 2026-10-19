<script>
  import '#lib/styles/theme.css';
  import '@fontsource-variable/inter/index.css';
  import '@fontsource-variable/source-serif-4/index.css';
  import '@fontsource-variable/jetbrains-mono/index.css';
  import { BookOpen, ChartNoAxesCombined, Library, ChevronDown } from '@lucide/svelte';
  import { resolve } from '$app/paths';
  import { page } from '$app/state';
  import { afterNavigate } from '$app/navigation';
  import { language } from '#lib/stores/language.js';
  import ThemeToggle from '#lib/components/ThemeToggle.svelte';
  import LanguageToggle from '#lib/components/LanguageToggle.svelte';

  let { children } = $props();
  /** @param {string} fr @param {string} en */
  const t = (fr,en) => $language === 'en' ? en : fr;
  const mainSite = 'https://rberrah.github.io/';
  const pmxExplain = 'https://rberrah.github.io/pharmacometrie/';
  /** @type {HTMLElement | undefined} */
  let navigation = $state();
  let menuOpen = $state(false);

  // The actual PMx Explain nav pattern: goal-based native <details> dropdowns.
  // Local routes remain under /multiomics/ in both BASE_PATH builds.
  let links = $derived([
    { key:'overview',href:'/multiomics/',label:t('Accueil','Home') },
    { key:'learn',href:'/multiomics/learn/',label:t('Cours et méthodes','Courses and methods') },
    { key:'design',href:'/multiomics/learn/',label:t('Plan d’étude et qualité','Study design and QC') },
    { key:'tool',href:'/multiomics/tool/',label:t('Lancer une analyse','Run an analysis') },
    { key:'workflow',href:'/multiomics/tool/',label:t('Découvrir en 1 clic','One-click demos') },
    { key:'disciplines',href:'/multiomics/disciplines/',label:t('Vue d’ensemble de l’écosystème','Ecosystem overview') },
    { key:'pmx',href:pmxExplain,label:'PMx Explain ↗',external:true },
    { key:'omics2pk',href:'https://rberrah.github.io/Omics2PK/',label:'Omics2PK ↗',external:true },
    { key:'resources',href:'/multiomics/resources/',label:t('Documentation','Documentation') },
    { key:'glossary',href:'/multiomics/glossary/',label:t('Glossaire','Glossary') },
    { key:'references',href:'/multiomics/references/',label:t('Références scientifiques','Scientific references') },
    { key:'source',href:'https://github.com/rberrah/MultiOmics',label:t('Code source et versions ↗','Source code and versions ↗'),external:true },
    { key:'main-site',href:mainSite,label:'Racym Berrah ↗',external:true }
  ]);
  let groups = $derived([
    {key:'learn',label:t('Apprendre','Learn'),icon:BookOpen,keys:['learn','design']},
    {key:'analyze',label:t('Analyser','Analyze'),icon:ChartNoAxesCombined,keys:['tool','workflow']},
    {key:'disciplines',label:t('Écosystème','Ecosystem'),icon:Library,keys:['disciplines','pmx','omics2pk']},
    {key:'resources',label:t('Ressources','Resources'),icon:Library,keys:['resources','glossary','references','source','main-site']}
  ]);
  // SvelteKit 3 only types concrete route IDs; menu hrefs are dynamic.
  // The prefix is resolved from an actual prerendered route, keeping both
  // standalone and PMx legacy base paths intact.
  /** @param {string} href */
  const internalPath = href => resolve('multiomics/') + href.slice('/multiomics/'.length);
  /** @param {string} href */
  const isActive = (href) => {
    if ((/^https?:/).test(href)) return false;
    const pathname = page.url.pathname;
    const route = internalPath(href);
    return href === '/multiomics/' ? pathname === route : pathname.startsWith(route);
  };
  function closeNavigation() {
    menuOpen = false;
    navigation?.querySelectorAll('details[open]').forEach((item) => item.removeAttribute('open'));
  }
  /** @param {MouseEvent} event */
  function handleOutsideClick(event) {
    if (navigation?.contains(event.target instanceof Node ? event.target : null)) return;
    if (event.target instanceof Element && event.target.closest('.burger')) return;
    closeNavigation();
  }
  /** @param {KeyboardEvent} event */
  function handleEscape(event) {
    if (event.key !== 'Escape') return;
    const opened = navigation?.querySelector('details[open] summary');
    if (opened instanceof HTMLElement) opened.focus();
    else if (menuOpen) {
      const burger = document.querySelector('.burger');
      if (burger instanceof HTMLElement) burger.focus();
    }
    closeNavigation();
  }

  afterNavigate(closeNavigation);
  $effect(() => { document.documentElement.lang = $language === 'en' ? 'en' : 'fr'; });
</script>

<svelte:window
  onclick={handleOutsideClick}
  onkeydown={handleEscape}
></svelte:window>

<a class="skip-link" href="#content">{t('Aller au contenu', 'Skip to content')}</a>

<div class="app">
  <header class="site-header" data-testid="site-header">
    <a
      class="logo"
      href={resolve('multiomics/')}
      aria-label="MultiOmics"
      data-testid="logo-link"
    >
      <span class="mark" aria-hidden="true">mO</span>
      <span class="word">Multi<em>Omics</em></span>
    </a>

    <button
      class="burger"
      type="button"
      aria-label={t('Ouvrir le menu', 'Open menu')}
      aria-expanded={menuOpen}
      aria-controls="site-nav"
      onclick={() => menuOpen = !menuOpen}
      data-testid="nav-toggle"
    ><span></span><span></span><span></span></button>

    <nav
      id="site-nav"
      bind:this={navigation}
      class:open={menuOpen}
      aria-label={t('Navigation principale', 'Main navigation')}
      data-testid="primary-nav"
    >
      {#each groups as group}
        <details
          class="goal"
          name="site-objective"
          data-testid={'goal-' + group.key}
        >
          <summary
            class:active={links.some((link) => group.keys.includes(link.key) && isActive(link.href))}
          >
            <group.icon size={16} />
            {group.label}
            <ChevronDown size={13} />
          </summary>
          <div class="goal-links">
            {#each group.keys as key}
              {@const link = links.find((item) => item.key === key)}
              {#if link}
                <a class:active={isActive(link.href)}
                  aria-current={isActive(link.href) ? 'page' : undefined}
                  href={link.external ? link.href : internalPath(link.href)}
                  target={link.external ? '_blank' : undefined}
                  rel={link.external ? 'noopener noreferrer' : undefined}
                  onclick={closeNavigation} data-testid={'nav-' + link.key}>{link.label}</a>
              {/if}
            {/each}
          </div>
        </details>
      {/each}
      <div class="tools" data-testid="header-controls">
        <LanguageToggle />
        <ThemeToggle />
      </div>
    </nav>
  </header>

  <main id="content" class="page-shell">
    {@render children()}
  </main>

  <footer class="site-footer">
    <div class="footer-shell">
      <div>
        <strong>MultiOmics</strong>
        <p>{t('Analyser les omiques dans un écosystème de recherche partagé avec PMx Explain et Omics2PK.', 'Analyze omics within a research ecosystem shared with PMx Explain and Omics2PK.')}</p>
      </div>
      <nav aria-label={t('Liens de pied de page', 'Footer links')}>
        <a href={pmxExplain}>PMx Explain</a>
        <a href="https://rberrah.github.io/Omics2PK/">Omics2PK</a>
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
  /* PMx Explain's editorial header/layout/dropdown geometry (same tokens). */
  .app { min-height: 100vh; display: flex; flex-direction: column; }
  .site-header {
    display: flex; align-items: center; justify-content: space-between;
    gap: var(--space-4); padding: var(--space-3) var(--space-6);
    min-height: var(--header-h);
    position: sticky; top: 0; z-index: 50;
    background: color-mix(in srgb, var(--bg-primary) 85%, transparent);
    backdrop-filter: blur(14px);
    border-bottom: 1px solid var(--border-subtle);
  }
  .logo { display: flex; align-items: center; gap: var(--space-3); text-decoration: none; color: var(--text-primary); white-space: nowrap; }
  .mark {
    font-family: var(--font-heading); font-weight: 900; color: var(--bg-primary);
    background: var(--accent-pk); width: 34px; height: 34px; border-radius: 8px;
    display: grid; place-items: center; font-size: .95rem; letter-spacing: -.04em;
  }
  .word { font-family: var(--font-heading); font-weight: 700; font-size: 1.05rem; letter-spacing: -.02em; }
  .word em { font-style: normal; color: var(--accent-pk); }
  nav#site-nav { display: flex; gap: var(--space-1); align-items: center; }
  .goal { position: relative; margin: 0; }
  .goal summary { list-style: none; display: flex; align-items: center; gap: 7px; padding: 10px; font-size: 13px; font-weight: 600; cursor: pointer; color: var(--text-secondary); white-space: nowrap; }
  .goal summary::-webkit-details-marker { display: none; }
  .goal summary.active { color: var(--accent-pk); }
  .goal[open] summary, .goal summary:hover { background: var(--bg-secondary); color: var(--text-primary); }
  .goal-links { position: absolute; top: 100%; left: 0; min-width: 220px; display: flex; flex-direction: column; background: var(--bg-primary); border: 1px solid var(--border-strong); padding: 6px; box-shadow: 0 8px 18px #0001; }
  .goal:last-of-type .goal-links { left: auto; right: 0; }
  .goal-links a { white-space: normal; }
  .tools { display: flex; align-items: center; gap: var(--space-2); margin-left: var(--space-2); }
  nav#site-nav :global(.language-toggle) { margin-left: 0; }
  nav#site-nav a {
    text-decoration: none; color: var(--text-secondary);
    padding: var(--space-2) var(--space-3); border-radius: 6px; font-weight: 600;
    font-size: var(--text-sm); transition: background-color .2s ease, color .2s ease;
  }

  nav#site-nav a:hover {
    color: var(--text-primary);
    background: var(--bg-secondary);
  }

  nav#site-nav a.active {
    color: var(--accent-pk);
    background: rgba(184,92,56,.1);
  }

  .burger {
    display: none;
    flex-direction: column;
    gap: 4px;
    background: none;
    border: none;
    cursor: pointer;
    padding: 8px;
  }

  .burger span {
    width: 22px;
    height: 2px;
    background: var(--text-primary);
    display: block;
  }

  .page-shell {
    flex: 1;
    width: 100%;
    max-width: var(--maxw);
    margin: 0 auto;
    padding: var(--space-8) var(--space-6) var(--space-24);
  }

  .site-footer {
    border-top: 1px solid var(--border-subtle);
    background: var(--bg-secondary);
    color: var(--text-secondary);
  }

  .footer-shell {
    max-width: var(--maxw);
    margin: 0 auto;
    padding: var(--space-8) var(--space-6);
    display: flex;
    justify-content: space-between;
    align-items: start;
    gap: var(--space-8);
  }

  .footer-shell strong {
    color: var(--text-primary);
  }

  .footer-shell p {
    font-size: var(--text-sm);
    max-width: 62ch;
    margin: var(--space-2) 0 0;
  }

  .footer-shell nav {
    display: flex;
    flex-wrap: wrap;
    justify-content: flex-end;
    gap: var(--space-4);
  }

  .footer-shell a {
    font-size: var(--text-sm);
  }

  .disclaimer {
    text-align: center;
    padding: 6px var(--space-4);
    background: var(--disclaimer-bg);
    color: var(--disclaimer-text);
    font-family: var(--font-mono);
    font-size: var(--text-xs);
  }

  .logo:focus-visible,
  .burger:focus-visible,
  .goal summary:focus-visible,
  nav#site-nav a:focus-visible {
    outline: 2px solid var(--focus-ring);
    outline-offset: 3px;
  }

  @media (max-width: 1100px) {
    .burger { display: flex; }
    nav#site-nav { position: absolute; top: 100%; right: 0; left: 0; flex-direction: column; align-items: stretch;
      background: var(--bg-primary); border-bottom: 1px solid var(--border-subtle);
      padding: var(--space-3) var(--space-6); display: none;
      max-height: calc(100dvh - var(--header-h)); overflow-y: auto; }
    nav#site-nav.open { display: flex; }
    .goal-links, .goal:last-of-type .goal-links { position: static; min-width: 0; box-shadow: none; border: 0; border-left: 2px solid var(--border-strong); margin-left: 18px; }
    .tools { margin: var(--space-2) 0 0; }
    nav#site-nav :global(.language-toggle) { width: max-content; }
  }
  @media (max-width: 760px) {
    .page-shell { padding: var(--space-6) var(--space-4) var(--space-24); }
    .footer-shell { flex-direction: column; gap: var(--space-4); padding: var(--space-8) var(--space-4); }
    .footer-shell nav { justify-content: flex-start; }
  }
  @media (max-width: 390px) {
    .site-header { padding-inline: var(--space-4); }
    .word { font-size: .98rem; }
  }
</style>
