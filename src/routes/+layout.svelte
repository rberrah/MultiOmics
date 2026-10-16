<script>
  import '$lib/styles/theme.css';
  import '@fontsource-variable/inter';
  import '@fontsource-variable/source-serif-4';
  import '@fontsource-variable/jetbrains-mono';
  import { base } from '$app/paths';
  import { language } from '$lib/stores/language';
  import ThemeToggle from '$lib/components/ThemeToggle.svelte';
  import LanguageToggle from '$lib/components/LanguageToggle.svelte';

  let { children } = $props();
  const t = (fr, en) => $language === 'en' ? en : fr;
  const mainSite = 'https://rberrah.github.io/';
</script>

<a class="skip-link" href="#content">{t('Aller au contenu','Skip to content')}</a>

<header class="site-header">
  <nav class="nav-shell" aria-label={t('Navigation principale','Main navigation')}>
    <div class="identity">
      <a class="brand" href={base + '/multiomics/'}>MultiOmics</a>
      <a class="main-site compact-main-link" href={mainSite}>
        {t('Site principal ↗','Main site ↗')}
      </a>
    </div>

    <div class="primary-links">
      <a href={base + '/multiomics/tool/'}>{t('Analyser','Analyse')}</a>
      <a href={base + '/multiomics/learn/'}>{t('Cours','Learn')}</a>
      <a href={base + '/multiomics/disciplines/'}>{t('Autres disciplines','Other disciplines')}</a>
      <a href={base + '/multiomics/resources/'}>{t('Ressources','Resources')}</a>
    </div>

    <div class="utilities">
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
      <p>{t('Analyse mono- et multi-omique reproductible, reliée à l’écosystème PMx.', 'Reproducible single- and multi-omics analysis connected to the PMx ecosystem.')}</p>
    </div>
    <nav aria-label={t('Liens de pied de page','Footer links')}>
      <a href="https://rberrah.github.io/pharmacometrie/">PMx Explain</a>
      <a href={mainSite}>{t('Site principal','Main site')} ↗</a>
      <a href="https://github.com/rberrah/MultiOmics">GitHub ↗</a>
    </nav>
  </div>
  <p class="fine-print">{t('Outil de recherche exploratoire — les calculs ne constituent pas une certification de validité biologique ou clinique.', 'Exploratory research tool — calculations do not certify biological or clinical validity.')}</p>
</footer>

<style>
  .site-header{background:color-mix(in srgb,var(--bg-tertiary) 94%,transparent);border-bottom:1px solid var(--border-subtle);position:sticky;top:0;z-index:40;backdrop-filter:blur(10px)}
  .nav-shell{width:min(100%,1440px);margin:0 auto;padding:.72rem clamp(1rem,3vw,2.5rem);display:grid;grid-template-columns:auto minmax(0,1fr) auto;align-items:center;gap:1.25rem}
  .identity{display:flex;align-items:baseline;gap:.8rem;white-space:nowrap}
  .brand{font-family:var(--font-heading);font-size:1.3rem;font-weight:700;text-decoration:none;color:var(--text-primary)}
  .main-site{color:var(--text-secondary);font-size:.78rem;text-decoration:none}
  .main-site:hover{color:var(--accent-pk)}
  .primary-links{display:flex;align-items:center;justify-content:center;gap:clamp(.65rem,1.5vw,1.25rem);min-width:0}
  .primary-links a{color:var(--text-primary);font-size:.9rem;text-decoration:none;white-space:nowrap;padding:.35rem .1rem;border-bottom:2px solid transparent}
  .primary-links a:hover{color:var(--accent-pk);border-bottom-color:var(--accent-pk)}
  .utilities{display:flex;align-items:center;gap:.5rem}
  .page-shell{width:min(100%,1440px);margin:0 auto;padding-inline:clamp(1rem,3vw,2.5rem);min-height:65vh}
  .site-footer{margin-top:4rem;border-top:1px solid var(--border-subtle);background:var(--bg-secondary);color:var(--text-secondary)}
  .footer-shell{width:min(100%,1440px);margin:0 auto;padding:1.6rem clamp(1rem,3vw,2.5rem) 1rem;display:flex;justify-content:space-between;gap:2rem;align-items:start}
  .footer-shell strong{color:var(--text-primary)}
  .footer-shell p{margin:.35rem 0 0;max-width:62ch;font-size:.86rem}
  .footer-shell nav{display:flex;flex-wrap:wrap;justify-content:flex-end;gap:.9rem}
  .footer-shell a{font-size:.84rem}
  .fine-print{width:min(100%,1440px);margin:0 auto;padding:0 clamp(1rem,3vw,2.5rem) 1.4rem;font-size:.76rem;color:var(--text-muted)}
  @media(max-width:1050px){
    .nav-shell{grid-template-columns:1fr auto}
    .primary-links{grid-column:1/-1;grid-row:2;justify-content:flex-start;overflow-x:auto;padding-bottom:.15rem}
    .utilities{grid-column:2;grid-row:1}
  }
  @media(max-width:640px){
    .nav-shell{padding-block:.6rem;gap:.55rem .75rem}
    .compact-main-link{display:none}
    .primary-links{gap:1rem;scrollbar-width:none}
    .primary-links::-webkit-scrollbar{display:none}
    .primary-links a{font-size:.82rem}
    .footer-shell{flex-direction:column;gap:1rem}
    .footer-shell nav{justify-content:flex-start}
  }
</style>
