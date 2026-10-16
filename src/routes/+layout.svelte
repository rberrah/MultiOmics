<script>
  import '$lib/styles/theme.css';
  import '@fontsource-variable/inter';
  import '@fontsource-variable/source-serif-4';
  import '@fontsource-variable/jetbrains-mono';
  import {base} from '$app/paths';
  import {language} from '$lib/stores/language';
  let {children}=$props();
  const t=(fr,en)=>$language==='en'?en:fr;
</script>

<a class="skip-link" href="#content">{t('Aller au contenu','Skip to content')}</a>
<header class="site-header">
  <nav aria-label={t('Navigation principale','Main navigation')}>
    <a class="brand" href={base+'/multiomics/'}>MultiOmics</a>
    <div class="links">
      <a href={base+'/multiomics/tool/'}>{t('Analyser','Analyse')}</a>
      <a href={base+'/multiomics/learn/'}>{t('Cours','Learn')}</a>
      <button type="button" class="lang" onclick={()=>language.set($language==='en'?'fr':'en')} aria-label={t('Switch to English','Passer en français')}>{$language==='en'?'FR':'EN'}</button>
    </div>
  </nav>
</header>
<main id="content">{@render children()}</main>
<footer class="site-footer">
  <p>MultiOmics · {t('Outil de recherche exploratoire — résultats R selon disponibilité du backend','Exploratory research tool — native R inference depends on backend availability')}.</p>
  <p>{t('Les calculs exploratoires ne constituent pas une certification de validité biologique ou clinique.','Exploratory calculations do not certify biological or clinical validity.')}</p>
</footer>
<style>
.site-header{background:var(--bg-tertiary);border-bottom:1px solid var(--border-subtle);position:relative;z-index:10}
nav{max-width:1280px;margin:auto;padding:.75rem 1.1rem;display:flex;align-items:center;justify-content:space-between;gap:1rem;flex-wrap:wrap}
.brand{font-size:1.25rem;font-weight:750;text-decoration:none;color:var(--text-primary)}
.links{display:flex;align-items:center;gap:1rem}
.links a{color:var(--text-primary);font-size:.93rem}
.lang{cursor:pointer;padding:.35rem .7rem;background:var(--bg-secondary);border:1px solid var(--border-strong);color:var(--text-primary);border-radius:5px}
main{min-height:65vh}
.site-footer{padding:1.7rem 1rem;margin-top:3rem;text-align:center;border-top:1px solid var(--border-subtle);color:var(--text-secondary);font-size:.84rem}
.site-footer p{margin:.4rem}
@media(max-width:420px){.links{gap:.55rem}.links a{font-size:.83rem}}
</style>
