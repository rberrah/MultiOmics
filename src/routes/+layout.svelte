<script>
  import '$lib/styles/global.css';
  import { base } from '$app/paths';
  import { language } from '$lib/stores/language';
  let { children } = $props();
  let dark = $state(false);
  function toggleTheme() { dark = !dark; document.documentElement.dataset.theme = dark ? 'dark' : 'light'; }
</script>
<svelte:head><meta name="color-scheme" content="light dark" /></svelte:head>
<a class="skip-link" href="#main-content">{$language === 'en' ? 'Skip to content' : 'Aller au contenu'}</a>
<header class="standalone-header">
  <a class="brand" href={base + '/multiomics/'}>MultiOmics <small>Research</small></a>
  <nav aria-label={$language === 'en' ? 'Main navigation' : 'Navigation principale'}>
    <a href={base + '/multiomics/'}>{$language === 'en' ? 'Overview' : 'Présentation'}</a>
    <a href={base + '/multiomics/tool/'}>{$language === 'en' ? 'Analyze' : 'Analyser'}</a>
    <a href={base + '/multiomics/learn/'}>{$language === 'en' ? 'Learn' : 'Comprendre'}</a>
    <button type="button" aria-label="Switch language" onclick={() => language.set($language === 'en' ? 'fr' : 'en')}>{$language === 'en' ? 'FR' : 'EN'}</button>
    <button type="button" aria-label={$language === 'en' ? 'Toggle theme' : 'Changer le thème'} onclick={toggleTheme}>{dark ? '☀' : '◐'}</button>
  </nav>
</header>
<main id="main-content" class="standalone-main">{@render children()}</main>
<footer class="standalone-footer">MultiOmics · {$language === 'en' ? 'Exploratory research software; statistical validation depends on method and study design.' : 'Logiciel de recherche exploratoire ; validation statistique spécifique à chaque méthode et plan d’étude.'}
<a href="https://github.com/rberrah/MultiOmics">Source GitHub</a></footer>
