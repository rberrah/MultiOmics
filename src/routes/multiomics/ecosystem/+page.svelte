<script>
  import {language} from '#lib/stores/language.js';
  import ecosystem from '#lib/content/ecosystem.json';
  $: t=$language==='en'?(/** @type {string} */ fr,/** @type {string} */ en)=>en:(/** @type {string} */ fr,/** @type {string} */ en)=>fr;
</script>
<svelte:head><title>{t('Écosystème — Omics Explain','Ecosystem — Omics Explain')}</title></svelte:head>
<section class="connections-page" data-testid="ecosystem-page">
  <h1>{t('L’écosystème Explain','The Explain ecosystem')}</h1>
  <p>{t('Des plateformes complémentaires, accessibles indépendamment. Choisissez selon votre question scientifique.','Complementary platforms with independent entry points. Choose according to your scientific question.')}</p>
  <nav aria-label={t('Plateformes publiques','Public platforms')} class="roles" data-testid="ecosystem-links">
    {#each ecosystem as item}
      <a href={item.url+(item.localized?'?lang='+$language:'')} data-sveltekit-reload aria-current={item.key==='omics'?'true':undefined}>{t(item.fr,item.en)}</a>
    {/each}
  </nav>
  <section class="handoff">
    <h2>{t('Une application des omiques : Omics2PK','One application of omics: Omics2PK')}</h2>
    <p>{t('MultiOmics fournit un échange moléculaire traçable. Une question PK exige aussi des observations PK, des doses, des chronologies et des covariables adaptées chez les mêmes sujets. La compatibilité technique ne démontre pas une valeur ajoutée PK.','MultiOmics provides traceable molecular exchange. A PK question also requires appropriate PK observations, doses, timing and covariates in the same subjects. Technical compatibility does not demonstrate added PK value.')}</p>
    <p><a href={'https://rberrah.github.io/ecosystem/?lang='+$language} data-sveltekit-reload>{t('Voir l’écosystème','Open ecosystem')}</a> · <a href={'https://rberrah.github.io/pharmacometrie/domaines/?lang='+$language} data-sveltekit-reload>{t('Cours Omics → PK/PD','Omics → PK/PD lessons')}</a></p>
  </section>
  <p role="note">{t('Usage recherche et pédagogie : aucune page ne certifie la causalité, la validité clinique ou une recommandation de dose.','Research and education: no page certifies causality, clinical validity or a dosing recommendation.')}</p>
</section>
<style>.connections-page{padding-block:clamp(2rem,6vw,4rem)}.roles{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,210px),1fr));gap:1rem;margin:2rem 0}.roles a{padding:1rem;min-height:44px;border:1px solid var(--border-strong);border-radius:8px}.roles a[aria-current]{border:2px solid var(--accent-pk)}.handoff{padding:1rem;border:1px solid var(--border-strong);border-radius:8px}</style>
