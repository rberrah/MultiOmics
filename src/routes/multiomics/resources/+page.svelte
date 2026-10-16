<script>
  import { base } from '$app/paths';
  import { language } from '$lib/stores/language';
  const t=(fr,en)=>$language==='en'?en:fr;
</script>

<svelte:head>
  <title>{t('Ressources — MultiOmics','Resources — MultiOmics')}</title>
  <meta name="description" content={t('Documentation, formats, reproductibilité, code source et citation de MultiOmics.','Documentation, file formats, reproducibility, source code and citation resources for MultiOmics.')} />
</svelte:head>

<section class="resource-page">
  <p class="eyebrow">MultiOmics · {t('Ressources','Resources')}</p>
  <h1>{t('Retrouver les méthodes, formats et preuves de reproductibilité.','Find methods, formats and reproducibility evidence.')}</h1>
  <p class="lede">{t('Les ressources restent séparées de l’outil d’analyse afin de garder l’interface débutant simple tout en permettant un audit scientifique complet.','Resources stay separate from the analysis interface so the beginner workflow remains simple while full scientific auditing stays available.')}</p>

  <div class="resource-grid">
    <article class="card">
      <h2>{t('Méthodes et interprétation','Methods and interpretation')}</h2>
      <p>{t('Cours bilingue sur le plan d’étude, le QC, les tests, l’intégration et les limites.', 'Bilingual guidance on study design, QC, statistical testing, integration and limitations.')}</p>
      <a href={base + '/multiomics/learn/'}>{t('Ouvrir le cours','Open learning resources')} →</a>
    </article>
    <article class="card">
      <h2>{t('Formats de fichiers','File formats')}</h2>
      <p>{t('Structure attendue des matrices et du tableau des échantillons.', 'Expected structure for matrices and sample metadata.')}</p>
      <a href={base + '/multiomics/README.txt'} target="_blank">{t('Voir les formats','View formats')} ↗</a>
    </article>
    <article class="card">
      <h2>{t('Code et versions','Code and releases')}</h2>
      <p>{t('Dépôt autonome, historique de développement et workflows de qualification.', 'Standalone repository, development history and qualification workflows.')}</p>
      <a href="https://github.com/rberrah/MultiOmics">GitHub ↗</a>
    </article>
    <article class="card">
      <h2>{t('Citation et réutilisation','Citation and reuse')}</h2>
      <p>{t('Métadonnées de citation, licence non commerciale et règles d’attribution.', 'Citation metadata, noncommercial licensing and attribution rules.')}</p>
      <a href="https://github.com/rberrah/MultiOmics/blob/main/CITATION.cff">CITATION.cff ↗</a>
    </article>
    <article class="card">
      <h2>{t('PMx Explain','PMx Explain')}</h2>
      <p>{t('Cours et outils de pharmacométrie reliés à cet environnement omique.', 'Pharmacometrics courses and tools connected to this omics environment.')}</p>
      <a href="https://rberrah.github.io/pharmacometrie/">{t('Ouvrir PMx Explain','Open PMx Explain')} ↗</a>
    </article>
    <article class="card">
      <h2>{t('Site principal','Main site')}</h2>
      <p>{t('Publications, communications et autres outils de recherche.', 'Publications, communications and other research tools.')}</p>
      <a href="https://rberrah.github.io/">{t('Retour au site','Back to main site')} ↗</a>
    </article>
  </div>
</section>

<style>
  .resource-page{padding-block:clamp(2.5rem,6vw,5rem)}
  h1{max-width:24ch}
  .lede{max-width:72ch;color:var(--text-secondary);font-size:var(--text-lg)}
  .resource-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:1rem;margin-top:2.3rem}
  .card{min-height:210px}
  .card h2{font-size:var(--text-xl);margin-top:0}
  .card p{color:var(--text-secondary)}
  @media(max-width:950px){.resource-grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
  @media(max-width:620px){.resource-grid{grid-template-columns:1fr}}
</style>
