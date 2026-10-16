<script>
  import { language } from '$lib/stores/language';
  /** @param {string} fr @param {string} en */
  const t=(fr,en)=>$language==='en'?en:fr;
</script>

<svelte:head>
  <title>{t('Autres disciplines — MultiOmics','Other disciplines — MultiOmics')}</title>
  <meta name="description" content={t('Relier l’analyse omique à la pharmacométrie avec PMx Explain, puis à terme Omics2PK.','Connect omics analysis to pharmacometrics through PMx Explain and, later, Omics2PK.')} />
</svelte:head>

<section class="connections-page">
  <p class="eyebrow">{t('Liens avec les autres disciplines','Connections with other disciplines')}</p>
  <h1>{t('MultiOmics ne s’arrête pas à l’analyse omique.','MultiOmics does not stop at omics analysis.')}</h1>
  <p class="lede">{t('Cette page explicite les outils complémentaires sans mélanger leurs niveaux de preuve ni leurs responsabilités scientifiques.', 'This page makes complementary tools explicit without mixing their evidence levels or scientific responsibilities.')}</p>

  <article class="connection primary">
    <div>
      <p class="eyebrow">PMx Explain</p>
      <h2>{t('Comprendre comment un signal biologique peut entrer dans un raisonnement PK/PD.', 'Understand how a biological signal can enter PK/PD reasoning.')}</h2>
      <p>{t('PMx Explain est la plateforme pédagogique de pharmacométrie : PK, PD, PopPK, covariables, MIPD et cours dédiés aux liens entre omiques et PK/PD.', 'PMx Explain is the pharmacometrics learning platform: PK, PD, PopPK, covariates, MIPD and dedicated lessons connecting omics with PK/PD.')}</p>
      <div class="actions">
        <a class="btn btn-primary" href="https://rberrah.github.io/pharmacometrie/">{t('Ouvrir PMx Explain','Open PMx Explain')} ↗</a>
        <a class="btn btn-outline" href="https://rberrah.github.io/pharmacometrie/domaines/">{t('Cours Omics → PK/PD','Omics → PK/PD lessons')} ↗</a>
      </div>
    </div>
    <div class="flow" aria-label={t('Chaîne conceptuelle','Conceptual chain')}>
      <span>Omics</span><b>→</b><span>{t('Hypothèse biologique','Biological hypothesis')}</span><b>→</b><span>PK / PD</span>
    </div>
  </article>

  <article class="connection future">
    <p class="eyebrow">Omics2PK · {t('à venir','planned')}</p>
    <h2>{t('Un pont logiciel dédié sera ajouté plus tard.', 'A dedicated software bridge will be added later.')}</h2>
    <p>{t('Omics2PK aura pour rôle de traduire des résultats omiques en hypothèses pharmacométriques candidates de façon traçable. Il n’est pas présenté comme un outil disponible tant que son dépôt et ses contrats ne sont pas qualifiés.', 'Omics2PK will translate omics results into traceable candidate pharmacometric hypotheses. It is not presented as an available tool until its repository and contracts are qualified.')}</p>
  </article>
</section>

<style>
  .connections-page{padding-block:clamp(2.5rem,6vw,5rem)}
  h1{max-width:25ch}
  .lede{max-width:74ch;color:var(--text-secondary);font-size:var(--text-lg)}
  .connection{margin-top:2rem;padding:clamp(1.25rem,3vw,2rem);border:1px solid var(--border-subtle);border-radius:var(--radius-lg);background:var(--bg-tertiary)}
  .connection.primary{display:grid;grid-template-columns:minmax(0,1.3fr) minmax(260px,.7fr);gap:2rem;align-items:center;border-top:3px solid var(--accent-pk)}
  .connection p{color:var(--text-secondary)}
  .actions{display:flex;flex-wrap:wrap;gap:.7rem;margin-top:1.3rem}
  .flow{display:grid;grid-template-columns:1fr auto 1fr auto 1fr;align-items:center;gap:.55rem;text-align:center}
  .flow span{padding:.9rem .7rem;border:1px solid var(--border-strong);background:var(--bg-secondary);border-radius:var(--radius)}
  .flow b{color:var(--accent-pk)}
  .future{border-style:dashed;background:var(--bg-secondary)}
  @media(max-width:820px){.connection.primary{grid-template-columns:1fr}.flow{grid-template-columns:1fr}.flow b{transform:rotate(90deg)}}
</style>
