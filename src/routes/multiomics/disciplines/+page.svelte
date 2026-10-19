<script>
  import { asset } from '$app/paths';
  import { language } from '#lib/stores/language.js';
  /** @param {string} fr @param {string} en */
  const t=(fr,en)=>$language==='en'?en:fr;
</script>

<svelte:head>
  <title>{t('Écosystème — MultiOmics','Ecosystem — MultiOmics')}</title>
  <meta name="description" content={t(
    'MultiOmics s’inscrit dans un écosystème cohérent avec PMx Explain et Omics2PK : apprendre la pharmacométrie, analyser les omiques puis tester leur apport pour la PK/PD.',
    'MultiOmics belongs to a coherent ecosystem with PMx Explain and Omics2PK: learn pharmacometrics, analyze omics, then test their value for PK/PD research.'
  )} />
</svelte:head>

<section class="connections-page" data-testid="ecosystem-page">
  <header class="intro">
    <p class="eyebrow">{t('ÉCOSYSTÈME DE RECHERCHE','RESEARCH ECOSYSTEM')}</p>
    <h1>{t('Analyser les omiques sans perdre le lien avec la pharmacométrie.','Analyze omics without losing the pharmacometrics context.')}</h1>
    <p class="lede">{t(
      'MultiOmics est la brique « Analyser ». PMx Explain apporte le cadre pédagogique et Omics2PK teste ensuite, séparément, si des informations moléculaires individuelles apportent réellement quelque chose à une question PK/PD.',
      'MultiOmics is the Analyze layer. PMx Explain provides the learning framework, while Omics2PK separately tests whether individual molecular information genuinely adds value to a PK/PD question.'
    )}</p>
  </header>

  <figure class="ecosystem-figure">
    <img src={asset('/multiomics/ecosystem-graphical-abstract.svg')} alt={t(
      'Écosystème : PMx Explain pour apprendre, MultiOmics pour analyser, Omics2PK pour intégrer les signaux moléculaires à la recherche PK.',
      'Ecosystem: PMx Explain for learning, MultiOmics for analysis, Omics2PK for integrating molecular signals into PK research.'
    )} />
    <figcaption>{t(
      'Même logique éditoriale, dépôts séparés et frontières de preuve explicites.',
      'Shared editorial logic, separate repositories and explicit evidence boundaries.'
    )}</figcaption>
  </figure>

  <div class="roles">
    <article class="connection">
      <p class="role">{t('01 · Apprendre','01 · Learn')}</p>
      <h2>PMx Explain</h2>
      <p>{t(
        'Cours, animations et outils pour comprendre PK/PD, PopPK, covariables, MIPD et les hypothèses nécessaires avant d’utiliser un signal biologique dans un modèle.',
        'Courses, animations and tools for understanding PK/PD, PopPK, covariates, MIPD and the assumptions needed before using a biological signal in a model.'
      )}</p>
      <div class="actions">
        <a class="btn btn-outline" href="https://rberrah.github.io/pharmacometrie/ecosystem/">{t('Voir l’écosystème','Open ecosystem')} ↗</a>
        <a class="btn btn-outline" href="https://rberrah.github.io/pharmacometrie/domaines/">{t('Cours Omics → PK/PD','Omics → PK/PD lessons')} ↗</a>
      </div>
    </article>

    <article class="connection current">
      <p class="role">{t('02 · Analyser','02 · Analyze')}</p>
      <h2>MultiOmics</h2>
      <p>{t(
        'QC, analyses mono- et multi-omiques, méthodes natives R lorsque qualifiées, niveaux de preuve visibles et export reproductible. MultiOmics ne transforme pas un résultat de groupe en covariable patient.',
        'QC, single- and multi-omics analyses, qualified native-R methods, visible evidence levels and reproducible export. MultiOmics does not turn a group-level result into a patient covariate.'
      )}</p>
      <a class="btn btn-primary" href="../tool/">{t('Lancer une analyse','Run an analysis')} →</a>
    </article>

    <article class="connection">
      <p class="role">{t('03 · Intégrer','03 · Integrate')}</p>
      <h2>Omics2PK</h2>
      <p>{t(
        'Consomme un échange individuel traçable et le confronte à des données cliniques, pharmacogénétiques et PK authentiques. La compatibilité technique n’est jamais présentée comme une preuve de bénéfice PK.',
        'Consumes a traceable individual exchange and confronts it with genuine clinical, pharmacogenetic and PK observations. Technical compatibility is never presented as evidence of PK benefit.'
      )}</p>
      <a class="btn btn-outline" href="https://rberrah.github.io/Omics2PK/">{t('Ouvrir Omics2PK','Open Omics2PK')} ↗</a>
    </article>
  </div>

  <section class="handoff">
    <div>
      <p class="eyebrow">{t('FRONTIÈRE MULTIOMICS → OMICS2PK','MULTIOMICS → OMICS2PK BOUNDARY')}</p>
      <h2>{t('Échange structurel vérifié ≠ donnée PK prête.','Verified structural exchange ≠ PK-ready data.')}</h2>
      <p>{t(
        'Un export peut préserver sujets, essais, matrices, unités, temps moléculaires et provenance. Omics2PK doit encore disposer des mêmes sujets avec doses, concentrations, chronologie partagée et covariables de référence avant toute comparaison pharmacocinétique.',
        'An export can preserve subjects, assays, matrices, units, molecular timing and provenance. Omics2PK still requires the same subjects with doses, concentrations, a shared time origin and reference covariates before any pharmacokinetic comparison.'
      )}</p>
    </div>
    <div class="flow" aria-label={t('Flux de recherche','Research flow')}>
      <span>{t('Données omiques','Omics data')}</span><b>→</b>
      <span>MultiOmics</span><b>→</b>
      <span>{t('Échange individuel','Individual exchange')}</span><b>→</b>
      <span>Omics2PK</span><b>→</b>
      <span>{t('Validation PK','PK validation')}</span>
    </div>
  </section>

  <p class="boundary" role="note">{t(
    'Usage recherche : aucune page de cet écosystème ne certifie automatiquement causalité biologique, validité clinique ou recommandation de dose.',
    'Research use: no page in this ecosystem automatically certifies biological causality, clinical validity or a dosing recommendation.'
  )}</p>
</section>

<style>
  .connections-page{padding-block:clamp(2.5rem,6vw,5rem)}
  .intro{max-width:920px}
  .intro h1{max-width:25ch}
  .lede{max-width:76ch;color:var(--text-secondary);font-size:var(--text-lg)}
  .ecosystem-figure{margin:2.5rem 0 3rem;padding:1rem;border:1px solid var(--border-subtle);border-radius:var(--radius-lg);background:var(--bg-secondary)}
  .ecosystem-figure img{display:block;width:100%;height:auto;border-radius:calc(var(--radius-lg) - 6px)}
  .ecosystem-figure figcaption{text-align:center;margin-top:.75rem;color:var(--text-secondary);font-size:var(--text-sm)}
  .roles{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:1rem}
  .connection{padding:1.4rem;border:1px solid var(--border-subtle);border-radius:var(--radius-lg);background:var(--bg-tertiary)}
  .connection.current{border-top:3px solid var(--accent-pk)}
  .connection p{color:var(--text-secondary)}
  .role{font-family:var(--font-mono);font-size:var(--text-xs);font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--accent-pk)!important}
  .actions{display:flex;flex-wrap:wrap;gap:.6rem}
  .handoff{margin-top:2rem;padding:clamp(1.3rem,3vw,2rem);border:1px solid var(--border-subtle);border-radius:var(--radius-lg);background:var(--bg-secondary)}
  .handoff p{color:var(--text-secondary)}
  .flow{display:grid;grid-template-columns:1fr auto 1fr auto 1fr auto 1fr auto 1fr;gap:.5rem;align-items:center;margin-top:1.2rem;text-align:center}
  .flow span{padding:.75rem .55rem;border:1px solid var(--border-strong);border-radius:var(--radius);background:var(--bg-primary)}
  .flow b{color:var(--accent-pk)}
  .boundary{margin-top:2rem;padding:1rem;border-left:3px solid var(--accent-pk);background:var(--bg-secondary);color:var(--text-secondary)}
  @media(max-width:900px){.roles{grid-template-columns:1fr}.flow{grid-template-columns:1fr}.flow b{transform:rotate(90deg)}}
</style>
