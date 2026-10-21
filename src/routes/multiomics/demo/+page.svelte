<script>
  import { resolve } from '$app/paths';
  import { language } from '#lib/stores/language.js';
  /** @param {string} fr @param {string} en */
  const t=(fr,en)=>$language==='en'?en:fr;
  const launchHref = resolve('multiomics/tool/') + '?demo=multi';
</script>

<svelte:head>
  <title>{t('Démo guidée — MultiOmics','Guided demo — MultiOmics')}</title>
  <meta name="description" content={t(
    'Cas pratique guidé utilisant le jeu synthétique MultiOmics : plan longitudinal, QC, lecture mono-omique, intégration et interprétation.',
    'Guided practical case using the MultiOmics synthetic dataset: longitudinal design, QC, single-omic reading, integration and interpretation.'
  )} />
</svelte:head>

<section class="guided">
  <header class="hero">
    <p class="eyebrow">MultiOmics · {t('Démo guidée','Guided demo')}</p>
    <h1>{t('Cas pratique : comprendre une réponse au traitement dans le temps','Practical case: understand a treatment response over time')}</h1>
    <p class="lede">{t(
      'Cette page reprend exactement la démonstration multi-omique disponible dans l’outil et l’explique étape par étape. Le but n’est pas d’obtenir une conclusion biologique réelle, mais d’apprendre à lire correctement un workflow multi-omique.',
      'This page uses the same multi-omics demonstration available in the tool and explains it step by step. The goal is not to produce a real biological conclusion, but to learn how to read a multi-omics workflow correctly.'
    )}</p>
    <div class="actions">
      <a class="primary" href={launchHref} data-testid="guided-demo-launch">{t('Lancer cette démo dans l’outil','Run this demo in the tool')}</a>
      <a class="secondary" href={resolve('multiomics/learn/')}>{t('Voir la méthodologie générale','Read the general methodology')}</a>
    </div>
  </header>

  <section class="case" data-testid="guided-demo-design">
    <div>
      <span class="step">01</span>
      <p class="eyebrow">{t('Commencer par le plan','Start with the design')}</p>
      <h2>{t('Quelle est la question ?','What is the question?')}</h2>
    </div>
    <div class="facts">
      <article><strong>8</strong><span>{t('sujets','subjects')}</span></article>
      <article><strong>4 + 4</strong><span>{t('contrôle / traitement','control / treatment')}</span></article>
      <article><strong>T0 · T12</strong><span>{t('deux temps par sujet','two time points per subject')}</span></article>
      <article><strong>3</strong><span>{t('omiques','omics layers')}</span></article>
    </div>
    <p>{t(
      'La question simulée est : le traitement modifie-t-il l’état moléculaire entre T0 et T12, et ces changements sont-ils cohérents entre transcriptome, protéome et métabolome ? Le même sujet est mesuré deux fois : les observations longitudinales ne sont donc pas indépendantes.',
      'The simulated question is: does treatment alter the molecular state between T0 and T12, and are those changes coherent across transcriptomics, proteomics and metabolomics? Each subject is measured twice, so longitudinal observations are not independent.'
    )}</p>
    <p class="interpretation"><strong>{t('À retenir :','Take-home:')}</strong> {t(
      'l’unité biologique reste le sujet. Les 49 assays du fichier de démonstration ne représentent pas 49 sujets indépendants.',
      'the biological unit remains the subject. The 49 assays in the demo file do not represent 49 independent subjects.'
    )}</p>
  </section>

  <section class="case">
    <div>
      <span class="step">02</span>
      <p class="eyebrow">{t('Vérifier avant d’interpréter','Check before interpreting')}</p>
      <h2>{t('Contrôle qualité et structure des données','Quality control and data structure')}</h2>
    </div>
    <div class="two-col">
      <div>
        <h3>{t('Ce que contient la démo','What the demo contains')}</h3>
        <ul>
          <li>{t('16 occasions biologiques : 8 sujets × 2 temps.', '16 biological occasions: 8 subjects × 2 time points.')}</li>
          <li>{t('Un même prélèvement alimente les trois couches omiques.', 'The same specimen is represented across all three omics layers.')}</li>
          <li>{t('Un réplicat technique RNA supplémentaire existe pour P007 à T12 : R015A / R015B.', 'One extra RNA technical replicate exists for P007 at T12: R015A / R015B.')}</li>
          <li>{t('Les batches sont déclarés dans les métadonnées.', 'Batches are explicitly declared in metadata.')}</li>
        </ul>
      </div>
      <div>
        <h3>{t('Comment l’interpréter','How to interpret it')}</h3>
        <p>{t(
          'Un réplicat technique ne crée pas un nouveau sujet. Avant tout test, il faut donc préserver l’appariement sujet–temps et distinguer répétition technique et répétition biologique.',
          'A technical replicate does not create a new subject. Any analysis must preserve subject–time matching and distinguish technical from biological replication.'
        )}</p>
        <p>{t(
          'Des signaux stables comme GAPDH dans la démo servent surtout de repères pédagogiques ; ils ne constituent pas à eux seuls une preuve de qualité expérimentale.',
          'Stable signals such as GAPDH in the demo are pedagogical anchors; they are not by themselves proof of experimental quality.'
        )}</p>
      </div>
    </div>
  </section>

  <section class="case" data-testid="guided-demo-interpretation">
    <div>
      <span class="step">03</span>
      <p class="eyebrow">{t('Lire chaque couche séparément','Read each layer separately')}</p>
      <h2>{t('Quels signaux sont visibles ?','Which signals are visible?')}</h2>
    </div>
    <div class="omics-grid">
      <article>
        <span class="omic">RNA</span>
        <h3>{t('Transcriptomique','Transcriptomics')}</h3>
        <p>{t(
          'Dans les sujets traités à T12, STAT1, IDO1 et KYNU augmentent nettement, alors que CYP3A5 diminue. GAPDH reste relativement stable.',
          'In treated subjects at T12, STAT1, IDO1 and KYNU increase clearly, while CYP3A5 decreases. GAPDH remains comparatively stable.'
        )}</p>
        <small>{t('Lecture : signal cohérent, mais encore descriptif dans le navigateur.', 'Reading: coherent signal, but still descriptive in the browser.')}</small>
      </article>
      <article>
        <span class="omic">PROT</span>
        <h3>{t('Protéomique','Proteomics')}</h3>
        <p>{t(
          'Les protéines correspondant à STAT1, IDO1 et KYNU suivent la même direction générale à T12 chez les sujets traités, tandis que GAPDH reste stable.',
          'Proteins corresponding to STAT1, IDO1 and KYNU follow the same overall direction at T12 in treated subjects, while GAPDH remains stable.'
        )}</p>
        <small>{t('Lecture : concordance inter-couches, pas preuve de causalité.', 'Reading: cross-layer concordance, not proof of causality.')}</small>
      </article>
      <article>
        <span class="omic">MET</span>
        <h3>{t('Métabolomique','Metabolomics')}</h3>
        <p>{t(
          'Le L-tryptophane diminue tandis que la L-kynurénine augmente chez les sujets traités à T12. Le lactate et le succinate augmentent aussi dans ce jeu synthétique.',
          'L-tryptophan decreases while L-kynurenine increases in treated subjects at T12. Lactate and succinate also increase in this synthetic dataset.'
        )}</p>
        <small>{t('Lecture : la direction est compatible avec le motif simulé, sans validation biologique réelle.', 'Reading: the direction is compatible with the simulated pattern, without real biological validation.')}</small>
      </article>
    </div>
  </section>

  <section class="case">
    <div>
      <span class="step">04</span>
      <p class="eyebrow">{t('Intégrer sans surinterpréter','Integrate without over-interpreting')}</p>
      <h2>{t('Que raconte l’ensemble des omiques ?','What does the combined dataset suggest?')}</h2>
    </div>
    <div class="flow" aria-label={t('Flux interprétatif de la démonstration','Interpretive flow of the demo')}>
      <span>STAT1 ↑</span><b>→</b><span>IDO1 / KYNU ↑</span><b>→</b><span>{t('Tryptophane ↓','Tryptophan ↓')}</span><b>→</b><span>{t('Kynurénine ↑','Kynurenine ↑')}</span>
    </div>
    <p>{t(
      'Dans cette démonstration construite, plusieurs niveaux moléculaires racontent volontairement une histoire cohérente autour de la voie du tryptophane–kynurénine. C’est précisément ce que MultiOmics doit aider à voir : une convergence entre couches, tout en conservant la provenance de chaque signal.',
      'In this constructed demonstration, several molecular levels intentionally tell a coherent story around the tryptophan–kynurenine pathway. This is what MultiOmics is meant to help reveal: convergence across layers while preserving the provenance of each signal.'
    )}</p>
    <p class="interpretation"><strong>{t('À retenir :','Take-home:')}</strong> {t(
      'une convergence transcriptome–protéome–métabolome renforce la cohérence d’une hypothèse, mais ne transforme pas une association en mécanisme causal.',
      'transcriptome–proteome–metabolome convergence strengthens the coherence of a hypothesis, but does not turn an association into a causal mechanism.'
    )}</p>
  </section>

  <section class="case">
    <div>
      <span class="step">05</span>
      <p class="eyebrow">{t('Lire les figures','Read the figures')}</p>
      <h2>{t('Dans quel ordre regarder les sorties ?','In what order should the outputs be read?')}</h2>
    </div>
    <ol class="reading">
      <li><strong>{t('Résumé du plan et QC.','Design summary and QC.')}</strong> {t('Vérifier sujets, temps, couches, manquants et réplicats avant toute histoire biologique.', 'Check subjects, time points, layers, missingness and replicates before any biological narrative.')}</li>
      <li><strong>{t('Heatmap / structure globale.','Heatmap / global structure.')}</strong> {t('Repérer les motifs dominants, sans confondre séparation visuelle et significativité statistique.', 'Identify dominant patterns without confusing visual separation with statistical significance.')}</li>
      <li><strong>{t('Effets par couche.','Layer-specific effects.')}</strong> {t('Examiner direction, amplitude et stabilité des signaux dans chaque omique.', 'Inspect direction, magnitude and stability of signals in each layer.')}</li>
      <li><strong>{t('Réseau et voies.','Network and pathways.')}</strong> {t('Vérifier si les signaux convergent vers une même zone biologique et si les identifiants sont correctement résolus.', 'Check whether signals converge on the same biological region and whether identifiers are correctly resolved.')}</li>
      <li><strong>{t('Intégration multi-omique.','Multi-omics integration.')}</strong> {t('Chercher une cohérence entre couches plutôt qu’un unique score spectaculaire.', 'Look for cross-layer coherence rather than one spectacular score.')}</li>
    </ol>
  </section>

  <section class="case final">
    <div>
      <span class="step">06</span>
      <p class="eyebrow">{t('Passer de la démo à une étude réelle','Move from the demo to a real study')}</p>
      <h2>{t('Exploration ici, inférence de référence ensuite','Explore here, then run reference inference')}</h2>
    </div>
    <p>{t(
      'La démonstration navigateur sert à apprendre le workflow et à visualiser les tendances. Pour une étude réelle destinée à un article, il faut exporter le dossier reproductible, vérifier le plan expérimental et exécuter le code R de référence adapté aux données.',
      'The browser demonstration is for learning the workflow and visualizing trends. For a real study intended for publication, export the reproducible dossier, verify the study design and run the appropriate reference R code.'
    )}</p>
    <div class="actions">
      <a class="primary" href={launchHref}>{t('Rejouer la démo','Replay the demo')}</a>
      <a class="secondary" href={resolve('multiomics/learn/')}>{t('Comprendre les méthodes en détail','Understand the methods in detail')}</a>
    </div>
  </section>

  <aside class="limit" data-testid="guided-demo-limit">
    <strong>{t('Limite essentielle','Essential limitation')}</strong>
    <p>{t(
      'Toutes les valeurs de ce cas sont synthétiques et volontairement construites pour produire un motif lisible. Elles ne valident ni un biomarqueur, ni une voie biologique, ni une décision clinique.',
      'All values in this case are synthetic and deliberately constructed to create an interpretable pattern. They validate neither a biomarker, nor a biological pathway, nor a clinical decision.'
    )}</p>
  </aside>
</section>

<style>
  .guided{max-width:1080px;margin:0 auto;padding:clamp(1.5rem,4vw,3rem) 0 5rem}
  .hero{max-width:860px;margin-bottom:3rem}
  .eyebrow{font-size:.8rem;font-weight:750;letter-spacing:.09em;text-transform:uppercase;color:var(--accent-pk)}
  h1{font-size:clamp(2rem,5vw,3.35rem);line-height:1.08;margin:.4rem 0 1rem}
  .lede{font-size:var(--text-lg);line-height:1.7;color:var(--text-secondary);max-width:72ch}
  .actions{display:flex;flex-wrap:wrap;gap:.75rem;margin-top:1.4rem}
  .actions a{display:inline-block;padding:.72rem 1rem;border-radius:7px;text-decoration:none;font-weight:700}
  .primary{background:var(--accent-pk);color:var(--bg-primary)}
  .secondary{border:1px solid var(--border-strong);color:var(--text-primary);background:var(--bg-secondary)}
  .case{margin-top:1.2rem;padding:clamp(1.2rem,3vw,2rem);border:1px solid var(--border-subtle);border-radius:14px;background:var(--bg-primary)}
  .case>div:first-child{display:grid;grid-template-columns:auto 1fr;column-gap:.8rem;align-items:start}
  .case>div:first-child .eyebrow,.case>div:first-child h2{grid-column:2}
  .step{grid-row:1 / span 2;font-family:var(--font-mono);font-weight:800;color:var(--accent-pk);padding-top:.08rem}
  h2{font-size:clamp(1.35rem,3vw,1.8rem);margin:.15rem 0 1rem}
  h3{font-size:1rem;margin-top:0}
  .facts{display:grid!important;grid-template-columns:repeat(4,minmax(0,1fr))!important;gap:.7rem;margin:1rem 0}
  .facts article{padding:1rem;background:var(--bg-secondary);border-radius:10px;border:1px solid var(--border-subtle)}
  .facts strong,.facts span{display:block}.facts strong{font-size:1.35rem}.facts span{font-size:.85rem;color:var(--text-secondary)}
  .two-col{display:grid!important;grid-template-columns:1fr 1fr!important;gap:1rem;margin-top:1rem}
  .two-col>div{padding:1rem;background:var(--bg-secondary);border-radius:10px}
  .two-col li{margin-bottom:.45rem}
  .interpretation{padding:.9rem 1rem;border-left:3px solid var(--accent-pk);background:var(--bg-secondary);color:var(--text-secondary)}
  .omics-grid{display:grid!important;grid-template-columns:repeat(3,minmax(0,1fr))!important;gap:1rem;margin-top:1rem}
  .omics-grid article{padding:1rem;border:1px solid var(--border-subtle);border-radius:10px;background:var(--bg-secondary)}
  .omics-grid small{color:var(--text-muted);line-height:1.5}
  .omic{display:inline-block;font-family:var(--font-mono);font-size:.72rem;font-weight:800;color:var(--accent-pk);margin-bottom:.6rem}
  .flow{display:flex!important;grid-template-columns:none!important;flex-wrap:wrap;align-items:center;gap:.55rem;margin:1rem 0 1.2rem}
  .flow span{padding:.65rem .8rem;border:1px solid var(--border-strong);border-radius:8px;background:var(--bg-secondary);font-weight:650}
  .flow b{color:var(--accent-pk)}
  .reading{padding-left:1.2rem}.reading li{margin:.8rem 0;line-height:1.6}
  .final{border-color:var(--border-strong);background:var(--bg-secondary)}
  .limit{margin-top:1.2rem;padding:1rem 1.2rem;border-left:4px solid var(--accent-ai);background:var(--bg-secondary)}
  .limit p{margin-bottom:0;color:var(--text-secondary)}
  @media(max-width:760px){
    .facts,.omics-grid,.two-col{grid-template-columns:1fr 1fr!important}
  }
  @media(max-width:520px){
    .facts,.omics-grid,.two-col{grid-template-columns:1fr!important}
    .actions{align-items:stretch}.actions a{width:100%;text-align:center}
  }
</style>
