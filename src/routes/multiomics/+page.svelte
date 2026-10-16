<script>
  import { base } from '$app/paths';
  import { language } from '$lib/stores/language';
  /** @param {string} fr @param {string} en */
  const t = (fr, en) => $language === 'en' ? en : fr;

  const branches = [
    {
      fr: 'Explorer les données',
      en: 'Explore the data',
      methodFr: 'ACP d’une omique ou ACP multi-blocs',
      methodEn: 'Single-layer PCA or multi-block PCA',
      detailFr: 'Une seule omique : l’ACP explore les variations de cette couche. Plusieurs omiques : l’ACP multi-blocs équilibrée explore les structures communes. Dans les deux cas, il ne s’agit pas d’un test de différences entre groupes.',
      detailEn: 'One omic: PCA explores variation within that layer. Multiple omics: balanced multi-block PCA explores shared structure. Neither is a statistical test of group differences.'
    },
    {
      fr: 'Comparer des groupes',
      en: 'Compare groups',
      methodFr: 'Comparaison adaptée au plan d’étude',
      methodEn: 'Study-design-aware comparison',
      detailFr: 'L’outil tient compte des répétitions techniques, des séries expérimentales, des facteurs d’ajustement et du caractère indépendant, apparié ou longitudinal des données.',
      detailEn: 'The tool accounts for technical replicates, experimental batches, adjustment factors and whether data are independent, paired or longitudinal.'
    },
    {
      fr: 'Relier les données à un critère clinique',
      en: 'Link data to an endpoint',
      methodFr: 'Association avec un critère',
      methodEn: 'Association with an endpoint',
      detailFr: 'Selon le type de critère, l’outil choisit automatiquement le modèle statistique adapté : continu, oui/non, comptage, plusieurs catégories ou survie.',
      detailEn: 'Depending on the endpoint type, the tool automatically selects the appropriate statistical model: continuous, yes/no, count, multiclass or survival.'
    },
    {
      fr: 'Étudier l’évolution dans le temps',
      en: 'Study changes over time',
      methodFr: 'Suivi des mesures chez un même sujet',
      methodEn: 'Repeated-measure follow-up',
      detailFr: 'Les mesures répétées d’une même personne, d’un même animal ou d’une même culture sont reliées afin d’étudier leur évolution dans le temps.',
      detailEn: 'Repeated measurements from the same person, animal or culture are linked to study their evolution over time.'
    }
  ];

  const glossary = [
    ['Metadata', 'Tableau des échantillons', 'Sample sheet'],
    ['Outcome', 'Critère étudié / résultat clinique', 'Endpoint / outcome'],
    ['Batch', 'Série technique ou expérimentale', 'Technical or experimental batch'],
    ['Covariable', 'Facteur à prendre en compte dans l’analyse', 'Adjustment variable'],
    ['Feature', 'Variable biologique mesurée : gène, protéine ou métabolite', 'Measured biological feature'],
    ['q BH / FDR', 'q BH : p-value ajustée par Benjamini–Hochberg ; FDR : proportion attendue de fausses découvertes parmi les découvertes déclarées, dans les hypothèses du modèle', 'BH-adjusted p-value; FDR: expected proportion of false discoveries among declared findings under model assumptions']
  ];
</script>

<svelte:head>
  <title>{t('MultiOmics — Omiques et pharmacométrie PK/PD', 'MultiOmics — Omics for PK/PD research')}</title>
  <meta name="description" content={t('MultiOmics : analyses traçables des données transcriptomiques, protéomiques et métabolomiques, pour explorer des hypothèses en PK/PD.', 'MultiOmics: reproducible transcriptomics, proteomics and metabolomics analyses to explore PK/PD hypotheses.')} />
</svelte:head>

<section class="hero">
  <p class="eyebrow">{t('Omiques → hypothèses PK/PD · outil de recherche', 'Omics → PK/PD hypotheses · research tool')}</p>
  <h1>{t('Analyser une ou plusieurs omiques avec des règles explicites.', 'Analyze one or multiple omics with explicit rules.')}</h1>
  <p class="lede">{t(
    'Choisissez RNA, protéines, métabolites ou leur combinaison. Le moteur confronte vos données à votre plan d’étude et à des bases scientifiques publiques, selon des règles auditables, sans LLM. Une seule omique suffit pour explorer, comparer ou prédire lorsque le plan le permet.',
    'Choose RNA, proteins, metabolites, or a combination. The engine follows auditable rules informed by your study design and public scientific databases, without an LLM. One omics layer is enough to explore, compare or predict when the design supports it.'
  )}</p>
  <div class="actions">
    <a class="btn btn-primary" href={`${base}/multiomics/tool`}>{t('Ouvrir l’outil', 'Open the tool')}</a>
    <a class="btn btn-outline" href={`${base}/multiomics/README.txt`} target="_blank">{t('Voir le format des fichiers', 'See file formats')}</a>
    <a class="btn btn-outline" href={`${base}/multiomics/learn`}>{t('Cours : méthodes et résultats', 'Learn methods and interpretation')}</a>
  </div>
</section>

<section class="panel sample-sheet-intro">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Avant de commencer', 'Before you start')}</p>
      <h2>{t('Le tableau des échantillons relie vos fichiers entre eux', 'The sample sheet links your files together')}</h2>
    </div>
    <p>{t('Ce tableau est parfois appelé “metadata”. Il ne contient pas les valeurs RNA, protéines ou métabolites : il indique à quoi correspond chaque colonne des matrices. Pour explorer une seule omique, un assistant peut créer un tableau minimal, uniquement après confirmation que chaque colonne représente un sujet indépendant.', 'This is sometimes called “metadata”. It describes what matrix columns represent. For one-omic exploratory PCA, an assistant can create a minimal sheet only after you confirm that every column is one independent subject.')}</p>
  </div>

  <div class="data-kind-grid">
    <article class="sample-sheet-card">
      <span class="kind-tag">{t('Tableau des échantillons', 'Sample sheet')}</span>
      <h3>{t('Il décrit le contexte', 'It describes the context')}</h3>
      <p>{t('Qui ? Quel prélèvement ? Quelle mesure ? Quel groupe ? Quel moment ? Quelle série technique ? Il ne contient pas les milliers de valeurs moléculaires.', 'Who? Which specimen? Which assay? Which group? Which time point? Which technical batch? It does not contain the thousands of molecular values.')}</p>
      <code>subject_id · sample_id · assay_id · omic · condition · timepoint</code>
    </article>
    <article>
      <span class="kind-tag">{t('Matrices de mesures', 'Measurement matrices')}</span>
      <h3>{t('Elles contiennent les valeurs biologiques', 'They contain biological measurements')}</h3>
      <p>{t('Une matrice RNA, une matrice protéomique et/ou une matrice métabolomique. Les lignes correspondent aux gènes, protéines ou métabolites ; les colonnes correspondent aux mesures de vos échantillons.', 'An RNA, proteomics and/or metabolomics matrix. Rows correspond to genes, proteins or metabolites; columns correspond to measurements from your samples.')}</p>
      <code>feature_id · RNA001 · RNA002 · …</code>
    </article>
  </div>

  <div class="link-example">
    <div class="example-title">
      <strong>{t('Exemple : un même prélèvement mesuré dans trois omiques', 'Example: one specimen measured in three omics layers')}</strong>
      <span>{t('assay_id indique quelle colonne de matrice correspond à chaque mesure.', 'assay_id indicates which matrix column corresponds to each assay.')}</span>
    </div>
    <div class="example-table" role="table" aria-label={t('Exemple de tableau des échantillons', 'Example sample sheet')}>
      <div class="example-head"><b>subject_id</b><b>sample_id</b><b>assay_id</b><b>omic</b><b>condition</b><b>timepoint</b></div>
      <div><code>P001</code><code>P001_T0</code><code>RNA001</code><span>transcriptomics</span><span>control</span><span>T0</span></div>
      <div><code>P001</code><code>P001_T0</code><code>PROT001</code><span>proteomics</span><span>control</span><span>T0</span></div>
      <div><code>P001</code><code>P001_T0</code><code>MET001</code><span>metabolomics</span><span>control</span><span>T0</span></div>
    </div>
    <p class="note">{t('Lecture : P001 est la même personne ; P001_T0 est le même prélèvement ; RNA001, PROT001 et MET001 sont trois mesures différentes de ce prélèvement. Chaque assay_id doit correspondre exactement à une colonne de la matrice concernée.', 'Read it as follows: P001 is the same person; P001_T0 is the same specimen; RNA001, PROT001 and MET001 are three different assays of that specimen. Each assay_id must match exactly one column in the relevant matrix.')}</p>
  </div>
</section>

<section class="panel" aria-label={t('Bien interpréter les résultats', 'Responsible interpretation')}>
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Avant toute conclusion', 'Before drawing conclusions')}</p>
      <h2>{t('Un résultat calculé n’est pas automatiquement une découverte validée', 'A computed result is not automatically a validated discovery')}</h2>
    </div>
    <p>{t(
      'Le navigateur propose des calculs exploratoires et le moteur R fournit des méthodes de référence quand elles sont adaptées. Une q-value faible ne prouve pas de mécanisme biologique ni un contrôle universel des faux positifs. Vérifiez l’unité biologique indépendante, l’appariement, les confusions, la qualité des données et les limites de chaque méthode.',
      'The browser provides exploratory calculations, while the R engine uses native reference methods where appropriate. A low adjusted p-value does not prove a biological mechanism or universally controlled false discoveries. Verify independent biological units, pairing, confounding, quality control and method-specific limitations.'
    )}</p>
  </div>
</section>

<section class="panel glossary-panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Vocabulaire', 'Vocabulary')}</p>
      <h2>{t('Les termes techniques sont traduits en langage courant', 'Technical terms are translated into plain language')}</h2>
    </div>
    <p>{t('Dans l’outil, le terme simple doit être affiché en premier. Le terme statistique reste disponible dans l’aide “?” ou entre parenthèses lorsqu’il est nécessaire pour la reproductibilité.', 'In the tool, the plain-language term should be shown first. The statistical term remains available in “?” help or in parentheses when needed for reproducibility.')}</p>
  </div>
  <div class="glossary-grid">
    {#each glossary as item}
      <div>
        <code>{item[0]}</code>
        <strong>{$language === 'en' ? item[2] : item[1]}</strong>
      </div>
    {/each}
  </div>
</section>

<section class="panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Principe', 'Principle')}</p>
      <h2>{t('Votre question scientifique choisit l’analyse', 'Your scientific question selects the analysis')}</h2>
    </div>
    <p>{t('Vous décrivez l’étude et les données. Le moteur choisit ensuite une branche d’analyse selon des règles fixes et explicables.', 'You describe the study and data. The engine then selects an analysis branch using fixed, explainable rules.')}</p>
  </div>

  <div class="branch-grid">
    {#each branches as branch}
      <article>
        <span>{$language === 'en' ? branch.methodEn : branch.methodFr}</span>
        <h3>{$language === 'en' ? branch.en : branch.fr}</h3>
        <p>{$language === 'en' ? branch.detailEn : branch.detailFr}</p>
      </article>
    {/each}
  </div>
</section>

<section class="panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Comment lire le tableau des échantillons', 'How to read the sample sheet')}</p>
      <h2>{t('Quatre colonnes expliquent l’essentiel', 'Four columns explain the essentials')}</h2>
    </div>
    <p>{t('Pensez simplement : qui a été étudié, quel prélèvement a été réalisé, quelle mesure correspond à la colonne de matrice, et de quel type d’omique il s’agit.', 'Think simply: who was studied, which specimen was collected, which assay matches the matrix column, and which omics layer it belongs to.')}</p>
  </div>
  <div class="flow">
    <article><b>subject_id</b><p>{t('Qui ? La personne, l’animal ou l’unité biologique.', 'Who? The person, animal or biological unit.')}</p></article>
    <article><b>sample_id</b><p>{t('Quel prélèvement ? Par exemple un tube ou un fragment de tissu à une visite donnée.', 'Which specimen? For example a tube or tissue sample at a given visit.')}</p></article>
    <article><b>assay_id</b><p>{t('Quelle mesure ? L’identifiant exact de la colonne correspondante dans la matrice.', 'Which assay? The exact identifier of the corresponding matrix column.')}</p></article>
    <article><b>omic</b><p>{t('Quel type de mesure ? RNA, protéomique ou métabolomique.', 'Which measurement type? RNA, proteomics or metabolomics.')}</p></article>
  </div>
</section>

<section class="panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Facteurs à prendre en compte', 'Factors to account for')}</p>
      <h2>{t('Les différences techniques et cliniques peuvent être intégrées à l’analyse', 'Technical and clinical differences can be included in the analysis')}</h2>
    </div>
  </div>
  <div class="two-col">
    <article>
      <h3>{t('Série technique (batch)', 'Technical batch')}</h3>
      <p>{t('Exemple : des échantillons analysés sur des plaques, jours ou séries différentes. Si une série technique correspond exactement à un seul groupe, l’outil signale que leurs effets ne peuvent pas être séparés. Sinon, cet effet technique est pris en compte dans le modèle.', 'Example: samples analyzed on different plates, days or runs. If a technical batch exactly matches one biological group, the tool reports that their effects cannot be separated. Otherwise, the technical effect is included in the model.')}</p>
    </article>
    <article>
      <h3>{t('Facteurs d’ajustement (covariables)', 'Adjustment factors (covariates)')}</h3>
      <p>{t('Exemples : âge, sexe, centre, traitement concomitant. Vous choisissez les colonnes qui doivent être prises en compte ; l’outil ne décide pas seul de ce qui constitue un facteur de confusion.', 'Examples: age, sex, centre, concomitant treatment. You choose which columns should be accounted for; the tool does not decide on its own what constitutes a confounder.')}</p>
    </article>
  </div>
</section>

<section class="panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Interprétation', 'Interpretation')}</p>
      <h2>{t('Chaque résultat doit répondre à une question simple', 'Each result should answer a simple question')}</h2>
    </div>
    <p>{t('Les détails statistiques restent disponibles, mais l’interface doit d’abord expliquer ce que le nombre signifie biologiquement.', 'Statistical details remain available, but the interface should first explain what the number means biologically.')}</p>
  </div>
  <div class="interpret">
    <div><b>{t('Amplitude de l’effet', 'Effect size')}</b><span>{t('de combien la variable change et dans quel sens', 'how much the feature changes and in which direction')}</span></div>
    <div><b>{t('Résultat corrigé (q BH)', 'Adjusted result (BH q)')}</b><span>{t('tient compte du grand nombre de gènes, protéines ou métabolites testés', 'accounts for the large number of genes, proteins or metabolites tested')}</span></div>
    <div><b>{t('Lien entre omiques', 'Cross-omics link')}</b><span>{t('uniquement si au moins deux omiques sont mesurées et compatibles', 'only when at least two compatible omics layers are measured')}</span></div>
    <div><b>Reactome</b><span>{t('replace les résultats dans des voies biologiques connues ; cette information vient d’une base externe', 'places results in known biological pathways; this information comes from an external database')}</span></div>
  </div>
</section>

<section class="panel validation">
  <div>
    <p class="eyebrow">{t('Validation publique', 'Public validation')}</p>
    <h2>{t('Le pipeline est testé sur plusieurs jeux de données publics', 'The pipeline is tested on several public datasets')}</h2>
  </div>
  <p>{t('La suite automatisée couvre notamment Nutrimouse, TCGA breast, IntLIM NCI-60/BRCA, AgingHFCD, STATegra, LRRK2, PaintOmics à signal planté et missRows.', 'The automated suite covers Nutrimouse, TCGA breast, IntLIM NCI-60/BRCA, AgingHFCD, STATegra, LRRK2, a planted-signal PaintOmics dataset and missRows.')}</p>
  <a class="btn btn-primary" href={`${base}/multiomics/tool`}>{t('Lancer une analyse', 'Run an analysis')}</a>
</section>

<style>
  .hero { max-width: 930px; padding: var(--space-12) 0 var(--space-10); }
  h1 { font-size: clamp(2.6rem, 6vw, 5rem); line-height: .98; max-width: 15ch; margin: var(--space-3) 0 var(--space-6); letter-spacing: -.045em; }
  h2 { margin: 0; font-size: var(--text-2xl); }
  h3 { margin: 0 0 var(--space-2); }
  p { line-height: 1.65; }
  .lede { max-width: 72ch; color: var(--text-secondary); font-size: var(--text-lg); }
  .actions { display:flex; flex-wrap:wrap; gap:var(--space-3); margin-top:var(--space-6); }
  .panel { border-top:1px solid var(--border-strong); padding-top:var(--space-6); margin-top:var(--space-10); }
  .section-head { display:flex; justify-content:space-between; gap:var(--space-6); align-items:end; margin-bottom:var(--space-5); }
  .section-head > p { color:var(--text-secondary); max-width:48ch; margin:0; }
  .sample-sheet-intro { border-top:3px solid var(--accent-pk); }
  .data-kind-grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:var(--space-4); }
  .data-kind-grid article { padding:var(--space-5); border:1px solid var(--border-subtle); border-radius:var(--radius); background:var(--bg-secondary); }
  .sample-sheet-card { border-color:var(--accent-pk) !important; }
  .kind-tag { display:inline-block; margin-bottom:var(--space-3); font-family:var(--font-mono); font-size:var(--text-xs); color:var(--accent-pk); }
  .data-kind-grid code { display:block; margin-top:var(--space-3); color:var(--text-secondary); white-space:normal; }
  .link-example { margin-top:var(--space-5); padding:var(--space-5); border:1px solid var(--border-strong); border-radius:var(--radius); overflow-x:auto; }
  .example-title { display:flex; justify-content:space-between; gap:var(--space-4); align-items:baseline; margin-bottom:var(--space-3); }
  .example-title span { color:var(--text-secondary); font-size:var(--text-sm); }
  .example-table { min-width:720px; border:1px solid var(--border-subtle); }
  .example-table > div { display:grid; grid-template-columns:.8fr 1fr .9fr 1.2fr .9fr .7fr; }
  .example-table > div > * { padding:8px 10px; border-right:1px solid var(--border-subtle); border-bottom:1px solid var(--border-subtle); overflow:hidden; text-overflow:ellipsis; }
  .example-head { background:var(--bg-secondary); }
  .glossary-panel { border-top:3px solid var(--accent-pd); }
  .glossary-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:var(--space-3); }
  .glossary-grid > div { display:grid; gap:5px; padding:var(--space-4); border:1px solid var(--border-subtle); background:var(--bg-secondary); }
  .glossary-grid code { color:var(--accent-pk); font-weight:700; }
  .branch-grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:var(--space-4); }
  .branch-grid article, .two-col article { padding:var(--space-5); border:1px solid var(--border-subtle); border-radius:var(--radius); background:var(--bg-secondary); }
  .branch-grid article > span { font-family:var(--font-mono); font-size:var(--text-xs); color:var(--accent-pk); }
  .branch-grid p, .two-col p, .note { color:var(--text-secondary); }
  .flow { display:grid; grid-template-columns:repeat(4,1fr); gap:var(--space-3); }
  .flow article { padding:var(--space-4); border-top:2px solid var(--accent-pd); background:var(--bg-secondary); }
  .flow b { font-family:var(--font-mono); color:var(--accent-pk); }
  .two-col { display:grid; grid-template-columns:repeat(2,1fr); gap:var(--space-4); }
  .interpret { display:grid; grid-template-columns:repeat(2,1fr); gap:var(--space-3); }
  .interpret div { display:grid; grid-template-columns:.9fr 2fr; gap:var(--space-3); padding:var(--space-4); border:1px solid var(--border-subtle); }
  .interpret span { color:var(--text-secondary); }
  .validation { margin-bottom:var(--space-12); }
  .validation > p { max-width:70ch; color:var(--text-secondary); }
  @media (max-width:900px) {
    .glossary-grid { grid-template-columns:repeat(2,minmax(0,1fr)); }
  }
  @media (max-width:700px) {
    .section-head, .example-title { align-items:start; flex-direction:column; }
    .branch-grid, .flow, .two-col, .interpret, .data-kind-grid, .glossary-grid { grid-template-columns:1fr; }
  }
</style>