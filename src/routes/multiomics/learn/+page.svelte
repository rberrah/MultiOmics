<script>
  import { base } from '$app/paths';
  import { language } from '$lib/stores/language';
  /** @param {string} fr @param {string} en */
  const t = (fr,en) => $language === 'en' ? en : fr;
  const lessons = [
    {
      title:'PCA / ACP',
      fr:'Explore les différences entre profils. Ne teste pas des effets différentiels et ne prouve pas une séparation biologique.',
      en:'Explores profile structure. Does not test differential effects or prove a biological separation.'
    },
    {
      title:'DESeq2 / edgeR-QL',
      fr:'Pour des comptages RNA-seq bruts entiers : modèles adaptés aux comptes et à la dispersion, avec contrastes, lots et bloc sujet si apparié.',
      en:'For integer raw RNA-seq counts: count-appropriate models and dispersion, with declared contrasts, batches and subject blocking if paired.'
    },
    {
      title:'limma / voom',
      fr:'Pour des données continues correctement transformées, ou RNA-seq avec voom. Ne pas appliquer directement limma à des comptages bruts.',
      en:'For appropriately transformed continuous data, or RNA-seq with voom. Do not use ordinary limma directly on raw counts.'
    },
    {
      title:'fgsea',
      fr:'Recherche des enrichissements de voies biologiques sur un classement de toutes les variables éligibles, en figeant voies et identifiants.',
      en:'Tests pathway enrichment on a signed ranking of all eligible features, with frozen pathways and identifier mappings.'
    },
    {
      title:'MOFA2',
      fr:'Modélisation latente non supervisée de plusieurs omiques ; vérifier appariement des sujets, normalisation et contribution de chaque omique.',
      en:'Unsupervised multi-omic latent modeling; verify subject matching, normalization and contribution of each omics block.'
    },
    {
      title:'DIABLO',
      fr:'Classification supervisée multi-omique avec validation sur sujets distincts. Prétraitement et sélection uniquement dans les données d’entraînement.',
      en:'Supervised multi-omic classification with independent-subject testing. Fit preprocessing and feature selection only on training data.'
    }
  ];
</script>

<svelte:head>
  <title>{t('Cours : méthodes et résultats omiques','Omics lessons: methods and interpretation')}</title>
  <meta name="description" content={t('Cours débutant sur les statistiques, les méthodes et l’export R Markdown de MultiOmics.','Beginner guide to omics statistics, methods and reproducible R Markdown.')} />
</svelte:head>

<div class="learn">
  <header class="intro">
    <p class="eyebrow">MultiOmics · {t('Apprendre','Learn')}</p>
    <h1>{t('Choisir une méthode et comprendre ses résultats','Choose a method and interpret its results')}</h1>
    <p>{t('Des fiches courtes pour comprendre la démarche, sans programmation et sans multiplier les réglages.', 'Short lessons explaining each decision without coding or unnecessary settings.')}</p>
    <a class="primary" href={base + '/multiomics/tool'}>{t('Ouvrir MultiOmics','Open MultiOmics')}</a>
  </header>

  <section>
    <h2>{t('1. Choisir selon sa question','1. Choose based on your question')}</h2>
    <div class="grid">
      <article><h3>{t('Explorer','Explore')}</h3><p>{t('ACP mono-omique ou exploration multi-omique, puis MOFA2 si les sujets sont réellement appariés.', 'Single-omic PCA or multi-omic exploration, then MOFA2 when subjects truly match.')}</p></article>
      <article><h3>{t('Comparer','Compare')}</h3><p>{t('Comptages RNA-seq : DESeq2/edgeR. Valeurs continues : limma selon les prétraitements.', 'RNA-seq counts: DESeq2/edgeR. Continuous values: limma with suitable preprocessing.')}</p></article>
      <article><h3>{t('Interpréter','Interpret')}</h3><p>{t('fgsea pour les voies, avec un classement signé complet et une base de voies versionnée.', 'fgsea for pathways, with a complete signed ranking and versioned pathway definitions.')}</p></article>
      <article><h3>{t('Prédire','Predict')}</h3><p>{t('DIABLO avec entraînement et validation séparés, sans fuite de données.', 'DIABLO with separate training and validation, without data leakage.')}</p></article>
    </div>
  </section>

  <section>
    <h2>{t('2. Les méthodes expliquées','2. Methods explained')}</h2>
    <p>{t('Cliquez sur une méthode pour lire sa condition d’utilisation.', 'Expand a method for its main usage conditions.')}</p>
    {#each lessons as lesson (lesson.title)}
      <details class="lesson">
        <summary>{lesson.title}</summary>
        <p>{$language === 'en' ? lesson.en : lesson.fr}</p>
      </details>
    {/each}
  </section>

  <section>
    <h2>{t('3. Comment lire les figures et statistiques','3. How to interpret figures and statistics')}</h2>
    <div class="grid">
      <article><h3>log2FC</h3><p>{t('Un effet +1 signifie un rapport 2 ; −1 signifie 0,5, si la mesure est bien une différence log2. Vérifier le sens du contraste.', 'A log2 effect +1 means a ratio of 2; −1 means 0.5 when the reported statistic is a log2 difference. Check contrast direction.')}</p></article>
      <article><h3>p-value / q-value</h3><p>{t('Une p-value ne donne pas la probabilité que l’hypothèse soit vraie. La q-value correspond à une correction multiple ; aucune ne prouve la causalité.', 'A p-value is not the probability the hypothesis is true. A q-value incorporates multiple-testing adjustment; neither proves causality.')}</p></article>
      <article><h3>Heatmap</h3><p>{t('Un z-score par variable montre sa variation relative. Une couleur forte ne représente pas nécessairement une quantité absolue élevée.', 'Per-feature z-scores describe relative variation. A strong color does not necessarily mean a high absolute abundance.')}</p></article>
      <article><h3>{t('Manquants et zéros','Missing values and zero')}</h3><p>{t('Un zéro mesuré diffère d’une valeur manquante. En MS, la non-détection et le signal sous la limite de quantification demandent des vérifications.', 'A measured zero differs from a missing value. In MS, non-detection and below-quantification signals require explicit handling.')}</p></article>
    </div>
    <p class="caution">{t('Les calculs exploratoires JavaScript peuvent utiliser des modèles simplifiés et ne sont pas toujours égaux aux résultats R. Pour une inférence statistique, préférer les méthodes natives R adaptées au plan et au type de données. Les résultats restent soumis aux hypothèses et au contrôle qualité.', 'Exploratory JavaScript calculations can use simplified models and need not equal native R results. For formal inference, use a design-appropriate native R method and evaluate assumptions and QC.')}</p>
  </section>

  <section>
    <h2>{t('4. Le code R commenté et reproductible','4. Annotated reproducible R code')}</h2>
    <p>{t('Après l’analyse, téléchargez le ZIP. Il comprend ANALYSE_REPRODUCTIBLE.Rmd, avec explications du code, des méthodes, des tableaux et des figures, ainsi que le code R d’origine.', 'Download the research ZIP after analysis. It includes ANALYSE_REPRODUCTIBLE.Rmd with explanations of methods, source code, tables and figures, plus the original R implementation.')}</p>
    <pre><code>Rscript render_notebook.R</code></pre>
    <p>{t('Les dépendances de rapport manquantes peuvent être installées automatiquement. En revanche, les versions des packages scientifiques ne sont jamais changées silencieusement : cela empêcherait de garantir les mêmes résultats.', 'Missing report packages can be installed automatically. Scientific package versions are never silently changed, since that would compromise numerical reproducibility.')}</p>
    <p>{t('Le rapport HTML est généré sur votre ordinateur. Attention : les fichiers peuvent contenir des données confidentielles.', 'The HTML report is generated locally. Files may contain confidential information.')}</p>
  </section>

  <section>
    <h2>{t('5. Et vers la PK/PD ?','5. Linking omics to PK/PD')}</h2>
    <p>{t('Le premier travail consiste à associer avec certitude les sujets, prélèvements, tissus et temps aux données d’exposition, puis à tester sans fuite de données les nouvelles covariables PK/PD. Cette extension n’est pas encore validée.', 'First reliably join patients, specimens, tissues and collection times to drug exposure data, then evaluate new PK/PD covariates without data leakage. That integration is not yet validated.')}</p>
  </section>
</div>

<style>
  .learn{max-width:1060px;margin:auto;padding:2rem 1rem 4rem;line-height:1.6}
  .intro,.lesson,.grid article{border:1px solid var(--border-subtle,#d1d9e5);border-radius:.75rem;background:var(--bg-secondary,#f7f9fc);padding:1rem}
  .intro{padding:1.5rem}.eyebrow{font-weight:700;letter-spacing:.08em;text-transform:uppercase;font-size:.85rem}
  h1{font-size:clamp(1.8rem,4vw,2.5rem);line-height:1.2}h2{font-size:1.45rem;margin-bottom:.8rem}h3{font-size:1.05rem}
  section{margin-top:2.5rem}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,230px),1fr));gap:1rem}
  .lesson{margin:.55rem 0}.lesson summary{cursor:pointer;font-weight:700}
  .primary{display:inline-block;background:#234f9e;color:white;border-radius:.5rem;padding:.7rem 1rem;text-decoration:none}
  .primary:focus-visible,.lesson summary:focus-visible{outline:3px solid #2563eb;outline-offset:3px}
  .caution{border-left:4px solid #bb740c;padding:1rem;background:var(--bg-secondary,#f8fafc)}
  pre{overflow-x:auto;background:#16243a;color:white;border-radius:.5rem;padding:1rem}
</style>
