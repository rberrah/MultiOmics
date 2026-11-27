<script>
  import { resolve } from '$app/paths';
  import { language } from '#lib/stores/language.js';
  $: t = $language === 'en' ? (/** @type {string} */ fr, /** @type {string} */ en) => en : (/** @type {string} */ fr, /** @type {string} */ en) => fr;
  const exercises = [
    { level:['Débutant','Beginner'], title:['Unité biologique','Biological unit'], prompt:['Airway comporte quatre lignées cellulaires, mesurées chacune sans puis avec dexaméthasone. Combien d’unités biologiques indépendantes pour la comparaison appariée ?','Airway contains four cell lines, each measured untreated and with dexamethasone. How many independent biological units are available for the paired comparison?'], solution:['Quatre, et non huit. Les huit mesures forment quatre paires ; la lignée est le bloc biologique.','Four, not eight. The eight assays form four pairs; cell line is the biological block.'] },
    { level:['Intermédiaire','Intermediate'], title:['Amplitude d’effet','Effect size'], prompt:['Exemple simulé : un gène a une expression moyenne normalisée de 20 dans le contrôle et 80 sous traitement. Calculez le log2 du rapport des moyennes (traitement / contrôle).','Simulated example: a gene has mean normalized expression of 20 in control and 80 under treatment. Calculate the log2 ratio of means (treatment / control).'], solution:['log2(80/20) = 2. Cela représente un rapport de quatre, pas une preuve de significativité ni une estimation DESeq2 ajustée.','log2(80/20) = 2. This is a four-fold ratio, not evidence of significance or an adjusted DESeq2 estimate.'] },
    { level:['Avancé','Advanced'], title:['Validation multi-omique','Multi-omics validation'], prompt:['Nutrimouse contient 40 animaux avec deux couches mesurées par animal. Pour un test sur 10 animaux, combien restent pour apprendre la normalisation et sélectionner les variables ? Où placer les deux couches d’un même animal ?','Nutrimouse has 40 animals with two measured layers per animal. With 10 held-out animals, how many remain for fitting scaling and feature selection? Where should both layers of one animal go?'], solution:['30 animaux d’apprentissage. Les deux couches de chaque animal restent dans le même ensemble. Les 10 animaux de test ne participent ni à la sélection ni au prétraitement appris.','30 training animals. Both layers from each animal stay in the same split. The 10 test animals do not contribute to selection or learned preprocessing.'] }
  ];
</script>
<svelte:head><title>{t('Exercices — Omics Explain','Practice — Omics Explain')}</title></svelte:head>
<a href={resolve('multiomics/')}>{t('Retour à Omics Explain','Back to Omics Explain')}</a>
<h1>{t('Exercices','Practice')}</h1>
<p>{t('Résolvez chaque problème indépendamment, puis ouvrez sa correction. Les cas réels et simulés sont identifiés.','Solve each problem independently, then reveal its worked answer. Real and simulated cases are identified.')}</p>
{#each exercises as exercise}
  <section class="exercise">
    <p>{t(exercise.level[0],exercise.level[1])}</p>
    <h2>{t(exercise.title[0],exercise.title[1])}</h2>
    <p>{t(exercise.prompt[0],exercise.prompt[1])}</p>
    <details><summary>{t('Voir la correction','Reveal worked answer')}</summary><p>{t(exercise.solution[0],exercise.solution[1])}</p></details>
  </section>
{/each}
<p><a href={resolve('multiomics/learn/')}>{t('Revoir le cours','Review the course')}</a> · <a href={resolve('multiomics/game/')}>{t('Jouer à Mo Game','Play Mo Game')}</a> · <a href={resolve('multiomics/tool/')}>{t('Ouvrir l’outil','Open the tool')}</a></p>
<style>
  .exercise{border-top:1px solid var(--border-strong);padding:1.5rem 0;max-width:800px}
  summary{cursor:pointer;min-height:44px;display:list-item;padding:.6rem 0}
</style>
