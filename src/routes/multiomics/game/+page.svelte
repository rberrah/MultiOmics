<script>
  import { resolve } from '$app/paths';
  import { language } from '#lib/stores/language.js';
  $: t = $language === 'en' ? (/** @type {string} */ fr, /** @type {string} */ en) => en : (/** @type {string} */ fr, /** @type {string} */ en) => fr;
  // Each pair is French / English. Cases teach design and interpretation; no numerical study results are fabricated.
  const paths = [
    {name:'Airway', title:['Transcriptomique : du plan au résultat','Transcriptomics: from design to results'], steps:[
      {q:['Des cellules provenant de mêmes donneurs sont mesurées avec et sans dexaméthasone. Quel plan retenir ?','Cells from the same donors are measured with and without dexamethasone. Which design is appropriate?'], a:[['Groupes indépendants','Independent groups'],['Comparaison appariée par donneur','Donor-paired comparison']],correct:1,why:['Le donneur est un bloc : ignorer l’appariement confond la variabilité inter-donneurs et l’effet du traitement.','Donor is a block: ignoring pairing conflates between-donor variation and treatment effects.']},
      {q:['Que faut-il vérifier avant une analyse différentielle RNA-seq ?','What should be checked before RNA-seq differential analysis?'],a:[['Correspondance des échantillons, comptages et plan d’étude','Sample alignment, counts and study design'],['Seulement les gènes les plus abondants','Only the most abundant genes']],correct:0,why:['Les colonnes doivent correspondre aux métadonnées et le modèle doit être adapté aux comptages et au plan.','Columns must match metadata, and the model must suit counts and design.']},
      {q:['Une ACP sépare visuellement les groupes. Cela suffit-il à prouver un effet significatif ?','A PCA visually separates groups. Does this prove a significant effect?'],a:[['Oui','Yes'],['Non, il faut une inférence adaptée','No, appropriate inference is needed']],correct:1,why:['L’ACP décrit la variance ; elle ne teste pas l’effet du traitement.','PCA describes variance; it does not test the treatment effect.']},
      {q:['Des milliers de gènes sont testés. Quel résultat faut-il considérer avec l’amplitude d’effet ?','Thousands of genes are tested. What should be considered alongside effect size?'],a:[['P-value brute uniquement','Raw p-value only'],['P-value ajustée pour tests multiples','Multiple-testing-adjusted p-value']],correct:1,why:['La correction BH vise le contrôle du taux de fausses découvertes sous ses hypothèses ; l’effet et son incertitude restent essentiels.','BH adjustment aims to control false discovery rate under its assumptions; effect size and uncertainty remain essential.']}
    ]},
    {name:'Nutrimouse',title:['Multi-omiques : intégrer et interpréter','Multi-omics: integrate and interpret'],steps:[
      {q:['Des mesures transcriptomiques et lipidiques concernent les mêmes souris. Comment les relier ?','Transcript and lipid measurements concern the same mice. How should they be linked?'],a:[['Par identifiant unique de souris','By unique mouse identifier'],['Par ordre des colonnes sans vérification','By column order without checking']],correct:0,why:['Un appariement erroné crée des associations artificielles entre les couches.','Incorrect matching creates artificial cross-layer associations.']},
      {q:['Les matrices ont des unités et dimensions différentes. Quelle précaution prendre ?','Matrices have different units and dimensions. What precaution is needed?'],a:[['Concaténer sans préparation','Concatenate without preparation'],['Prétraiter chaque couche et évaluer son poids','Preprocess each layer and assess its weight']],correct:1,why:['Une couche à variance ou dimension dominante peut masquer les autres ; le choix de mise à l’échelle doit être explicite.','A high-variance or high-dimensional layer can dominate; scaling choices must be explicit.']},
      {q:['Quel outil correspond à la recherche non supervisée de facteurs latents partagés ?','Which tool addresses unsupervised discovery of shared latent factors?'],a:[['MOFA','MOFA'],['DIABLO supervisé','Supervised DIABLO']],correct:0,why:['MOFA est non supervisé ; DIABLO utilise les classes pour sélectionner des signatures discriminantes.','MOFA is unsupervised; DIABLO uses class labels to select discriminative signatures.']},
      {q:['Une signature DIABLO discrimine les régimes sur les données d’apprentissage. Que vérifier ?','A DIABLO signature discriminates diets in training data. What should be checked?'],a:[['La performance hors apprentissage, sans fuite de données','Out-of-training performance without data leakage'],['Seulement la séparation sur les données d’apprentissage','Only training-set separation']],correct:0,why:['La sélection et le réglage doivent rester dans les plis d’apprentissage. Une association prédictive ne démontre pas une causalité.','Feature selection and tuning must stay within training folds. Predictive association does not establish causality.']}
    ]}
  ];
  let selected=0, step=0, answer=-1;
  let completed=[0,0];
  $: journey=paths[selected];
  $: question=journey.steps[step];
  $: finished=step===journey.steps.length;
  /** @param {number} index */
  function start(index){ selected=index; step=0; answer=-1; completed=[...completed]; completed[index]=0; }
  function next(){ if(answer!==question.correct)return; completed=[...completed]; completed[selected]=step+1; step+=1; answer=-1; }
</script>
<svelte:head><title>Mo Game — MultiOmics</title></svelte:head>
<a class="back" href={resolve('multiomics/')}>{t('Retour à MultiOmics','Back to MultiOmics')}</a>
<h1>Mo Game</h1>
<p class="lede">{t('Deux parcours interactifs pour apprendre à choisir et interpréter une analyse omique.','Two interactive journeys to learn how to choose and interpret an omics analysis.')}</p>
<p>{t('Cas pédagogiques inspirés de jeux publics. Les questions ne constituent pas une réanalyse des données et ne présentent aucun résultat numérique inventé.','Educational cases inspired by public datasets. The questions do not reanalyze those data or present invented numerical results.')}</p>
<div class="journeys">
{#each paths as path,i}
  <button class:active={selected===i} aria-pressed={selected===i} onclick={()=>start(i)}><strong>{path.name}</strong><span>{t(path.title[0],path.title[1])}</span><small>{completed[i]}/{path.steps.length}</small></button>
{/each}
</div>
<section class="mission" aria-labelledby="mission-title">
  <h2 id="mission-title">{journey.name} · {finished?t('Parcours terminé','Journey complete'):t('Étape','Step')+' '+(step+1)+' / '+journey.steps.length}</h2>
  {#if finished}
    <p>{t('Vous avez validé les quatre décisions. Poursuivez avec les méthodes et les exemples du logiciel.','You completed all four decisions. Continue with the methods and software examples.')}</p>
    <a href={resolve('multiomics/learn/')}>{t('Comprendre les méthodes','Learn the methods')}</a> · <a href={resolve('multiomics/demo/')}>{t('Voir la démonstration','View the demo')}</a> · <a href={resolve('multiomics/tool/')}>{t('Ouvrir l’outil','Open the tool')}</a>
  {:else}
    <h3>{t(question.q[0],question.q[1])}</h3>
    <div class="answers">
      {#each question.a as item,i}<button aria-pressed={answer===i} onclick={()=>answer=i}>{t(item[0],item[1])}</button>{/each}
    </div>
    {#if answer>=0}
      <p role="status"><strong>{answer===question.correct?t('Correct.','Correct.'):t('À revoir.','Try again.')}</strong> {t(question.why[0],question.why[1])}</p>
    {/if}
    <button class="next" disabled={answer!==question.correct} onclick={next}>{t('Continuer','Continue')}</button>
  {/if}
</section>
<style>
.back{display:inline-flex;align-items:center;min-height:44px}
h1{font-size:clamp(2rem,6vw,4rem);margin:.5rem 0}.lede{font-size:1.2rem}
.journeys{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,260px),1fr));gap:1rem;margin:2rem 0}
button{font:inherit;color:var(--text-primary);background:var(--bg-secondary);border:1px solid var(--border-strong);border-radius:8px;min-height:44px;padding:1rem;text-align:left;cursor:pointer}
button:focus-visible,a:focus-visible{outline:2px solid var(--accent-pk);outline-offset:3px}
button[aria-pressed=true]{border:2px solid var(--accent-pk);box-shadow:inset 0 0 0 1px var(--accent-pk)}
.journeys button{display:grid;gap:.5rem}.mission{border-top:1px solid var(--border-strong);padding:1.5rem 0 3rem;max-width:900px}
.answers{display:grid;gap:.75rem}.next{margin-top:1rem}.next:disabled{opacity:.55;cursor:not-allowed}
</style>
