<script>
  import {language} from '#lib/stores/language.js';
  import {batchExample} from '#lib/multiomics/game-batch.js';
  $: t=$language==='en'?(/** @type {string} */ fr,/** @type {string} */ en)=>en:(/** @type {string} */ fr,/** @type {string} */ en)=>fr;
  let crossed=false,batchEffect=2,biologicalEffect=0;
  $: result=batchExample(crossed,batchEffect,biologicalEffect);
  $: fmt=new Intl.NumberFormat($language==='en'?'en':'fr',{maximumFractionDigits:2});
  const y=(/** @type {number} */ value)=>206-(value+.5)*27;
</script>
<section class="batch-lab" data-testid="game-batch-lab">
  <h3>{t('Simulation : le lot peut imiter un effet biologique','Simulation: batch can mimic biological effect')}</h3>
  <p>{t('Mesures fictives d’une variable, deux groupes et deux lots. La dispersion est fixée ; aucun résultat de PCA réelle.','Fictional measurements of one feature, two groups and two batches. Spread is fixed; this is not a real PCA result.')}</p>
  <label><input type="checkbox" bind:checked={crossed} data-testid="cross-batches"/>{t('Représenter les deux groupes dans chaque lot','Represent both groups in each batch')}</label>
  <label>{t('Effet du lot','Batch effect')} : {fmt.format(batchEffect)}<input type="range" min="0" max="4" step=".25" bind:value={batchEffect}/></label>
  <label>{t('Effet biologique simulé','Simulated biological effect')} : {fmt.format(biologicalEffect)}<input type="range" min="0" max="2" step=".25" bind:value={biologicalEffect}/></label>
  <svg viewBox="0 0 430 250" role="img" aria-label={t('Mesures par groupe et lot','Measurements by group and batch')}>
    {#each [0,2,4,6] as tick}<line x1="42" x2="412" y1={y(tick)} y2={y(tick)} stroke="currentColor" opacity=".18"/><text x="35" y={y(tick)+4} text-anchor="end">{tick}</text>{/each}
    {#each result.points as p,i}
      {@const cx=125+p.group*180+(p.batch===0?-14:14)+(i%4-1.5)*5}
      {#if p.batch===0}<circle cx={cx} cy={y(p.value)} r="5" class="batch-a"/>{:else}<rect x={cx-4} y={y(p.value)-4} width="8" height="8" class="batch-b"/>{/if}
    {/each}
    <text x="125" y="235" text-anchor="middle">{t('Groupe A','Group A')}</text><text x="305" y="235" text-anchor="middle">{t('Groupe B','Group B')}</text>
  </svg>
  <p>{t('Cercle : lot 1. Carré : lot 2.','Circle: batch 1. Square: batch 2.')}</p>
  <p>{t('Différence brute B − A','Raw difference B − A')}: <strong data-testid="batch-difference">{fmt.format(result.difference)}</strong></p>
  <p data-testid="batch-identifiability">{result.within?t('Comparaisons intra-lot possibles : ','Within-batch comparisons possible: ')+result.within.map(x=>fmt.format(x)).join(' / '):t('Effets non séparables dans ce plan.','Effects are not separable in this design.')}</p>
</section>
<style>
.batch-lab{padding:1rem;border:1px solid var(--border-strong);border-radius:8px;display:grid;gap:.65rem}.batch-lab label{display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;min-height:44px}.batch-lab input[type="range"]{flex:1;min-width:100px;min-height:24px}svg{display:block;width:100%;max-width:620px;height:auto;color:var(--text-primary)}svg text{font-size:12px;fill:currentColor}.batch-a{fill:var(--accent-pk)}.batch-b{fill:var(--bg-primary);stroke:var(--text-primary);stroke-width:2}
</style>
