<script>
  import {language} from '#lib/stores/language.js';
  import data from '#lib/multiomics/public-game-evidence/airway-volcano.json';
  $: t=$language==='en'?(/** @type {string} */ fr,/** @type {string} */ en)=>en:(/** @type {string} */ fr,/** @type {string} */ en)=>fr;
  let qIndex=1, effectIndex=0;
  $: threshold=data.q[qIndex];
  $: effect=data.fx[effectIndex];
  $: count=data.counts[qIndex][effectIndex];
  $: fmt=new Intl.NumberFormat($language==='en'?'en':'fr',{maximumFractionDigits:2});
  const x=(/** @type {number} */ value)=>48+(value+6)*26;
  const y=(/** @type {number} */ value)=>248-value*18;
  const max=Math.max(...data.bins.flat());
  const density=(/** @type {number} */ n)=>n===0?0:.15+.85*Math.log1p(n)/Math.log1p(max);
</script>
<section class="volcano" data-testid="airway-volcano">
  <h3>{t('Lire un volcano plot : résultats réels DESeq2','Read a volcano plot: real DESeq2 results')}</h3>
  <p>{t('Les 17 199 gènes du contraste apparié sont représentés par densité. L’axe horizontal indique le log2FC ajusté ; l’axe vertical −log10(q BH). Les cases plus foncées contiennent davantage de gènes (échelle logarithmique).','All 17,199 genes in the paired contrast are shown as density. The horizontal axis is fitted log2FC; the vertical axis is −log10(BH-adjusted p-value). Darker cells contain more genes (logarithmic scale).')}</p>
  <div class="controls">
    <label>{t('Seuil q BH','BH-adjusted p-value threshold')}<select data-testid="volcano-q" bind:value={qIndex}>{#each data.q as value,i}<option value={i}>{value}</option>{/each}</select></label>
    <label>{t('Amplitude minimale |log2FC|','Minimum absolute log2FC')}<select data-testid="volcano-effect" bind:value={effectIndex}>{#each data.fx as value,i}<option value={i}>{value}</option>{/each}</select></label>
  </div>
  <figure>
    <svg viewBox="0 0 490 296" role="img" aria-label={t('Densité des effets et q-values des 17 199 gènes, avec seuils de sélection','Density of effects and adjusted p-values for 17,199 genes, with selection thresholds')}>
      {#each data.bins as row,j}
        {#each row as n,i}
          {#if n>0}<rect x={x(i-6)} y={y(j+1)} width="26" height="18" class="density" fill-opacity={density(n)}><title>{n} {t('gènes','genes')}</title></rect>{/if}
        {/each}
      {/each}
      {#each [-6,-4,-2,0,2,4,6,8,10] as tick}<text x={x(tick)} y="264" text-anchor="middle">{tick}</text>{/each}
      {#each [0,2,4,6,8,10] as tick}<text x="42" y={y(tick)+4} text-anchor="end">{tick}</text>{/each}
      <text x="42" y={y(11.5)+4} text-anchor="end">≥11</text>
      <path d="M48 32 V248 H464" fill="none" stroke="currentColor"/>
      <line x1="48" x2="464" y1={y(-Math.log10(threshold))} y2={y(-Math.log10(threshold))} class="limit"/>
      {#if effect>0}{#each [-effect,effect] as v}<line x1={x(v)} x2={x(v)} y1="32" y2="248" class="limit"/>{/each}{/if}
      <circle cx={x(data.marker.log2FC)} cy={y(11.5)} r="5" class="marker"/>
      <text x={x(data.marker.log2FC)+9} y={y(11.5)+4}>CRISPLD2</text>
      <text x="256" y="286" text-anchor="middle">log2FC</text>
      <text x="12" y="145" transform="rotate(-90 12 145)" text-anchor="middle">−log10(q BH)</text>
    </svg>
    <figcaption>{t('Agrégation en cases d’une unité ; la dernière rangée regroupe −log10(q) ≥ 11. CRISPLD2 est repéré dans cette rangée, sa valeur exacte étant q = 2,05 × 10⁻⁴⁵. Les lignes en pointillés sont les seuils choisis.','One-unit bins; the top row pools −log10(q) ≥ 11. CRISPLD2 is marked in that row, with exact q = 2.05 × 10⁻⁴⁵. Dashed lines are the selected thresholds.')}</figcaption>
  </figure>
  <p><strong data-testid="volcano-selected">{fmt.format(count)}</strong> / {fmt.format(data.n)} {t('gènes répondent aux deux critères (comptage exact, non estimé à partir des cases).','genes meet both criteria (exact count, not estimated from the bins).')}</p>
  <p>{t('Changer les seuils ne réestime pas le modèle. Le filtre d’amplitude aide à prioriser, mais ne constitue pas un test de seuil d’effet ni une garantie de FDR pour le sous-ensemble filtré. Une association d’expression ne démontre ni mécanisme ni utilité clinique.','Changing thresholds does not refit the model. The effect-size filter helps prioritize but is not an effect-threshold hypothesis test or an FDR guarantee for the filtered subset. Expression association does not establish mechanism or clinical utility.')}</p>
  <p><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778">GSE52778</a> · {t('Rejeu natif apparié qualifié ; pas de nouvelle analyse navigateur.','Qualified paired native replay; no new browser fit.')}</p>
</section>
<style>
.volcano{margin:1rem 0;padding:1rem;border:1px solid var(--border-strong);border-radius:8px;display:grid;gap:.7rem;min-width:0}.controls{display:flex;flex-wrap:wrap;gap:1rem}.controls label{display:grid;gap:.3rem;min-width:0}.controls select{min-height:44px;background:var(--bg-secondary);color:var(--text-primary);border:1px solid var(--border-strong)}figure{margin:0;min-width:0}svg{display:block;width:100%;max-width:680px;height:auto;color:var(--text-primary)}svg text{font-size:10px;fill:currentColor}.density{fill:var(--accent-pk);stroke:var(--bg-primary);stroke-width:.4}.limit{stroke:currentColor;stroke-width:1.5;stroke-dasharray:4 3}.marker{fill:var(--bg-primary);stroke:var(--text-primary);stroke-width:2}figcaption{font-size:.875rem}
</style>
