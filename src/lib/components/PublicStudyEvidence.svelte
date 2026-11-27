<script>
  import { language } from '#lib/stores/language.js';
  import data from '#lib/multiomics/public-game-evidence/measurements.json';
  /** @type {'airway'|'nutrimouse'} */
  export let study = 'airway';
  $: t = $language === 'en' ? (/** @type {string} */ fr, /** @type {string} */ en) => en : (/** @type {string} */ fr, /** @type {string} */ en) => fr;
  let qIndex = 1;
  let subset = 'all';
  const thresholds = [0.01,0.05,0.10];
  const heldout = new Set(data.heldoutScores.map(row => row[0]));
  const pairs = data.airwayPairs.map(row => ({id:String(row[0]),before:Number(row[1]),after:Number(row[2])}));
  const mice = data.nutrimousePoints.map(([id,x,y])=>({id,x,y,held:heldout.has(id)}));
  $: visible = mice.filter(p=>subset==='all'||(subset==='heldout'?p.held:!p.held));
  const yCount = (/** @type {number} */ value)=>210-180*value/6000;
  const xMouse = (/** @type {number} */ value)=>64+346*(value+1.1)/1.2;
  const yMouse = (/** @type {number} */ value)=>210-180*value/3.5;
  $: discoveryCount = data.qCounts[qIndex];
</script>

<section class="evidence" data-testid={'public-study-'+study}>
  {#if study==='airway'}
    <h3>{t('Données réelles : quatre paires Airway','Real data: four Airway pairs')}</h3>
    <p>{t('CRISPLD2 : comptages RNA-seq bruts, avant normalisation. Les traits relient les deux mesures de la même lignée.','CRISPLD2: raw RNA-seq counts, before normalization. Lines connect measurements from the same cell line.')}</p>
    <figure>
      <svg viewBox="0 0 450 260" role="img" aria-label={t('Comptages bruts CRISPLD2, quatre lignées, avant et après dexaméthasone','Raw CRISPLD2 counts, four cell lines, untreated and dexamethasone')}>
        {#each [0,1500,3000,4500,6000] as v}
          <line class="gridline" x1="62" y1={yCount(v)} x2="412" y2={yCount(v)} />
          <text x="54" y={yCount(v)+4} text-anchor="end">{v}</text>
        {/each}
        <path class="axis" d="M62 30 V210 H412"/>
        {#each pairs as p}
          <line class="paired" x1="140" x2="330" y1={yCount(p.before)} y2={yCount(p.after)} />
          <circle class="dot" cx="140" cy={yCount(p.before)} r="4"/>
          <circle class="dot" cx="330" cy={yCount(p.after)} r="4"/>
        {/each}
        <text x="140" y="234" text-anchor="middle">{t('Sans traitement','Untreated')}</text>
        <text x="330" y="234" text-anchor="middle">Dexamethasone</text>
      </svg>
      <figcaption>{t('Comptages bruts : ne pas interpréter le rapport comme un effet ajusté. DESeq2 apparié : log2FC = 2,630 ; q = 2,05 × 10⁻⁴⁵.','Raw counts: do not interpret the ratio as an adjusted effect. Paired DESeq2: log2FC = 2.630; q = 2.05 × 10⁻⁴⁵.')}</figcaption>
    </figure>
    <label class="control">{t('Seuil de p-value ajustée (BH)','Adjusted p-value threshold (BH)')}
      <select bind:value={qIndex} data-testid="public-q-threshold">
        {#each thresholds as q,i}<option value={i}>{q}</option>{/each}
      </select>
    </label>
    <p data-testid="public-discovery-count"><strong>{discoveryCount.toLocaleString($language==='en'?'en-US':'fr-FR')}</strong> / 17 199 {t('gènes déclarés au seuil choisi','genes called at the selected threshold')}</p>
    <div class="fraction" aria-hidden="true"><span style:width={100*discoveryCount/17199+'%'}></span></div>
    <p>{t('Changer le seuil ne réestime pas le modèle ni les q-values. La FDR n’est pas le nombre exact de faux positifs de cette étude.','Changing the threshold does not refit the model or recompute adjusted p-values. FDR is not the exact number of false positives in this study.')}</p>
    <details><summary>{t('Voir les comptages exacts','View exact counts')}</summary><table><thead><tr><th scope="col">{t('Lignée','Cell line')}</th><th scope="col">{t('Sans traitement','Untreated')}</th><th scope="col">Dex.</th></tr></thead><tbody>{#each pairs as p}<tr><th scope="row">{p.id}</th><td>{p.before}</td><td>{p.after}</td></tr>{/each}</tbody></table></details>
    <p class="source"><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778">GSE52778</a> · <a href="https://bioconductor.org/packages/airway">airway / Himes et al.</a> · {t('Rejeu natif DESeq2 qualifié, pas une validation externe.','Qualified native DESeq2 replay, not external validation.')}</p>
  {:else}
    <h3>{t('Données réelles : deux couches par animal','Real data: two layers per animal')}</h3>
    <p>{t('Un point = une souris. Expression macroarray hépatique CYP3A11 sur son échelle source, et C14.0 en % des acides gras mesurés. C14.0 est la première variable lipidique du fichier, pas une sélection par corrélation.','One point = one mouse. Liver CYP3A11 macroarray expression on its source scale and C14.0 as a percentage of measured fatty acids. C14.0 is the first lipid variable in the file, not selected by correlation.')}</p>
    <label class="control">{t('Animaux affichés','Displayed animals')}
      <select bind:value={subset} data-testid="public-mouse-subset">
        <option value="all">{t('Tous (40)','All (40)')}</option>
        <option value="train">{t('Apprentissage (30)','Training (30)')}</option>
        <option value="heldout">{t('Test interne (10)','Internal holdout (10)')}</option>
      </select>
    </label>
    <figure>
      <svg viewBox="0 0 450 260" role="img" aria-label={t('Expression CYP3A11 et pourcentage C14.0 par souris','CYP3A11 expression and C14.0 percentage by mouse')}>
        {#each [0,1,2,3] as v}<line class="gridline" x1="64" x2="410" y1={yMouse(v)} y2={yMouse(v)}/><text x="56" y={yMouse(v)+4} text-anchor="end">{v}</text>{/each}
        {#each [-1,-.5,0] as v}<text x={xMouse(v)} y="230" text-anchor="middle">{v}</text>{/each}
        <path class="axis" d="M64 30 V210 H410"/>
        <text x="64" y="18">C14.0 (%)</text>
        <text x="237" y="252" text-anchor="middle">CYP3A11</text>
        {#each visible as p}
          {#if p.held}<rect class="held" x={xMouse(p.x)-4} y={yMouse(p.y)-4} width="8" height="8" data-mouse={p.id}/>{:else}<circle class="dot" cx={xMouse(p.x)} cy={yMouse(p.y)} r="3.5" data-mouse={p.id}/>{/if}
        {/each}
      </svg>
      <figcaption><span class="legend-dot"></span> {t('Apprentissage','Training')} · <span class="legend-square"></span> {t('Test interne','Internal holdout')} · <strong data-testid="public-mouse-count">{visible.length}</strong> {t('animaux','animals')}</figcaption>
    </figure>
    <p>{t('Les paramètres de standardisation de CYP3A11 sont calculés sur 30 animaux seulement : moyenne −0,514 ; écart-type 0,261634. Les 10 autres sont projetés sans réestimation.','CYP3A11 standardization parameters are computed on 30 animals only: mean −0.514; SD 0.261634. The other 10 are projected without refitting.')}</p>
    <details><summary>{t('Voir les dix scores gelés','View the ten frozen scores')}</summary><table><thead><tr><th scope="col">{t('Souris','Mouse')}</th><th scope="col">z</th></tr></thead><tbody>{#each data.heldoutScores as [id,z]}<tr><th scope="row">{String(id).padStart(3,'0')}</th><td>{z.toFixed(3)}</td></tr>{/each}</tbody></table></details>
    <p class="source"><a href="https://github.com/mixOmics-org/mixOmics/blob/ef3e760526623d9e91e945e4b50be48d764efc40/data/nutrimouse.rda">mixOmics · Nutrimouse / Martin et al.</a> · {t('Un seul jeu animal ; aucune mesure PK. Pas de preuve causale ni de validation clinique.','One animal dataset; no PK measurements. No causal evidence or clinical validation.')}</p>
  {/if}
</section>
<style>
  .evidence{border:1px solid var(--border-strong);border-radius:8px;padding:clamp(.75rem,2vw,1.25rem);margin:1rem 0;min-width:0;overflow-wrap:anywhere}
  figure{margin:1rem 0;max-width:620px}
  svg{display:block;width:100%;height:auto;color:var(--text-primary)}
  svg text{font-size:12px;fill:currentColor}
  .axis{fill:none;stroke:currentColor;stroke-width:1}
  .gridline{stroke:currentColor;stroke-opacity:.18}
  .paired{stroke:var(--accent-pk);stroke-width:2;stroke-opacity:.8}
  .dot,.legend-dot{fill:var(--accent-pk);background:var(--accent-pk)}
  .held{fill:var(--bg-primary);stroke:var(--text-primary);stroke-width:2}
  .legend-dot,.legend-square{display:inline-block;width:.7em;height:.7em;vertical-align:baseline}
  .legend-dot{border-radius:50%}.legend-square{border:2px solid currentColor}
  .control{display:flex;flex-wrap:wrap;gap:.5rem;align-items:center;margin:.75rem 0}
  select,summary{min-height:44px;font:inherit;cursor:pointer}
  summary{padding:.5rem 0}
  select{max-width:100%;padding:.4rem;background:var(--bg-secondary);color:var(--text-primary);border:1px solid var(--border-strong)}
  .fraction{height:12px;background:var(--bg-secondary);border:1px solid var(--border-strong);border-radius:6px;overflow:hidden}
  .fraction span{display:block;height:100%;background:var(--accent-pk)}
  figcaption,.source{font-size:.875rem}
  table{width:100%;max-width:460px;border-collapse:collapse;font-variant-numeric:tabular-nums}th,td{text-align:left;padding:.35rem;border-bottom:1px solid var(--border-strong)}
</style>
