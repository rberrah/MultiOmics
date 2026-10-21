<script>
  import { onMount } from 'svelte';
  import { resolve, asset } from '$app/paths';
  import { language } from '#lib/stores/language.js';
  import { runDeterministicAnalysis, parseDelimited } from '#lib/multiomics/deterministic.js';
  import { buildMultiomicsVisualizationData } from '#lib/multiomics/visualization-data.js';
  import { assessScientificAssurance } from '#lib/multiomics/scientific-assurance.js';
  import MultiomicsVisualizations from '#lib/components/MultiomicsVisualizations.svelte';
  import {
    MULTIOMICS_DEMO_COLUMN_MAPPING,
    MULTIOMICS_DEMO_DATA_TYPES,
    MULTIOMICS_DEMO_IDENTIFIER_TYPES,
    MULTIOMICS_DEMO_PROTOCOL,
    MULTIOMICS_DEMO_UI_STATE
  } from '#lib/multiomics/demo-case.js';

  /** @param {string} fr @param {string} en */
  const t=(fr,en)=>$language==='en'?en:fr;
  const launchHref = resolve('multiomics/tool/') + '?demo=multi';
  const demoAssets = {
    metadata: asset('multiomics/demo_metadata.csv'),
    transcriptomics: asset('multiomics/demo_transcriptomics.csv'),
    proteomics: asset('multiomics/demo_proteomics.csv'),
    metabolomics: asset('multiomics/demo_metabolomics.csv')
  };
  const layerOrder=['transcriptomics','proteomics','metabolomics'];

  let status='loading';
  let error='';
  /** @type {any} */
  let analysisResult=null;
  /** @type {Record<string,string>[]} */
  let metadataRows=[];
  /** @type {Record<string,{features:number,assays:number}>} */
  let inputSummary={};

  /** @param {string} layer */
  function layerLabel(layer){
    if(layer==='transcriptomics') return t('Transcriptomique','Transcriptomics');
    if(layer==='proteomics') return t('Protéomique','Proteomics');
    return t('Métabolomique','Metabolomics');
  }
  /** @param {unknown} value */
  function sig(value){
    return Number.isFinite(Number(value)) ? Number(value).toPrecision(3) : '—';
  }
  /** @param {unknown} value */
  function effect(value){
    return Number.isFinite(Number(value)) ? Number(value).toPrecision(4) : '—';
  }
  /** @param {any} row @param {any} layer */
  function rowReading(row,layer){
    const qualified=Number.isFinite(row?.qValue) && row.qValue<=0.10;
    const direction=Number(row?.effect)>0
      ? t('changement plus positif dans le groupe traitement','more positive change in the treatment group')
      : Number(row?.effect)<0
        ? t('changement plus négatif dans le groupe traitement','more negative change in the treatment group')
        : t('pas de différence de changement détectée','no detected difference in change');
    const multiplicity=qualified
      ? t('La variable franchit le seuil exploratoire q BH ≤ 0,10 utilisé par l’outil.','The feature passes the exploratory BH q ≤ 0.10 display threshold used by the tool.')
      : t('La variable ne franchit pas le seuil exploratoire q BH ≤ 0,10.','The feature does not pass the exploratory BH q ≤ 0.10 threshold.');
    const ratio=layer?.effectScale==='log2' && Number.isFinite(row?.foldRatio)
      ? t(` Ici, 2^effet = ${row.foldRatio.toFixed(2)}× : c’est le ratio des changements T12/T0 entre traitement et contrôle sur l’échelle log2.`,
          ` Here, 2^effect = ${row.foldRatio.toFixed(2)}×: this is the ratio of T12/T0 changes between treatment and control on the log2 scale.`)
      : '';
    return direction+'. '+multiplicity+ratio;
  }

  const questionRows=()=>[
    [t('Cadre de l’étude','Study setting'),t('Intervention clinique','Clinical intervention'),MULTIOMICS_DEMO_UI_STATE.studySetting],
    [t('Unité biologique','Biological unit'),t('Participant / sujet','Participant / subject'),MULTIOMICS_DEMO_UI_STATE.unitType],
    [t('Plan','Design'),t('Mesures répétées','Repeated measures'),MULTIOMICS_DEMO_UI_STATE.designType],
    [t('Objectif principal','Main objective'),t('Évolution dans le temps','Change over time'),MULTIOMICS_DEMO_UI_STATE.objective],
    [t('Groupes','Groups'),'2',MULTIOMICS_DEMO_UI_STATE.groupCount],
    [t('Temps','Time points'),'2 (T0, T12)',MULTIOMICS_DEMO_UI_STATE.timepointCount],
    [t('Même prélèvement entre omiques ?','Same specimen across omics?'),t('Oui','Yes'),MULTIOMICS_DEMO_UI_STATE.sampleOverlap],
    [t('Réplicats techniques attendus ?','Technical replicates expected?'),t('Oui','Yes'),MULTIOMICS_DEMO_UI_STATE.technicalReplicatesExpected],
    [t('Batch connu ?','Known batch?'),t('Oui','Yes'),MULTIOMICS_DEMO_UI_STATE.batchKnown],
    [t('Omiques partielles attendues ?','Partial omics expected?'),t('Non','No'),MULTIOMICS_DEMO_UI_STATE.partialOmicsExpected],
    [t('Transcriptome','Transcriptome'),t('Comptages bruts · symboles de gènes','Raw counts · gene symbols'),MULTIOMICS_DEMO_UI_STATE.transcriptomicsValues],
    [t('Protéome','Proteome'),t('Intensités log · UniProt','Log intensities · UniProt'),MULTIOMICS_DEMO_UI_STATE.proteomicsValues],
    [t('Métabolome','Metabolome'),t('Aires de pics · ChEBI','Peak areas · ChEBI'),MULTIOMICS_DEMO_UI_STATE.metabolomicsValues],
    [t('Niveau d’analyse','Analysis intent'),t('Exploratoire navigateur','Browser exploratory'),MULTIOMICS_DEMO_UI_STATE.analysisIntent],
    [t('Services externes','External services'),t('Désactivés dans la démo','Disabled in the demo'),'Reactome / mapping = off']
  ];

  async function loadExactDemo(){
    status='loading';
    error='';
    try{
      /** @type {Record<string,File>} */
      const files={};
      for(const layer of ['metadata',...layerOrder]){
        const response=await fetch(demoAssets[layer]);
        if(!response.ok) throw new Error('Demo file unavailable: '+layer);
        const text=await response.text();
        files[layer]=new File([text],`demo_${layer}.csv`,{type:'text/csv'});
        const parsed=parseDelimited(text);
        if(layer==='metadata') metadataRows=parsed.rows;
        else inputSummary={...inputSummary,[layer]:{features:parsed.rows.length,assays:Math.max(0,parsed.headers.length-1)}};
      }
      const result=await runDeterministicAnalysis({
        files,
        metadataRows,
        columnMapping:{...MULTIOMICS_DEMO_COLUMN_MAPPING},
        protocol:{...MULTIOMICS_DEMO_PROTOCOL},
        dataTypes:{...MULTIOMICS_DEMO_DATA_TYPES},
        identifierTypes:{...MULTIOMICS_DEMO_IDENTIFIER_TYPES},
        resolveIdentifiers:MULTIOMICS_DEMO_UI_STATE.resolveIdentifiers,
        useReactome:MULTIOMICS_DEMO_UI_STATE.useReactome
      });
      result.visualizations=await buildMultiomicsVisualizationData({
        files,
        metadataRows,
        columnMapping:{...MULTIOMICS_DEMO_COLUMN_MAPPING},
        dataTypes:{...MULTIOMICS_DEMO_DATA_TYPES},
        analysisResult:result,
        metaboliteAnnotations:[]
      });
      result.scientificAssurance=assessScientificAssurance(result,{demo:true});
      analysisResult=result;
      status='done';
    }catch(e){
      status='error';
      error=e instanceof Error?e.message:String(e);
    }
  }

  onMount(()=>{ void loadExactDemo(); });
</script>

<svelte:head>
  <title>{t('Cas pratique MultiOmics — sorties exactes','MultiOmics worked example — exact outputs')}</title>
  <meta name="description" content={t(
    'Cas pratique MultiOmics calculé avec les mêmes fichiers, paramètres, moteur statistique et visualisations que l’outil.',
    'MultiOmics worked example calculated with the same files, settings, statistical engine and visualizations as the analysis tool.'
  )} />
</svelte:head>

<section class="hero">
  <p class="eyebrow">{t('Cas pratique','Worked example')}</p>
  <h1>{t('MultiOmics — du questionnaire aux sorties exactes','MultiOmics — from questionnaire to exact outputs')}</h1>
  <p class="lede">{t(
    'Comme le cas pratique PMx, on suit ici un projet complet. La différence importante : les tableaux et figures ci-dessous ne sont pas recopiés à la main. Cette page relance directement les mêmes données synthétiques, la même configuration et le même moteur que le bouton « Multi-omiques » de l’outil.',
    'As in the PMx worked example, this page follows a complete project. The important difference is that the tables and figures below are not manually copied. This page directly reruns the same synthetic data, configuration and engine as the “Multi-omics” button in the analysis tool.'
  )}</p>
  <p class="warn">{t('Données synthétiques · exploration navigateur · aucune conclusion clinique ou biomarqueur validé.','Synthetic data · browser exploration · no clinical conclusion or validated biomarker claim.')}</p>
  <p class="steps-map">{t(
    'Données → réponses au questionnaire → méthode sélectionnée → sorties exactes → figures → interprétation → analyse R de référence',
    'Data → questionnaire answers → selected method → exact outputs → figures → interpretation → reference R analysis'
  )}</p>
  <div class="actions">
    <a class="primary" href={launchHref} data-testid="guided-demo-launch">{t('Ouvrir exactement ce cas dans Tool','Open this exact case in Tool')}</a>
    <a class="secondary" href={resolve('multiomics/learn/')}>{t('Méthodologie générale','General methodology')}</a>
  </div>
</section>

<section data-testid="guided-demo-design">
  <h2>{t('1 · Le jeu de données','1 · Dataset')}</h2>
  <p>{t(
    'Le cas simule une intervention avec 8 sujets, 4 contrôles et 4 traités. Chaque sujet est mesuré à T0 et T12 dans les trois couches omiques. Un réplicat technique RNA supplémentaire est présent pour P007 à T12.',
    'The case simulates an intervention with 8 subjects, 4 controls and 4 treated subjects. Each subject is measured at T0 and T12 in all three omics layers. One extra RNA technical replicate is present for P007 at T12.'
  )}</p>

  {#if status==='done'}
    <div class="facts" data-testid="guided-demo-live-counts">
      <article><strong>{analysisResult.metadataSummary.subjects}</strong><span>{t('sujets','subjects')}</span></article>
      <article><strong>{analysisResult.metadataSummary.samples}</strong><span>{t('prélèvements biologiques','biological samples')}</span></article>
      <article><strong>{analysisResult.metadataSummary.assays}</strong><span>assays</span></article>
      <article><strong>{Object.keys(analysisResult.layers).length}</strong><span>{t('couches analysées','analysed layers')}</span></article>
    </div>
    <div class="table-wrap">
      <table>
        <thead><tr><th>{t('Fichier','File')}</th><th>{t('Variables','Features')}</th><th>Assays</th><th>{t('Type déclaré','Declared type')}</th></tr></thead>
        <tbody>
          <tr><td>demo_transcriptomics.csv</td><td>{inputSummary.transcriptomics?.features}</td><td>{inputSummary.transcriptomics?.assays}</td><td>raw_counts</td></tr>
          <tr><td>demo_proteomics.csv</td><td>{inputSummary.proteomics?.features}</td><td>{inputSummary.proteomics?.assays}</td><td>log_intensity</td></tr>
          <tr><td>demo_metabolomics.csv</td><td>{inputSummary.metabolomics?.features}</td><td>{inputSummary.metabolomics?.assays}</td><td>peak_area</td></tr>
        </tbody>
      </table>
    </div>
  {:else if status==='loading'}
    <p class="note">{t('Calcul de la démonstration avec le moteur MultiOmics…','Running the demonstration through the MultiOmics engine…')}</p>
  {:else}
    <p class="error">{error}</p>
  {/if}

  <details>
    <summary>{t('Voir quelques lignes du tableau des échantillons','Show sample-sheet rows')}</summary>
    <div class="table-wrap compact">
      <table>
        <thead><tr><th>subject_id</th><th>sample_id</th><th>omic</th><th>condition</th><th>timepoint</th><th>batch</th><th>technical_replicate</th></tr></thead>
        <tbody>
          {#each metadataRows.slice(0,12) as row}
            <tr><td>{row.subject_id}</td><td>{row.sample_id}</td><td>{row.omic}</td><td>{row.condition}</td><td>{row.timepoint}</td><td>{row.batch}</td><td>{row.technical_replicate}</td></tr>
          {/each}
        </tbody>
      </table>
    </div>
  </details>
</section>

<section>
  <h2>{t('2 · Les réponses que l’utilisateur doit fournir','2 · Answers the user would provide')}</h2>
  <p>{t(
    'Le bouton de démonstration remplit automatiquement ces réponses. Dans une vraie étude, ce sont elles qui déterminent la branche statistique ; MultiOmics ne doit pas deviner le plan expérimental à partir des valeurs moléculaires.',
    'The demo button fills these answers automatically. In a real study, these answers determine the statistical branch; MultiOmics should not infer the experimental design from molecular values.'
  )}</p>
  <div class="table-wrap">
    <table data-testid="guided-demo-questionnaire">
      <thead><tr><th>{t('Question','Question')}</th><th>{t('Réponse affichée','Displayed answer')}</th><th>{t('Valeur utilisée par le moteur','Engine value')}</th></tr></thead>
      <tbody>
        {#each questionRows() as row}
          <tr><td>{row[0]}</td><td><strong>{row[1]}</strong></td><td><code>{row[2]}</code></td></tr>
        {/each}
      </tbody>
    </table>
  </div>
</section>

<section>
  <h2>{t('3 · Quel test l’outil choisit-il ?','3 · Which test does the tool choose?')}</h2>
  {#if status==='done'}
    <p>{t(
      'Avec deux temps et des mesures répétées, l’outil calcule d’abord, sujet par sujet, le changement Δ = T12 − T0. Il compare ensuite ces changements entre contrôle et traitement. Comme chaque couche contient moins de 500 variables, la branche navigateur utilise ici le test de permutation déterministe prévu par le moteur, puis ajuste les p-values par Benjamini–Hochberg dans chaque couche.',
      'With two time points and repeated measures, the tool first computes each subject’s change Δ = T12 − T0. It then compares those changes between control and treatment. Because each layer has fewer than 500 features, the browser branch uses the engine’s deterministic permutation test here, then adjusts p-values with Benjamini–Hochberg within each layer.'
    )}</p>
    <div class="method-grid" data-testid="guided-demo-methods">
      {#each layerOrder as layer}
        {@const r=analysisResult.layers[layer]}
        <article>
          <p class="eyebrow">{layerLabel(layer)}</p>
          <h3>{r.contrast}</h3>
          <dl>
            <div><dt>{t('Mode','Mode')}</dt><dd><code>{r.mode}</code></dd></div>
            <div><dt>{t('Test','Test')}</dt><dd>{r.inferenceMethod}</dd></div>
            <div><dt>n</dt><dd>{r.groupSizes?.[0]} vs {r.groupSizes?.[1]}</dd></div>
            <div><dt>{t('Échelle','Scale')}</dt><dd><code>{r.effectScale}</code></dd></div>
          </dl>
          <p class="note">{r.inferencePolicy}</p>
          <p class="note">{r.selectionRule}</p>
        </article>
      {/each}
    </div>
    <p class="cap">{t(
      'Le réplicat technique RNA est agrégé avant l’inférence ; il ne devient jamais un neuvième sujet. Les batchs déclarés sont audités avant le contraste.',
      'The RNA technical replicate is aggregated before inference; it never becomes a ninth subject. Declared batches are audited before the contrast.'
    )}</p>
  {:else}
    <p class="note">{t('Le choix du test apparaîtra après le calcul.','The selected test will appear after calculation.')}</p>
  {/if}
</section>

<section data-testid="guided-demo-exact-results">
  <h2>{t('4 · Sorties numériques exactes de Tool','4 · Exact numerical outputs from Tool')}</h2>
  <p>{t(
    'Les tableaux ci-dessous lisent directement result.layers du moteur. La précision d’affichage p/q est la même que dans Tool. « Effet » reste affiché en plus pour montrer le calcul sous-jacent ; lorsque l’échelle est log2, Tool présente principalement 2^effet.',
    'The tables below read result.layers directly from the engine. p/q display precision matches Tool. “Effect” is also shown here to expose the underlying calculation; on a log2 scale, Tool mainly displays 2^effect.'
  )}</p>
  {#if status==='done'}
    {#each layerOrder as layer}
      {@const r=analysisResult.layers[layer]}
      <div class="result-block" data-testid={'guided-demo-layer-'+layer}>
        <h3>{layerLabel(layer)} — {r.contrast}</h3>
        <div class="table-wrap">
          <table>
            <thead>
              <tr><th>{t('Variable','Feature')}</th><th>{t('Effet','Effect')}</th><th>{t('Échelle','Scale')}</th><th>2^effect</th><th>p</th><th>q BH</th><th>n ref</th><th>n comp</th></tr>
            </thead>
            <tbody>
              {#each r.rows as row}
                <tr data-feature={row.feature}>
                  <td><code>{row.feature}</code></td>
                  <td>{effect(row.effect)}</td>
                  <td><code>{row.effectScale}</code></td>
                  <td>{row.foldRatio==null?'—':row.foldRatio.toFixed(2)+'×'}</td>
                  <td>{sig(row.pValue)}</td>
                  <td class:qualified={Number.isFinite(row.qValue)&&row.qValue<=0.10}>{sig(row.qValue)}</td>
                  <td>{row.nReference ?? '—'}</td>
                  <td>{row.nComparison ?? '—'}</td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
      </div>
    {/each}
  {/if}
</section>

<section>
  <h2>{t('5 · Comment interpréter ces sorties ?','5 · How should these outputs be interpreted?')}</h2>
  {#if status==='done'}
    <p>{t(
      'Le contraste est toujours « changement T12−T0 chez les traités » moins « changement T12−T0 chez les contrôles ». Un effet positif signifie donc une augmentation relative plus forte sous traitement ; un effet négatif, une diminution relative plus forte. Ce n’est pas simplement une comparaison traitement vs contrôle à T12.',
      'The contrast is always “T12−T0 change in treated subjects” minus “T12−T0 change in controls”. A positive effect therefore means a stronger relative increase under treatment; a negative effect means a stronger relative decrease. This is not simply a treatment-versus-control comparison at T12.'
    )}</p>
    <div class="interpret-grid" data-testid="guided-demo-interpretation">
      {#each layerOrder as layer}
        {@const r=analysisResult.layers[layer]}
        <article>
          <h3>{layerLabel(layer)}</h3>
          {#each r.rows.slice(0,Math.min(4,r.rows.length)) as row}
            <div class="reading-row">
              <code>{row.feature}</code>
              <p><strong>{row.foldRatio==null?effect(row.effect):row.foldRatio.toFixed(2)+'×'}</strong> · p={sig(row.pValue)} · q={sig(row.qValue)}<br>{rowReading(row,r)}</p>
            </div>
          {/each}
        </article>
      {/each}
    </div>
    <p class="cap">{t(
      'q BH est une p-value ajustée pour les tests multiples de la couche. Ici q ≤ 0,10 sert au classement exploratoire. Sur 8 sujets synthétiques, même un petit q ne transforme pas ce cas en validation biologique externe.',
      'BH q is a multiple-testing-adjusted p-value within the layer. Here q ≤ 0.10 is used for exploratory ranking. With 8 synthetic subjects, even a small q does not turn this case into external biological validation.'
    )}</p>
  {/if}
</section>

<section>
  <h2>{t('6 · Les figures exactes rendues par MultiOmics','6 · Exact figures rendered by MultiOmics')}</h2>
  <p>{t(
    'Ce composant est le même MultiomicsVisualizations que celui utilisé par Tool ; les valeurs ne sont donc pas reconstituées dans une seconde implémentation pédagogique.',
    'This is the same MultiomicsVisualizations component used by Tool, so values are not reconstructed in a second educational implementation.'
  )}</p>
  {#if status==='done' && analysisResult.visualizations}
    <div class="viz" data-testid="guided-demo-live-visualizations">
      <MultiomicsVisualizations result={analysisResult} language={$language ?? 'fr'} />
    </div>
  {:else if status==='error'}
    <p class="error">{error}</p>
  {:else}
    <p class="note">{t('Préparation des figures…','Preparing figures…')}</p>
  {/if}
</section>

<section>
  <h2>{t('7 · Intégration inter-omique : que peut-on dire ?','7 · Cross-omics integration: what can be said?')}</h2>
  {#if status==='done'}
    <div class="facts small">
      <article><strong>{analysisResult.crossOmics?.testedPairs ?? 0}</strong><span>{t('paires testées','pairs tested')}</span></article>
      <article><strong>{analysisResult.crossOmics?.significantPairs ?? 0}</strong><span>{t('paires q ≤ 0,10','pairs q ≤ 0.10')}</span></article>
    </div>
    <p class="note"><strong>{t('Méthode publiée par le moteur :','Method reported by the engine:')}</strong> {analysisResult.crossOmics?.method || '—'}</p>
    {#if analysisResult.crossOmics?.pairs?.length}
      <div class="table-wrap">
        <table>
          <thead><tr><th>{t('Paire','Pair')}</th><th>Pattern</th><th>r ref</th><th>r comp</th><th>Δr</th><th>q BH</th></tr></thead>
          <tbody>
            {#each analysisResult.crossOmics.pairs.slice(0,15) as pair}
              <tr>
                <td><code>{pair.featureA}</code> ↔ <code>{pair.featureB}</code></td>
                <td>{pair.pattern?.replaceAll('_',' ')}</td>
                <td>{sig(pair.rReference)}</td><td>{sig(pair.rComparison)}</td><td>{sig(pair.deltaR)}</td><td>{sig(pair.qValue)}</td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    {:else}
      <p class="cap">{t('Aucune paire inter-omique n’est estimable dans cette configuration avec le critère du moteur.','No cross-omics pair is estimable under the engine criterion for this configuration.')}</p>
    {/if}
    <p class="cap">{t(
      'Un changement de corrélation décrit un changement de couplage statistique entre variables. Il ne prouve ni causalité, ni chaîne enzymatique, ni mécanisme.',
      'A correlation change describes a change in statistical coupling between variables. It proves neither causality, an enzymatic chain, nor mechanism.'
    )}</p>
  {/if}
</section>

<section>
  <h2>{t('8 · Ce qui est exploratoire, et ce qu’il faut refaire pour un article','8 · What is exploratory, and what should be rerun for a paper')}</h2>
  {#if status==='done'}
    <p><strong>{t('Statut de cette exécution :','Status of this run:')}</strong> <code>{analysisResult.engine.execution}</code> · <code>{analysisResult.protocol.analysisIntent}</code>.</p>
  {/if}
  <p>{t(
    'Cette démo sert à comprendre le flux et à vérifier visuellement les résultats. Pour une étude réelle, le ZIP de recherche fourni par Tool doit être conservé avec les données, paramètres et provenance, puis le code R de référence doit être exécuté pour l’inférence destinée à une publication. Le navigateur reste utile pour l’exploration ; il n’est pas présenté comme une validation confirmatoire.',
    'This demo is for understanding the workflow and visually checking results. For a real study, the research ZIP provided by Tool should be retained with data, settings and provenance, then the reference R code should be run for publication-oriented inference. The browser remains useful for exploration; it is not presented as confirmatory validation.'
  )}</p>
  <div class="actions">
    <a class="primary" href={launchHref}>{t('Rejouer ce cas dans Tool','Replay this case in Tool')}</a>
    <a class="secondary" href={resolve('multiomics/learn/')}>{t('Comprendre chaque méthode','Understand each method')}</a>
  </div>
</section>

<aside class="limit" data-testid="guided-demo-limit">
  <strong>{t('Point de sécurité scientifique','Scientific safety point')}</strong>
  <p>{t(
    'Les données sont synthétiques et volontairement structurées pour l’enseignement. La page affiche les sorties réellement calculées par MultiOmics, mais l’interprétation biologique reste un exercice pédagogique et non une preuve externe.',
    'The data are synthetic and deliberately structured for teaching. The page displays outputs actually calculated by MultiOmics, but biological interpretation remains a teaching exercise rather than external evidence.'
  )}</p>
</aside>

<style>
  :global(main){max-width:1120px}
  .hero{max-width:900px;margin-bottom:var(--space-10,3rem)}
  .eyebrow{font-family:var(--font-mono);font-size:var(--text-xs,.78rem);text-transform:uppercase;letter-spacing:.08em;color:var(--accent-pk);font-weight:700}
  h1{font-size:clamp(2rem,5vw,3.3rem);line-height:1.08;margin:.45rem 0 1rem}
  .lede{font-size:var(--text-lg,1.1rem);line-height:1.7;max-width:78ch;color:var(--text-secondary)}
  .warn{padding:.8rem 1rem;border-left:3px solid var(--accent-ai);background:var(--bg-secondary);max-width:80ch}
  .steps-map{font-family:var(--font-mono);font-size:.82rem;color:var(--text-secondary);line-height:1.6}
  section{margin:0 0 3.4rem}
  section>h2{font-size:clamp(1.35rem,3vw,1.85rem);margin-bottom:1rem}
  h3{margin:.4rem 0 .7rem}
  .actions{display:flex;flex-wrap:wrap;gap:.7rem;margin-top:1.25rem}
  .actions a{padding:.72rem 1rem;border-radius:6px;text-decoration:none;font-weight:700}
  .primary{background:var(--accent-pk);color:var(--bg-primary)}
  .secondary{border:1px solid var(--border-strong);color:var(--text-primary)}
  .facts{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:.7rem;margin:1.2rem 0}
  .facts.small{grid-template-columns:repeat(2,minmax(0,220px))}
  .facts article{padding:1rem;border:1px solid var(--border-subtle);border-radius:8px;background:var(--bg-secondary)}
  .facts strong,.facts span{display:block}.facts strong{font-size:1.45rem}.facts span{font-size:.84rem;color:var(--text-secondary)}
  .table-wrap{overflow-x:auto;margin:1rem 0}
  table{width:100%;border-collapse:collapse;font-size:.92rem}
  th,td{text-align:left;padding:.68rem .72rem;border-bottom:1px solid var(--border-subtle);vertical-align:top}
  th{font-size:.78rem;text-transform:uppercase;letter-spacing:.035em;color:var(--text-secondary);background:var(--bg-secondary)}
  td code{white-space:nowrap}
  .compact{max-height:340px;overflow:auto}
  details{margin-top:1rem;border:1px solid var(--border-subtle);border-radius:8px;padding:.75rem 1rem}
  summary{cursor:pointer;font-weight:700}
  .method-grid,.interpret-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:1rem;margin-top:1rem}
  .method-grid article,.interpret-grid article{border:1px solid var(--border-subtle);border-radius:9px;padding:1rem;background:var(--bg-secondary);min-width:0}
  dl{margin:.8rem 0}dl div{display:grid;grid-template-columns:90px 1fr;gap:.5rem;padding:.35rem 0;border-bottom:1px solid var(--border-subtle)}dt{color:var(--text-muted)}dd{margin:0;overflow-wrap:anywhere}
  .note,.cap{font-size:.88rem;color:var(--text-secondary);line-height:1.55}
  .cap{padding:.8rem 1rem;border-left:3px solid var(--accent-pk);background:var(--bg-secondary)}
  .result-block{margin:1.5rem 0 2.2rem}
  .qualified{font-weight:800}
  .reading-row{padding:.65rem 0;border-bottom:1px solid var(--border-subtle)}
  .reading-row:last-child{border-bottom:0}.reading-row p{margin:.3rem 0 0;line-height:1.5;font-size:.9rem}
  .viz{margin-top:1rem}
  .limit{margin:0 0 4rem;padding:1rem 1.2rem;border-left:4px solid var(--accent-ai);background:var(--bg-secondary)}
  .limit p{margin-bottom:0;color:var(--text-secondary)}
  .error{color:var(--danger,#b91c1c);font-weight:650}
  @media(max-width:820px){.method-grid,.interpret-grid{grid-template-columns:1fr}.facts{grid-template-columns:1fr 1fr}}
  @media(max-width:520px){.facts,.facts.small{grid-template-columns:1fr}.actions{align-items:stretch}.actions a{width:100%;text-align:center}}
</style>
