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
  /** @type {Record<string,string>} */
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
    if(value == null || value === '') return '—';
    const x=Number(value);
    return Number.isFinite(x) ? x.toPrecision(3) : '—';
  }
  /** @param {unknown} value */
  function probability(value){
    if(value == null || value === '') return t('non calculé','not computed');
    const x=Number(value);
    if(!Number.isFinite(x)) return t('non calculé','not computed');
    if(x===0) return '<1e−15';
    if(Math.abs(x)<0.001) return x.toExponential(2).replace('e-','e−');
    return x.toPrecision(3);
  }
  /** @param {unknown} value */
  function effect(value){
    return Number.isFinite(Number(value)) ? Number(value).toPrecision(4) : '—';
  }
  const demoNames = {
    transcriptomics: { STAT1:'STAT1', IDO1:'IDO1', KYNU:'KYNU', CYP3A5:'CYP3A5', GAPDH:'GAPDH' },
    proteomics: { P42224:'STAT1', P14902:'IDO1', Q16719:'KYNU', P04406:'GAPDH' },
    metabolomics: {
      'CHEBI:16828':t('L-tryptophane','L-tryptophan'),
      'CHEBI:16946':t('L-kynurénine','L-kynurenine'),
      'CHEBI:24996':t('Lactate','Lactate'),
      'CHEBI:30031':t('Succinate','Succinate')
    }
  };
  /** @param {string} layer @param {string} feature */
  function featureRow(layer,feature){
    return analysisResult?.layers?.[layer]?.rows?.find((row)=>row.feature===feature) || null;
  }
  /** @param {string} layer @param {string} feature */
  function featureName(layer,feature){
    return demoNames?.[layer]?.[feature] || feature;
  }
  /** @param {string} layer @param {string} feature */
  function featureValue(layer,feature){
    const row=featureRow(layer,feature);
    if(!row)return '—';
    return row.foldRatio==null ? effect(row.effect) : row.foldRatio.toFixed(2)+'×';
  }
  function topCrossPair(){
    return analysisResult?.crossOmics?.pairs?.[0] || null;
  }
  /** @param {any} row @param {any} layer */
  function rowReading(row,layer){
    const direction=Number(row?.effect)>0
      ? t('la trajectoire T0→T12 est plus positive sous traitement que sous contrôle','the T0→T12 trajectory is more positive under treatment than control')
      : Number(row?.effect)<0
        ? t('la trajectoire T0→T12 est plus négative sous traitement que sous contrôle','the T0→T12 trajectory is more negative under treatment than control')
        : t('aucune différence de trajectoire n’est estimée','no trajectory difference is estimated');
    const inference = layer?.inferenceStatus === 'reference_required_for_p_q_and_ci'
      ? t(
          `p/q de référence : non calculé dans le navigateur. Approximation exploratoire : p*=${probability(row?.approxPValue)}, q*=${probability(row?.approxQValue)}. Rejouez ${layer?.referenceMethodRequired || 'lmerTest'} sous R pour l’inférence de référence.`,
          `Reference p/q: not computed in the browser. Exploratory approximation: p*=${probability(row?.approxPValue)}, q*=${probability(row?.approxQValue)}. Rerun ${layer?.referenceMethodRequired || 'lmerTest'} in R for reference inference.`
        )
      : Number.isFinite(row?.qValue)
        ? (row.qValue<=0.10
            ? t('La variable franchit le seuil exploratoire q BH ≤ 0,10 utilisé par l’outil.','The feature passes the exploratory BH q ≤ 0.10 display threshold used by the tool.')
            : t('La variable ne franchit pas le seuil exploratoire q BH ≤ 0,10.','The feature does not pass the exploratory BH q ≤ 0.10 threshold.'))
        : t('Aucune q-value n’est publiée pour cette sortie.','No q-value is published for this output.');
    const ratio=layer?.effectScale==='log2' && Number.isFinite(row?.foldRatio)
      ? t(
          ` 2^effet = ${row.foldRatio.toFixed(2)}× : c’est le ratio des variations multiplicatives T0→T12, traitement par rapport au contrôle.`,
          ` 2^effect = ${row.foldRatio.toFixed(2)}×: this is the ratio of multiplicative T0→T12 changes, treatment relative to control.`
        )
      : '';
    return direction+'. '+inference+ratio;
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
    [t('Confiance d’identification métabolomique','Metabolomics identification confidence'),t('Mixte / inconnue','Mixed / unknown'),MULTIOMICS_DEMO_UI_STATE.metabolomicsIdentificationConfidence],
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
      /** @type {any} */
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
      'Le questionnaire décrit un plan longitudinal à mesures répétées. Tool sélectionne donc un modèle exploratoire à intercept aléatoire par sujet : feature ~ condition × temps + batch + (1|sujet). Le coefficient présenté est l’interaction condition × temps ramenée sur T0→T12. Les p/q de référence restent volontairement non calculées dans le navigateur car la calibration sous le nul était insuffisante ; l’outil expose désormais séparément p*/q* comme approximation exploratoire. Pour l’inférence de référence, le moteur demande lmerTest sous R.',
      'The questionnaire describes a repeated-measures longitudinal design. Tool therefore selects an exploratory subject random-intercept model: feature ~ condition × time + batch + (1|subject). The displayed coefficient is the condition × time interaction over T0→T12. Reference p/q remain intentionally not computed in the browser because null calibration was insufficient; the tool now exposes separate p*/q* as exploratory approximations. For reference inference, the engine requires lmerTest in R.'
    )}</p>
    <div class="method-grid" data-testid="guided-demo-methods">
      {#each layerOrder as layer}
        {@const r=analysisResult.layers[layer]}
        <article>
          <p class="eyebrow">{layerLabel(layer)}</p>
          <h3>{r.contrast}</h3>
          <dl>
            <div><dt>{t('Mode','Mode')}</dt><dd><code>{r.mode}</code></dd></div>
            <div><dt>{t('Méthode','Method')}</dt><dd>{r.inferenceMethod}</dd></div>
            <div><dt>n</dt><dd>{r.groupSizes?.[0]} vs {r.groupSizes?.[1]}</dd></div>
            <div><dt>{t('Échelle','Scale')}</dt><dd><code>{r.effectScale}</code></dd></div>
            <div><dt>{t('Statut inférentiel','Inferential status')}</dt><dd><code>{r.inferenceStatus || '—'}</code></dd></div>
            <div><dt>{t('Référence R','R reference')}</dt><dd><code>{r.referenceMethodRequired || '—'}</code></dd></div>
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
    'Les p/q de référence sont affichées « non calculé » lorsqu’elles ne sont pas validées dans le navigateur. Les colonnes p*/q* montrent séparément l’approximation navigateur, avec correction BH appliquée aux p*. Ces valeurs exploratoires sont affichées en notation scientifique lorsqu’elles sont petites et doivent être rejouées sous R pour l’inférence de référence.',
    'Reference p/q are shown as “not computed” when they are not validated in the browser. Separate p*/q* columns show the browser approximation, with BH applied to p*. Small exploratory values use scientific notation and must be rerun in R for reference inference.'
  )}</p>
  {#if status==='done'}
    {#each layerOrder as layer}
      {@const r=analysisResult.layers[layer]}
      <div class="result-block" data-testid={'guided-demo-layer-'+layer}>
        <h3>{layerLabel(layer)} — {r.contrast}</h3>
        <div class="table-wrap">
          <table>
            <thead>
              <tr><th>{t('Variable','Feature')}</th><th>{t('Effet','Effect')}</th><th>{t('Échelle','Scale')}</th><th>2^effect</th><th>{t('p référence','reference p')}</th><th>{t('q BH référence','reference BH q')}</th><th>p* {t('approx.','approx.')}</th><th>q* BH {t('approx.','approx.')}</th><th>n ref</th><th>n comp</th></tr>
            </thead>
            <tbody>
              {#each r.rows as row}
                <tr data-feature={row.feature}>
                  <td><code>{row.feature}</code></td>
                  <td>{effect(row.effect)}</td>
                  <td><code>{row.effectScale}</code></td>
                  <td>{row.foldRatio==null?'—':row.foldRatio.toFixed(2)+'×'}</td>
                  <td class="not-computed">{probability(row.pValue)}</td>
                  <td class="not-computed">{probability(row.qValue)}</td>
                  <td>{probability(row.approxPValue)}</td>
                  <td class:qualified={Number.isFinite(row.approxQValue)&&row.approxQValue<=0.10}>{probability(row.approxQValue)}</td>
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
    <p class="section-intro">{t(
      'Une bonne lecture ne commence pas par chercher la plus petite p-value. Elle suit quatre questions successives : quel contraste est estimé, quelle est l’amplitude, le signal se répète-t-il entre couches, puis quelle hypothèse biologique ou clinique mérite réellement d’être testée ?',
      'A sound interpretation does not start by hunting for the smallest p-value. It follows four questions in sequence: what contrast is estimated, how large is the effect, is the signal repeated across layers, and which biological or clinical hypothesis is actually worth testing?'
    )}</p>

    <ol class="reading-sequence" data-testid="guided-demo-interpretation">
      <li>
        <span>1</span>
        <div><strong>{t('Lire le contraste','Read the contrast')}</strong><p>{t(
          'Ici, chaque effet est une interaction condition × temps : différence entre la trajectoire T0→T12 du groupe traité et celle du groupe contrôle. Ce n’est pas une simple différence à T12.',
          'Here every effect is a condition × time interaction: the difference between the treated and control T0→T12 trajectories. It is not a simple T12 group difference.'
        )}</p></div>
      </li>
      <li>
        <span>2</span>
        <div><strong>{t('Lire l’amplitude avant la significativité','Read magnitude before significance')}</strong><p>{t(
          'Sur une échelle log2, 2^effet donne un ratio intuitif. Par exemple IDO1 vaut '+featureValue('transcriptomics','IDO1')+' au niveau RNA et '+featureValue('proteomics','P14902')+' au niveau protéique dans ce cas synthétique.',
          'On a log2 scale, 2^effect gives an intuitive ratio. In this synthetic case IDO1 is '+featureValue('transcriptomics','IDO1')+' at RNA level and '+featureValue('proteomics','P14902')+' at protein level.'
        )}</p></div>
      </li>
      <li>
        <span>3</span>
        <div><strong>{t('Vérifier la cohérence entre couches','Check cross-layer coherence')}</strong><p>{t(
          'Un même axe est plus crédible lorsqu’il est visible à plusieurs niveaux indépendants : ici STAT1, IDO1 et KYNU évoluent dans le même sens en transcriptomique et protéomique, tandis que le tryptophane et la kynurénine apportent une lecture métabolique complémentaire.',
          'A biological axis is more compelling when it appears across independent layers: here STAT1, IDO1 and KYNU move in the same direction in transcriptomics and proteomics, while tryptophan and kynurenine provide complementary metabolic readout.'
        )}</p></div>
      </li>
      <li>
        <span>4</span>
        <div><strong>{t('Séparer preuve et hypothèse','Separate evidence from hypothesis')}</strong><p>{t(
          'Les p*/q* navigateur servent à hiérarchiser les signaux. Elles ne remplacent pas le replay R. Une histoire biologique cohérente reste une hypothèse tant qu’elle n’est pas reproduite, reliée au phénotype et confrontée à des mesures fonctionnelles.',
          'Browser p*/q* help rank signals. They do not replace the R replay. A coherent biological story remains a hypothesis until it is replicated, linked to phenotype and challenged with functional measurements.'
        )}</p></div>
      </li>
    </ol>

    <div class="interpretation-story">
      <h3>{t('Interprétation scientifique du cas','Scientific interpretation of the case')}</h3>
      <p>{t(
        'Le motif dominant est compatible avec une activation coordonnée de l’axe interféron–kynurénine. STAT1 augmente au niveau RNA ('+featureValue('transcriptomics','STAT1')+') et protéique ('+featureValue('proteomics','P42224')+'). IDO1 augmente également dans les deux couches ('+featureValue('transcriptomics','IDO1')+' RNA ; '+featureValue('proteomics','P14902')+' protéine), de même que KYNU ('+featureValue('transcriptomics','KYNU')+' ; '+featureValue('proteomics','Q16719')+').',
        'The dominant pattern is compatible with coordinated activation of an interferon–kynurenine axis. STAT1 increases at RNA ('+featureValue('transcriptomics','STAT1')+') and protein ('+featureValue('proteomics','P42224')+') levels. IDO1 also increases in both layers ('+featureValue('transcriptomics','IDO1')+' RNA; '+featureValue('proteomics','P14902')+' protein), as does KYNU ('+featureValue('transcriptomics','KYNU')+'; '+featureValue('proteomics','Q16719')+').'
      )}</p>
      <p>{t(
        'Le métabolome va dans le même sens qualitatif : le '+featureName('metabolomics','CHEBI:16828')+' est à '+featureValue('metabolomics','CHEBI:16828')+' tandis que la '+featureName('metabolomics','CHEBI:16946')+' est à '+featureValue('metabolomics','CHEBI:16946')+'. Cette convergence est compatible avec une mobilisation de la voie kynurénine, mais elle ne démontre pas un flux métabolique : il faudrait notamment un ratio kynurénine/tryptophane validé, des concentrations quantitatives et idéalement des mesures fonctionnelles ou isotopiques.',
        'The metabolome points in the same qualitative direction: '+featureName('metabolomics','CHEBI:16828')+' is '+featureValue('metabolomics','CHEBI:16828')+' while '+featureName('metabolomics','CHEBI:16946')+' is '+featureValue('metabolomics','CHEBI:16946')+'. This convergence is compatible with engagement of the kynurenine pathway, but it does not demonstrate metabolic flux: a validated kynurenine/tryptophan ratio, quantitative concentrations and ideally functional or isotope measurements would be needed.'
      )}</p>
    </div>

    <div class="translation">
      <h3>{t('Interprétation clinique / translationnelle','Clinical / translational interpretation')}</h3>
      <p>{t(
        'Si un profil analogue était reproduit chez de vrais patients, il ferait d’abord suspecter un état inflammatoire ou immunomodulé associé à une activation de la voie kynurénine. Ce serait une raison de confronter les omiques à des marqueurs cliniques indépendants — CRP, cytokines, symptômes, réponse au traitement — plutôt que de conclure à partir des omiques seules.',
        'If a similar profile were replicated in real patients, it would first raise the hypothesis of an inflammatory or immunomodulated state associated with kynurenine-pathway activation. The omics signal should then be confronted with independent clinical markers—CRP, cytokines, symptoms and treatment response—rather than treated as a conclusion by itself.'
      )}</p>
      <p>{t(
        'Le signal CYP3A5 ajoute une hypothèse pharmacologique : son expression RNA est ici à '+featureValue('transcriptomics','CYP3A5')+'. Dans une vraie cohorte exposée à un médicament substrat de CYP3A, cela justifierait de tester une phénoconversion inflammatoire possible avec concentrations PK, génotype CYP3A5, activité/phénotype CYP3A et covariables cliniques. Cela ne justifie jamais, à lui seul, une adaptation posologique.',
        'CYP3A5 adds a pharmacological hypothesis: its RNA expression is '+featureValue('transcriptomics','CYP3A5')+' here. In a real cohort receiving a CYP3A substrate, this would justify testing possible inflammation-related phenoconversion using PK concentrations, CYP3A5 genotype, CYP3A activity/phenotype and clinical covariates. By itself it never justifies a dose change.'
      )}</p>
      <div class="translation-next">
        <strong>{t('Étape suivante logique','Logical next step')}</strong>
        <span>{t(
          'Relier ce profil individuel à la PK : exposition/AUC ou concentrations → score moléculaire ou variables omiques → comparaison contre clinique + pharmacogénétique. C’est précisément le type de pont destiné à Omics2PK.',
          'Link this individual profile to PK: exposure/AUC or concentrations → molecular score or omics variables → comparison against clinical + pharmacogenetic predictors. This is precisely the type of bridge intended for Omics2PK.'
        )}</span>
      </div>
    </div>

    <p class="source-note">{t('Contexte biologique externe :','External biological context:')}
      <a href="https://www.uniprot.org/uniprotkb/P42224" target="_blank" rel="noopener noreferrer">STAT1 / P42224 ↗</a>,
      <a href="https://www.uniprot.org/uniprotkb/P14902" target="_blank" rel="noopener noreferrer">IDO1 / P14902 ↗</a>,
      <a href="https://www.uniprot.org/uniprotkb/Q16719" target="_blank" rel="noopener noreferrer">KYNU / Q16719 ↗</a>,
      <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8637893/" target="_blank" rel="noopener noreferrer">{t('inflammation et CYP','inflammation and CYP')} ↗</a>.
    </p>
    <p class="cap">{t(
      'p/q de référence restent non calculées dans cette exécution navigateur ; p*/q* sont uniquement des approximations exploratoires. Les conclusions quantitatives d’un article doivent être établies avec le replay R lmerTest et une validation biologique indépendante.',
      'Reference p/q remain not computed in this browser run; p*/q* are exploratory approximations only. Quantitative conclusions for a paper should be established with the R lmerTest replay and independent biological validation.'
    )}</p>
  {/if}
</section>

<section>
  <h2>{t('6 · Comprendre chaque figure','6 · Understand each figure')}</h2>
  <p class="section-intro">{t(
    'Chaque figure répond à une question différente. Elles utilisent les mêmes résultats que Tool mais ne constituent pas trois tests statistiques supplémentaires.',
    'Each figure answers a different question. They reuse the same Tool results and are not three additional statistical tests.'
  )}</p>

  <div class="output-guide" data-testid="guided-demo-figure-guide">
    <div class="output-guide-head">
      <span>{t('Sortie','Output')}</span><span>{t('Exemple de sortie','Example output')}</span><span>{t('Méthode','Method')}</span><span>{t('Interprétation / utilité','Interpretation / utility')}</span><span>{t('Limite','Limit')}</span>
    </div>
    <div class="output-guide-row">
      <strong>Heatmap</strong>
      <span>{t('Sujets T0/T12 ordonnés avec les variables les plus informatives.','T0/T12 specimens ordered with the most informative features.')}</span>
      <span>{t('Valeurs d’affichage transformées ; réplicats techniques moyennés ; z-score par variable.','Display-transformed values; technical replicates averaged; per-feature z-score.')}</span>
      <span>{t('Vérifier la structure globale, les groupes/temps, les réplicats et les éventuels échantillons atypiques.','Check global structure, groups/time, replicates and potential outlier specimens.')}</span>
      <span>{t('La couleur n’est ni un fold change ni une p-value ; une séparation visuelle n’est pas un test.','Color is neither fold change nor p-value; visual separation is not a test.')}</span>
    </div>
    <div class="output-guide-row">
      <strong>{t('Metabologramme','Metabologram')}</strong>
      <span>STAT1 / IDO1 / KYNU ↔ {featureName('metabolomics','CHEBI:16828')} / {featureName('metabolomics','CHEBI:16946')}</span>
      <span>{t('Effets RNA et métabolites replacés dans des régions biologiques connues ; amplitude normalisée séparément par couche.','RNA and metabolite effects placed in known biological regions; magnitude scaled separately by layer.')}</span>
      <span>{t('Repérer une convergence de direction entre couches et choisir les voies à investiguer.','Identify directional convergence across layers and prioritize pathways for investigation.')}</span>
      <span>{t('Présence sur une voie ≠ enrichissement, activité enzymatique ou flux.','Pathway presence ≠ enrichment, enzyme activity or flux.')}</span>
    </div>
    <div class="output-guide-row">
      <strong>{t('Réseau ciblé','Focused network')}</strong>
      <span>{t('La région tryptophane–kynurénine est sélectionnée automatiquement ici car elle contient plusieurs signaux mesurés.','The tryptophan–kynurenine region is selected automatically here because it contains several measured signals.')}</span>
      <span>{t('Appariement strict des identifiants vers un catalogue curé, sélection automatique des régions mesurées, puis affichage du voisinage à 0–2 sauts.','Strict identifier matching to a curated catalog, automatic selection of measured regions, then 0–2-hop neighborhood display.')}</span>
      <span>{t('Comprendre où se situent les signaux mesurés et quels voisins biologiques donnent du contexte.','Understand where measured signals sit and which biological neighbors provide context.')}</span>
      <span>{t('Catalogue non exhaustif ; une voie absente n’est pas inventée. Protéomique non projetée directement à ce jour.','Non-exhaustive catalog; an absent pathway is never invented. Proteomics is not directly projected yet.')}</span>
    </div>
  </div>

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
    {@const pair=topCrossPair()}
    <p class="section-intro">{t(
      'Cette sortie ne cherche pas seulement à savoir quelles variables changent. Elle demande si la relation entre deux variables de couches différentes change elle aussi entre contrôle et traitement.',
      'This output does not only ask which features change. It asks whether the relationship between two variables from different layers also changes between control and treatment.'
    )}</p>

    <div class="output-explainer" data-testid="guided-demo-crossomics-explainer">
      <div>
        <span class="output-label">{t('Exemple de sortie','Example output')}</span>
        {#if pair}
          <p><code>{featureName(pair.layerA,pair.featureA)}</code> ↔ <code>{featureName(pair.layerB,pair.featureB)}</code><br>
          r<sub>contrôle</sub>={sig(pair.rReference)} · r<sub>traitement</sub>={sig(pair.rComparison)} · Δr={sig(pair.deltaR)} · p*={probability(pair.approxPValue)} · q*={probability(pair.approxQValue)}</p>
        {:else}
          <p>{t('Aucune paire estimable dans ce cas.','No estimable pair in this case.')}</p>
        {/if}
      </div>
      <div>
        <span class="output-label">{t('Méthode','Method')}</span>
        <p>{t(
          'Pour chaque paire candidate, MultiOmics résume la trajectoire individuelle puis calcule une corrélation de Spearman séparément dans les deux groupes. Δr = r_traitement − r_contrôle. L’approximation Fisher-z produit p*, puis BH produit q* sur les paires effectivement testées.',
          'For each candidate pair, MultiOmics summarizes the individual trajectory and calculates Spearman correlation separately in both groups. Δr = r_treatment − r_control. A Fisher-z approximation yields p*, then BH yields q* across the tested pairs.'
        )}</p>
      </div>
      <div>
        <span class="output-label">{t('Interprétation / utilité','Interpretation / utility')}</span>
        <p>{t(
          'Un |Δr| élevé suggère que deux variables deviennent plus ou moins couplées sous traitement. Cette sortie sert surtout à prioriser des relations inter-couches à vérifier expérimentalement, ou à construire une signature mécanistique candidate.',
          'A large |Δr| suggests that two variables become more or less coupled under treatment. This output is mainly useful for prioritizing cross-layer relationships for experimental follow-up or for building a candidate mechanistic signature.'
        )}</p>
      </div>
      <div>
        <span class="output-label">{t('Limites','Limitations')}</span>
        <p>{t(
          'Une corrélation différentielle ne prouve pas que A régule B. Ici les p/q de référence ne sont pas calculées dans le navigateur, l’effectif est faible et les paires sont présélectionnées parmi les variables les plus informatives. Les p*/q* servent à l’exploration, pas à une revendication confirmatoire.',
          'Differential correlation does not prove that A regulates B. Here reference p/q are not computed in the browser, sample size is small and pairs are preselected among informative features. p*/q* are exploratory, not confirmatory.'
        )}</p>
      </div>
    </div>

    <div class="facts small">
      <article><strong>{analysisResult.crossOmics?.testedPairs ?? 0}</strong><span>{t('paires testées','pairs tested')}</span></article>
      <article><strong>{analysisResult.crossOmics?.significantPairs ?? 0}</strong><span>{t('paires avec q référence ≤ 0,10','pairs with reference q ≤ 0.10')}</span></article>
      <article><strong>{analysisResult.crossOmics?.approximateSignificantPairs ?? 0}</strong><span>{t('paires avec q* approx. ≤ 0,10','pairs with approx. q* ≤ 0.10')}</span></article>
    </div>

    <details class="technical-pairs">
      <summary>{t('Voir le tableau détaillé des paires inter-omiques','Show detailed cross-omics pair table')}</summary>
      {#if analysisResult.crossOmics?.pairs?.length}
        <div class="table-wrap">
          <table>
            <thead><tr><th>{t('Paire','Pair')}</th><th>Pattern</th><th>r ref</th><th>r comp</th><th>Δr</th><th>{t('p référence','reference p')}</th><th>{t('q référence','reference q')}</th><th>p* {t('approx.','approx.')}</th><th>q* BH {t('approx.','approx.')}</th></tr></thead>
            <tbody>
              {#each analysisResult.crossOmics.pairs.slice(0,15) as pairRow}
                <tr>
                  <td><code>{featureName(pairRow.layerA,pairRow.featureA)}</code> ↔ <code>{featureName(pairRow.layerB,pairRow.featureB)}</code></td>
                  <td>{pairRow.pattern?.replaceAll('_',' ')}</td>
                  <td>{sig(pairRow.rReference)}</td><td>{sig(pairRow.rComparison)}</td><td>{sig(pairRow.deltaR)}</td><td>{probability(pairRow.pValue)}</td><td>{probability(pairRow.qValue)}</td><td>{probability(pairRow.approxPValue)}</td><td>{probability(pairRow.approxQValue)}</td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
      {:else}
        <p>{t('Aucune paire inter-omique estimable.','No estimable cross-omics pair.')}</p>
      {/if}
    </details>
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
  .facts.small{grid-template-columns:repeat(3,minmax(0,220px))}
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
  .method-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:1rem;margin-top:1rem}
  .method-grid article{border:1px solid var(--border-subtle);border-radius:9px;padding:1rem;background:var(--bg-secondary);min-width:0}
  .section-intro{max-width:82ch;color:var(--text-secondary);line-height:1.7}
  .reading-sequence{list-style:none;padding:0;margin:1.4rem 0 2rem;display:grid;gap:.9rem}
  .reading-sequence li{display:grid;grid-template-columns:42px 1fr;gap:.8rem;align-items:start}
  .reading-sequence li>span{display:grid;place-items:center;width:34px;height:34px;border:1px solid var(--border-strong);border-radius:50%;font-family:var(--font-mono);font-weight:800;color:var(--accent-pk)}
  .reading-sequence p{margin:.25rem 0 0;color:var(--text-secondary);line-height:1.6}
  .interpretation-story,.translation{max-width:900px;margin:1.5rem 0;padding:1.1rem 1.25rem;border-left:4px solid var(--accent-pk);background:var(--bg-secondary)}
  .translation{border-left-color:var(--accent-ai)}
  .interpretation-story p,.translation p{line-height:1.7}
  .translation-next{display:grid;gap:.25rem;margin-top:1rem;padding-top:.8rem;border-top:1px solid var(--border-subtle)}
  .source-note{font-size:.82rem;color:var(--text-muted);line-height:1.6}
  .source-note a{color:inherit}
  .output-guide{margin:1.2rem 0;border:1px solid var(--border-subtle);border-radius:10px;overflow:auto}
  .output-guide-head,.output-guide-row{display:grid;grid-template-columns:.8fr 1.2fr 1.4fr 1.5fr 1.3fr;min-width:980px}
  .output-guide-head{background:var(--bg-secondary);font-size:.75rem;text-transform:uppercase;letter-spacing:.04em;color:var(--text-secondary);font-weight:750}
  .output-guide-head span,.output-guide-row>*{padding:.8rem;border-right:1px solid var(--border-subtle);border-bottom:1px solid var(--border-subtle);line-height:1.5}
  .output-guide-row:last-child>*{border-bottom:0}
  .output-guide-row>*:last-child,.output-guide-head span:last-child{border-right:0}
  .output-guide-row span{font-size:.88rem;color:var(--text-secondary)}
  .output-explainer{display:grid;gap:0;border:1px solid var(--border-subtle);border-radius:10px;overflow:hidden;margin:1.2rem 0}
  .output-explainer>div{display:grid;grid-template-columns:180px 1fr;gap:1rem;padding:1rem;border-bottom:1px solid var(--border-subtle)}
  .output-explainer>div:last-child{border-bottom:0}
  .output-label{font-weight:750;color:var(--text-primary)}
  .output-explainer p{margin:0;line-height:1.65;color:var(--text-secondary)}
  .technical-pairs{margin-top:1.2rem}
  .not-computed{color:var(--text-muted);font-style:italic}
  dl{margin:.8rem 0}dl div{display:grid;grid-template-columns:90px 1fr;gap:.5rem;padding:.35rem 0;border-bottom:1px solid var(--border-subtle)}dt{color:var(--text-muted)}dd{margin:0;overflow-wrap:anywhere}
  .note,.cap{font-size:.88rem;color:var(--text-secondary);line-height:1.55}
  .cap{padding:.8rem 1rem;border-left:3px solid var(--accent-pk);background:var(--bg-secondary)}
  .result-block{margin:1.5rem 0 2.2rem}
  .qualified{font-weight:800}
  .viz{margin-top:1rem}
  .limit{margin:0 0 4rem;padding:1rem 1.2rem;border-left:4px solid var(--accent-ai);background:var(--bg-secondary)}
  .limit p{margin-bottom:0;color:var(--text-secondary)}
  .error{color:var(--danger,#b91c1c);font-weight:650}
  @media(max-width:820px){.method-grid{grid-template-columns:1fr}.facts{grid-template-columns:1fr 1fr}.output-explainer>div{grid-template-columns:1fr;gap:.35rem}}
  @media(max-width:520px){.facts,.facts.small{grid-template-columns:1fr}.actions{align-items:stretch}.actions a{width:100%;text-align:center}}
</style>
