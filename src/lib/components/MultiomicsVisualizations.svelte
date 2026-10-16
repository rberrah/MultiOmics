<script>
  // @ts-nocheck
  import FocusedMetabolicNetwork from './FocusedMetabolicNetwork.svelte';
  export let result = null;
  export let language = 'fr';

  const tr = (fr, en) => language === 'en' ? en : fr;
  let heatmapLayer = 'transcriptomics';
  let activeFigure = 'heatmap';
  let selectedPathway = 'all';

  $: visuals = result?.visualizations || null;
  $: heatmapLayers = visuals ? Object.keys(visuals.heatmaps || {}).filter((layer) => visuals.heatmaps[layer]?.rows?.length) : [];
  $: if (heatmapLayers.length && !heatmapLayers.includes(heatmapLayer)) heatmapLayer = heatmapLayers[0];
  $: heatmap = visuals?.heatmaps?.[heatmapLayer] || null;
  $: metabologram = visuals?.metabologram || null;
  $: central = visuals?.centralCarbon || null;
  $: pathwayChoices = visuals?.metabologramPathways || [];
  $: activeMetabologram = selectedPathway === 'all'
    ? metabologram
    : (pathwayChoices.find((pathway) => pathway.id === selectedPathway) || metabologram);
  $: pathwayCoverage = selectedPathway === 'all' ? null
    : pathwayChoices.find((pathway) => pathway.id === selectedPathway);
  $: metLimit = effectLimit(activeMetabologram?.metabolomics);
  $: rnaLimit = effectLimit(activeMetabologram?.transcriptomics);
  $: mapMetLimit = effectLimit(central?.metabolites?.map((node) => node.measurement).filter(Boolean));
  $: mapRnaLimit = effectLimit(central?.enzymes?.map((node) => node.measurement).filter(Boolean));

  function layerLabel(layer) {
    if (layer === 'transcriptomics') return tr('Transcriptomique', 'Transcriptomics');
    if (layer === 'proteomics') return tr('Protéomique', 'Proteomics');
    if (layer === 'metabolomics') return tr('Métabolomique', 'Metabolomics');
    return layer;
  }

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function effectLimit(rows) {
    const source = (rows || []).filter((item) => Number.isFinite(item?.effect));
    if (!source.length) return 3;
    if (source.every((item) => item.effectScale === 'log2')) return 3;
    return Math.max(1e-9, ...source.map((item) => Math.abs(item.effect)));
  }

  function scaleLabel(rows) {
    const scales = [...new Set((rows || []).map((item) => item.effectScale).filter(Boolean))];
    if (!scales.length) return '—';
    if (scales.length > 1) return tr('échelles mixtes', 'mixed scales');
    return scales[0] === 'log2' ? 'log2FC'
      : scales[0] === 'transformed_unknown'
        ? tr('différence sur échelle transformée inconnue', 'difference on unknown transformed scale')
        : tr('différence sur échelle fournie', 'difference on supplied scale');
  }

  function diverging(value, limit = 3) {
    if (!Number.isFinite(value)) return 'rgb(224,224,224)';
    const x = clamp(value / limit, -1, 1);
    const white = [247, 247, 247];
    const red = [204, 62, 62];
    const blue = [53, 103, 178];
    const target = x >= 0 ? red : blue;
    const amount = Math.abs(x);
    return 'rgb(' + white.map((channel, i) => Math.round(channel + (target[i] - channel) * amount)).join(',') + ')';
  }

  function qStroke(q) {
    if (Number.isFinite(q) && q <= 0.05) return { width: 3, dash: '' };
    if (Number.isFinite(q) && q <= 0.10) return { width: 2, dash: '' };
    return { width: 1, dash: '3 2' };
  }

  function fmt(value, digits = 2) {
    return Number.isFinite(value) ? Number(value).toFixed(digits) : '—';
  }

  function itemLabel(item) {
    return language === 'en' ? (item.labelEn || item.feature) : (item.labelFr || item.feature);
  }

  function shortFeature(value, max = 16) {
    const x = String(value || '');
    return x.length > max ? x.slice(0, max - 1) + '…' : x;
  }

  function polar(cx, cy, radius, angle) {
    const radians = angle * Math.PI / 180;
    return { x: cx + radius * Math.cos(radians), y: cy + radius * Math.sin(radians) };
  }

  function donutPath(cx, cy, innerRadius, outerRadius, startAngle, endAngle) {
    const outerStart = polar(cx, cy, outerRadius, startAngle);
    const outerEnd = polar(cx, cy, outerRadius, endAngle);
    const innerEnd = polar(cx, cy, innerRadius, endAngle);
    const innerStart = polar(cx, cy, innerRadius, startAngle);
    const large = Math.abs(endAngle - startAngle) > 180 ? 1 : 0;
    const sweep = endAngle > startAngle ? 1 : 0;
    const innerSweep = sweep ? 0 : 1;
    return [
      'M', outerStart.x, outerStart.y,
      'A', outerRadius, outerRadius, 0, large, sweep, outerEnd.x, outerEnd.y,
      'L', innerEnd.x, innerEnd.y,
      'A', innerRadius, innerRadius, 0, large, innerSweep, innerStart.x, innerStart.y,
      'Z'
    ].join(' ');
  }

  function ringSegments(entries, side) {
    if (!entries?.length) return [];
    const start = side === 'left' ? 90 : -90;
    const end = side === 'left' ? 270 : 90;
    const step = (end - start) / entries.length;
    return entries.map((entry, index) => {
      const a0 = start + index * step + 0.6;
      const a1 = start + (index + 1) * step - 0.6;
      return { ...entry, path: donutPath(250, 250, 112, 190, a0, a1) };
    });
  }

  $: metaboliteSegments = ringSegments(activeMetabologram?.metabolomics || [], 'left');
  $: transcriptSegments = ringSegments(activeMetabologram?.transcriptomics || [], 'right');

  function findCentralNode(id) {
    return central?.metabolites?.find((node) => node.id === id);
  }

  function measurementColor(measurement, limit = 3) {
    return measurement?.status === 'measured' ? diverging(measurement.effect, limit) : 'rgb(220,220,220)';
  }

  function measurementTitle(label, measurement) {
    if (!measurement) return label + ' · ' + tr('non mesuré / non reconnu', 'not measured / not matched');
    if (measurement.status !== 'measured') return label + ' · ' + measurement.feature + ' · ' + measurement.status;
    return label + ' · ' + measurement.feature + ' · ' + scaleLabel([measurement]) + '=' + fmt(measurement.effect) + ' · q=' + fmt(measurement.qValue, 3);
  }

  function downloadSvg(id, filename) {
    const source = document.getElementById(id);
    if (!(source instanceof SVGElement)) return;
    const copy = source.cloneNode(true);
    copy.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    const xml = new XMLSerializer().serializeToString(copy);
    const blob = new Blob([xml], { type: 'image/svg+xml;charset=utf-8' });
    const href = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = href;
    anchor.download = filename;
    anchor.click();
    URL.revokeObjectURL(href);
  }
</script>

{#if visuals}
<section class="visual-panel" data-testid="multiomics-visualizations">
  <div class="visual-head">
    <div>
      <p class="eyebrow">{tr('Visualisation intégrée', 'Integrated visualization')}</p>
      <h2>{tr('Voir le signal biologique sous plusieurs angles', 'View the biological signal from several angles')}</h2>
      <p>{tr(
        'Ces figures réutilisent les résultats calculés. Elles n’ajoutent aucun test statistique et ne changent pas les q-values.',
        'These figures reuse the computed results. They add no statistical test and do not change q-values.'
      )}</p>
    </div>
  </div>

  <nav class="figure-switcher" aria-label={tr('Choisir une visualisation', 'Choose a visualization')}>
    <button type="button" class:active={activeFigure === 'heatmap'} aria-current={activeFigure === 'heatmap' ? 'true' : undefined} on:click={() => activeFigure = 'heatmap'}>
      {tr('1. Heatmap', '1. Heatmap')}
    </button>
    <button type="button" class:active={activeFigure === 'pathways'} aria-current={activeFigure === 'pathways' ? 'true' : undefined} on:click={() => activeFigure = 'pathways'}>
      {tr('2. Voies biologiques', '2. Biological pathways')}
    </button>
    <button type="button" class:active={activeFigure === 'map'} aria-current={activeFigure === 'map' ? 'true' : undefined} on:click={() => activeFigure = 'map'}>
      {tr('3. Réseau ciblé', '3. Focused network')}
    </button>
  </nav>

  {#if activeFigure === 'heatmap' && heatmapLayers.length}
    <article class="figure-card" data-testid="multiomics-heatmap">
      <div class="figure-title">
        <div>
          <span class="figure-number">01</span>
          <div>
            <h3>{tr('Heatmap des variables principales', 'Heatmap of leading features')}</h3>
            <p>{tr(
              'Chaque ligne est centrée-réduite séparément : la couleur montre une valeur relative au profil moyen de cette variable, pas un fold change.',
              'Each row is standardized separately: color shows a value relative to that feature’s mean profile, not a fold change.'
            )}</p>
          </div>
        </div>
        <div class="figure-actions">
          <select bind:value={heatmapLayer} aria-label={tr('Couche de la heatmap', 'Heatmap layer')}>
            {#each heatmapLayers as layer}
              <option value={layer}>{layerLabel(layer)}</option>
            {/each}
          </select>
          <button type="button" on:click={() => downloadSvg('pmx-omics-heatmap-svg', 'multiomics_heatmap.svg')}>SVG</button>
        </div>
      </div>

      {#if heatmap?.rows?.length}
        {@const cellW = 23}
        {@const rowH = 21}
        {@const left = 185}
        {@const top = 118}
        {@const width = Math.max(760, left + heatmap.samples.length * cellW + 35)}
        {@const height = top + heatmap.rows.length * rowH + 55}
        <div class="svg-scroll">
          <svg id="pmx-omics-heatmap-svg" viewBox={'0 0 ' + width + ' ' + height} role="img" aria-label={tr('Heatmap multi-omique', 'Multi-omics heatmap')}>
            <rect width={width} height={height} fill="white"/>
            <text x="12" y="24" class="svg-heading">{layerLabel(heatmapLayer)} · z-score</text>
            <text x="12" y="43" class="svg-note">{heatmap.selectionRule}</text>

            {#each heatmap.samples as sample, column}
              <text
                x={left + column * cellW + cellW / 2}
                y={top - 9}
                transform={'rotate(-58 ' + (left + column * cellW + cellW / 2) + ' ' + (top - 9) + ')'}
                text-anchor="start"
                class="sample-label"
              >{sample.sampleId}</text>
              <text x={left + column * cellW + cellW / 2} y={top - 2} text-anchor="middle" class="condition-mark">
                {sample.condition ? sample.condition.slice(0, 1) : '·'}
              </text>
            {/each}

            {#each heatmap.rows as row, rowIndex}
              <text x={left - 9} y={top + rowIndex * rowH + 15} text-anchor="end" class="feature-label">{shortFeature(row.feature, 24)}</text>
              {#each row.z as value, column}
                <rect
                  x={left + column * cellW}
                  y={top + rowIndex * rowH}
                  width={cellW - 1}
                  height={rowH - 1}
                  fill={diverging(value, 2.5)}
                >
                  <title>{row.feature} · {heatmap.samples[column].sampleId} · z={fmt(value)} · {tr('effet', 'effect')}={fmt(row.effect)} · q={fmt(row.qValue, 3)}</title>
                </rect>
              {/each}
            {/each}

            <text x={left} y={height - 18} class="svg-note">−2.5</text>
            {#each [-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5] as value, i}
              <rect x={left + 32 + i * 16} y={height - 31} width="16" height="13" fill={diverging(value, 2.5)}/>
            {/each}
            <text x={left + 218} y={height - 18} class="svg-note">+2.5 z</text>
          </svg>
        </div>
        <p class="method-note">{heatmap.method} {heatmap.samplesCapped ? tr('Seuls les 48 premiers échantillons ordonnés sont affichés.', 'Only the first 48 ordered samples are displayed.') : ''}</p>
      {/if}
    </article>
  {/if}

  {#if activeFigure === 'heatmap' && !heatmapLayers.length}
    <p class="empty-note">{tr('Aucune matrice utilisable pour la heatmap. Les résultats numériques restent accessibles ci-dessous.', 'No matrix available for a heatmap. Numeric results remain available below.')}</p>
  {/if}

  {#if activeFigure === 'pathways'}
  <article class="figure-card" data-testid="multiomics-metabologram">
    <div class="figure-title">
      <div>
        <span class="figure-number">02</span>
        <div>
          <h3>{tr('Metabologramme transcriptome ↔ métabolome', 'Transcriptome ↔ metabolome metabologram')}</h3>
          <p>{tr(
            'À gauche : métabolites ; à droite : gènes exprimés. Bleu = baisse ; rouge = hausse. Sélectionnez une voie biologique si besoin.',
            'Metabolites on the left; gene transcripts on the right. Blue = decrease; red = increase. Select a biological pathway if needed.'
          )}</p>
        </div>
      </div>
      <div class="figure-actions">
        <label class="pathway-select">
          <span>{tr('Voie affichée', 'Displayed pathway')}</span>
          <select bind:value={selectedPathway} data-testid="multiomics-pathway-select" aria-label={tr('Choisir la voie biologique', 'Choose biological pathway')}>
            <option value="all">{tr('Toutes les variables (sans filtre de voie)', 'All features (no pathway filter)')}</option>
            {#each pathwayChoices as pathway}
              <option value={pathway.id}>{language === 'en' ? pathway.labelEn : pathway.labelFr}</option>
            {/each}
          </select>
        </label>
        <button type="button" on:click={() => downloadSvg('pmx-metabologram-svg', 'multiomics_metabologram.svg')}>SVG</button>
      </div>
    </div>

    {#if metaboliteSegments.length || transcriptSegments.length}
      <div class="metabologram-layout">
        <svg id="pmx-metabologram-svg" viewBox="0 0 500 500" role="img" aria-label={tr('Metabologramme circulaire', 'Circular metabologram')}>
          <rect width="500" height="500" fill="white"/>
          <circle cx="250" cy="250" r="194" fill="none" stroke="#cfd6dc" stroke-width="1"/>
          <line x1="250" y1="54" x2="250" y2="446" stroke="#b9c1c8" stroke-width="1.5"/>

          {#each metaboliteSegments as segment}
            {@const stroke = qStroke(segment.qValue)}
            <path d={segment.path} fill={diverging(segment.effect, metLimit)} stroke="#263238" stroke-width={stroke.width} stroke-dasharray={stroke.dash}>
              <title>{itemLabel(segment)} · {segment.feature} · {scaleLabel([segment])}={fmt(segment.effect)} · q={fmt(segment.qValue, 3)}</title>
            </path>
          {/each}
          {#each transcriptSegments as segment}
            {@const stroke = qStroke(segment.qValue)}
            <path d={segment.path} fill={diverging(segment.effect, rnaLimit)} stroke="#263238" stroke-width={stroke.width} stroke-dasharray={stroke.dash}>
              <title>{itemLabel(segment)} · {segment.feature} · {scaleLabel([segment])}={fmt(segment.effect)} · q={fmt(segment.qValue, 3)}</title>
            </path>
          {/each}

          <path d="M250 196 A54 54 0 0 0 250 304 L250 250 Z" fill={diverging(activeMetabologram?.meanMetabolomicLog2Fc, metLimit)} stroke="#263238"/>
          <path d="M250 196 A54 54 0 0 1 250 304 L250 250 Z" fill={diverging(activeMetabologram?.meanTranscriptomicLog2Fc, rnaLimit)} stroke="#263238"/>
          <text x="196" y="244" text-anchor="middle" class="center-label">MET</text>
          <text x="304" y="244" text-anchor="middle" class="center-label">RNA</text>
          <text x="196" y="263" text-anchor="middle" class="center-value">{fmt(activeMetabologram?.meanMetabolomicLog2Fc)}</text>
          <text x="304" y="263" text-anchor="middle" class="center-value">{fmt(activeMetabologram?.meanTranscriptomicLog2Fc)}</text>
          <text x="132" y="32" text-anchor="middle" class="half-label">{tr('Métabolites', 'Metabolites')}</text>
          <text x="368" y="32" text-anchor="middle" class="half-label">{tr('Transcrits', 'Transcripts')}</text>
        </svg>

        <div class="metabologram-keys" data-testid="multiomics-readable-legend">
          <p class="color-key"><i class="key-blue"></i>{tr('Baisse', 'Decrease')}
            <i class="key-white"></i>{tr('Proche de zéro', 'Near zero')}
            <i class="key-red"></i>{tr('Hausse', 'Increase')}
            <small>{tr('Couleur normalisée séparément par couche ; ne pas comparer directement les intensités RNA et MET.', 'Color scaled within each layer; do not compare RNA and MET color intensity.')}</small>
          </p>
          <div class="legend-section">
            <h4>{tr('Métabolites', 'Metabolites')} <small>({activeMetabologram?.metabolomics?.length || 0}) · {scaleLabel(activeMetabologram?.metabolomics)}</small></h4>
            {#if activeMetabologram?.metabolomics?.length}
              <ul class="legend-items">
                {#each activeMetabologram.metabolomics.slice(0, 15) as item}
                  <li class="legend-entry">
                    <i class="legend-swatch" style={'background:' + diverging(item.effect, metLimit)}></i>
                    <span class="legend-name"><strong>{itemLabel(item)}</strong>{#if itemLabel(item) !== item.feature}<small>{item.feature}</small>{/if}</span>
                    <b class="legend-effect">{item.effect > 0 ? '+' : ''}{fmt(item.effect)}</b>
                  </li>
                {/each}
              </ul>
              {#if activeMetabologram.metabolomics.length > 15}<p class="legend-more">{tr('Liste limitée aux 15 premiers éléments.', 'List limited to the first 15 features.')}</p>{/if}
            {:else}
              <p class="legend-empty">{tr('Aucune mesure reconnue dans cette sélection.', 'No matched measurements in this selection.')}</p>
            {/if}
          </div>
          <div class="legend-section">
            <h4>{tr('Gènes exprimés', 'Gene transcripts')} <small>({activeMetabologram?.transcriptomics?.length || 0}) · {scaleLabel(activeMetabologram?.transcriptomics)}</small></h4>
            {#if activeMetabologram?.transcriptomics?.length}
              <ul class="legend-items">
                {#each activeMetabologram.transcriptomics.slice(0, 15) as item}
                  <li class="legend-entry">
                    <i class="legend-swatch" style={'background:' + diverging(item.effect, rnaLimit)}></i>
                    <span class="legend-name"><strong>{itemLabel(item)}</strong>{#if itemLabel(item) !== item.feature}<small>{item.feature}</small>{/if}</span>
                    <b class="legend-effect">{item.effect > 0 ? '+' : ''}{fmt(item.effect)}</b>
                  </li>
                {/each}
              </ul>
              {#if activeMetabologram.transcriptomics.length > 15}<p class="legend-more">{tr('Liste limitée aux 15 premiers éléments.', 'List limited to the first 15 features.')}</p>{/if}
            {:else}
              <p class="legend-empty">{tr('Aucune mesure reconnue dans cette sélection.', 'No matched measurements in this selection.')}</p>
            {/if}
          </div>
        </div>
      </div>
      {#if pathwayCoverage}
        <p class="coverage-note" data-testid="multiomics-pathway-coverage">
          {language === 'en' ? pathwayCoverage.labelEn : pathwayCoverage.labelFr} ·
          {pathwayCoverage.coverageMetabolites} {tr('métabolite(s)', 'metabolite(s)')} ·
          {pathwayCoverage.coverageTranscripts} {tr('transcrit(s) enzymatique(s)', 'enzyme transcript(s)')}
        </p>
      {/if}
      <p class="method-note">{activeMetabologram.method} {tr('La présence sur cette carte ne signifie pas un enrichissement significatif de la voie.', 'Being on this map does not imply significant pathway enrichment.')}</p>
    {:else}
      <p class="empty-note">{tr(
        'Aucun effet exploitable pour cette sélection. La voie peut comporter des molécules non mesurées, non identifiées ou sans estimation valide.',
        'No usable fitted effects for this selection. Features can be missing, unidentified or have no valid effect estimate.'
      )}</p>
    {/if}
  </article>

  {/if}

  {#if activeFigure === 'map'}
  <FocusedMetabolicNetwork graph={visuals?.focusedMetabolicNetwork} {language} />
  <details class="legacy-map" data-testid="multiomics-central-carbon-details">
    <summary>{tr('Voir aussi le schéma simplifié du métabolisme central', 'Also view the original central metabolism schematic')}</summary>
  <article class="figure-card" data-testid="multiomics-central-carbon-map">
    <div class="figure-title">
      <div>
        <span class="figure-number">03</span>
        <div>
          <h3>{tr('Carte du métabolisme central', 'Central carbon metabolism map')}</h3>
          <p>{tr(
            'Ovales = métabolites ; carrés = gènes enzymatiques. Gris = donnée absente, non reconnue ou inexploitable. Les flèches ne sont pas des flux mesurés.',
            'Ovals = metabolites; squares = enzyme transcripts. Gray = missing, unmatched or unusable data. Arrows are not measured fluxes.'
          )}</p>
        </div>
      </div>
      <button type="button" on:click={() => downloadSvg('pmx-central-carbon-svg', 'central_carbon_multiomics.svg')}>SVG</button>
    </div>

    {#if central}
      <div class="map-summary">
        <span><b>{central.measuredMetabolites}</b> {tr('métabolites reconnus sur la carte', 'metabolites matched on the map')}</span>
        <span><b>{central.measuredTranscripts}</b> {tr('gènes reconnus sur la carte', 'genes matched on the map')}</span>
      </div>

      {#if central.offMapMetabolites?.length || central.offMapTranscripts?.length}
        <div class="off-map" data-testid="multiomics-off-map">
          <strong>{tr('D’autres molécules sont mesurées, mais hors de cette carte.', 'Other molecules are measured, but outside this map.')}</strong>
          <p>{tr('Cette figure montre uniquement le métabolisme central. Une molécule absente de la carte reste disponible dans les autres graphiques.', 'This figure shows only central metabolism. Features outside the map remain in the other charts.')}</p>
          {#if central.offMapMetabolites?.length}
            <div><b>{tr('Métabolites hors carte : ', 'Metabolites outside map: ')}</b>{central.offMapMetabolites.map((item) => itemLabel(item)).join(', ')}</div>
          {/if}
          {#if central.offMapTranscripts?.length}
            <details><summary>{tr('Autres gènes mesurés', 'Other measured genes')}</summary>
              <p>{central.offMapTranscripts.map((item) => itemLabel(item)).join(', ')}</p>
            </details>
          {/if}
        </div>
      {/if}
      <div class="svg-scroll">
        <svg id="pmx-central-carbon-svg" viewBox="0 0 920 610" role="img" aria-label={tr('Carte du métabolisme central annotée par les effets sur leur échelle d’origine', 'Central metabolism map annotated with effects on their original scale')}>
          <defs>
            <marker id="arrow-central" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto">
              <path d="M0,0 L0,6 L8,3 z" fill="#98a3aa"/>
            </marker>
          </defs>
          <rect width="920" height="610" fill="white"/>
          <text x="30" y="30" class="svg-heading">{tr('Métabolisme central — effets différentiels', 'Central metabolism — differential effects')}</text>

          {#each central.edges as edge}
            {@const from = findCentralNode(edge[0])}
            {@const to = findCentralNode(edge[1])}
            {#if from && to}
              <line x1={from.x} y1={from.y} x2={to.x} y2={to.y} stroke="#aab3b8" stroke-width="1.4" marker-end="url(#arrow-central)"/>
            {/if}
          {/each}

          <ellipse cx="620" cy="375" rx="180" ry="150" fill="none" stroke="#d3dadd" stroke-dasharray="7 5"/>
          <text x="605" y="232" class="pathway-label">TCA</text>
          <text x="90" y="125" class="pathway-label">PPP</text>
          <text x="360" y="105" class="pathway-label">{tr('Glycolyse', 'Glycolysis')}</text>

          {#each central.metabolites as node}
            {@const stroke = qStroke(node.measurement?.qValue)}
            <g>
              <ellipse
                cx={node.x}
                cy={node.y}
                rx="36"
                ry="20"
                fill={measurementColor(node.measurement, mapMetLimit)}
                stroke="#253238"
                stroke-width={stroke.width}
                stroke-dasharray={stroke.dash}
              >
                <title>{measurementTitle(node.label, node.measurement)}</title>
              </ellipse>
              <text x={node.x} y={node.y + 4} text-anchor="middle" class="metabolite-label">{node.label}</text>
            </g>
          {/each}

          {#each central.enzymes as enzyme}
            {@const stroke = qStroke(enzyme.measurement?.qValue)}
            <g>
              <rect
                x={enzyme.x - 17}
                y={enzyme.y - 11}
                width="34"
                height="22"
                rx="3"
                fill={measurementColor(enzyme.measurement, mapRnaLimit)}
                stroke="#253238"
                stroke-width={stroke.width}
                stroke-dasharray={stroke.dash}
              >
                <title>{measurementTitle(enzyme.label, enzyme.measurement)}</title>
              </rect>
              <text x={enzyme.x} y={enzyme.y + 4} text-anchor="middle" class="enzyme-label">{enzyme.label}</text>
            </g>
          {/each}

          <g transform="translate(38 565)">
            <rect x="0" y="0" width="42" height="15" fill={diverging(-2)}/><text x="48" y="12" class="legend-label">{tr('baisse', 'decrease')}</text>
            <rect x="135" y="0" width="42" height="15" fill={diverging(0)}/><text x="183" y="12" class="legend-label">0</text>
            <rect x="220" y="0" width="42" height="15" fill={diverging(2)}/><text x="268" y="12" class="legend-label">{tr('hausse', 'increase')}</text>
            <rect x="375" y="0" width="42" height="15" fill="rgb(220,220,220)"/><text x="423" y="12" class="legend-label">{tr('aucune valeur exploitable', 'no usable value')}</text>
          </g>
        </svg>
      </div>
      <p class="method-note">{central.method} {central.scope} {tr(
        'Contrairement à la figure 4D de Guyon et al., cette carte ne montre aucun marquage isotopique 13C. Les niveaux de métabolites et les transcrits ne suffisent pas pour déduire un flux.',
        'Unlike Figure 4D by Guyon et al., this map contains no 13C isotope tracing. Metabolite abundance and transcripts alone do not establish metabolic flux.'
      )}</p>
    {/if}
  </article>
  </details>

  {/if}

  <details class="interpretation-boundary">
    <summary>{tr('Limites et méthode des figures', 'Figure methods and limitations')}</summary>

    {#each visuals.methodologicalBoundary || [] as note}
      <span>• {note}</span>
    {/each}
  </details>
</section>
{/if}

<style>
  /* Plots are intentionally presented on a light background even when PMx
     Explain is in dark mode. Reset the full palette in this subtree, not only
     inherited color: global h2/h3/code styles use theme CSS variables. */
  .visual-panel {
    --text-primary:#182a34;
    --text-secondary:#3b5060;
    --text-muted:#526473;
    --bg-primary:#ffffff;
    --bg-secondary:#eff3f5;
    --bg-tertiary:#ffffff;
    --border-subtle:#dbe3e8;
    --border-strong:#bdcbd3;
    --accent-pk:#954021;
    --accent-ai:#15536a;
    --focus-ring:#954021;
    margin-top:24px; display:grid; gap:18px; padding:14px; border-radius:12px;
    background:#fff; color:#182a34; color-scheme:light;
  }
  .legacy-map { border:1px solid #bdcbd3; background:#fff; color:#182a34; padding:12px 14px; border-radius:12px; }
  .legacy-map > summary { cursor:pointer; font-weight:700; }
  .legacy-map > article { margin-top:12px; }
  .figure-switcher { display:flex; flex-wrap:wrap; gap:8px; padding:5px; border:1px solid var(--border,#d9e0e3); border-radius:13px; width:max-content; max-width:100%; }
  .figure-switcher button { padding:10px 16px; border:0; border-radius:9px; background:#fff; color:#304755; font-size:.88rem; font-weight:650; }
  .figure-switcher button.active { background:#12627b; color:#fff; }
  .figure-switcher button:focus-visible { outline:3px solid #e5b75b; outline-offset:2px; }
  .pathway-select { display:flex; align-items:center; flex-wrap:wrap; gap:7px; font-size:.82rem; }
  .pathway-select span { color:#304755; }
  .coverage-note { color:#243947; font-size:.83rem; margin:8px 0; font-weight:700; }
  .interpretation-boundary summary { cursor:pointer; font-weight:700; }
  .interpretation-boundary[open] { gap:6px; }
  .visual-head { display:flex; justify-content:space-between; gap:20px; align-items:flex-end; }
  .visual-head h2 { margin:.2rem 0 .35rem; font-size:clamp(1.3rem,2vw,1.8rem); }
  .visual-head p { margin:0; max-width:850px; color:#3b5060; }
  .eyebrow { margin:0; font-size:.76rem; font-weight:800; letter-spacing:.11em; text-transform:uppercase; color:#15536a; }
  /* Static figures use white plotting surfaces. Make ink independent of dark
     parent themes, including cards, controls, footnotes and legends. */
  .figure-card { border:1px solid #c3d0d8; border-radius:18px; background:#fff; color:#182a34; color-scheme:light; padding:18px; overflow:hidden; }
  .figure-card :global(p), .figure-card :global(summary), .figure-card :global(span) { color:inherit; }
  .figure-title { display:flex; justify-content:space-between; gap:18px; align-items:flex-start; margin-bottom:14px; }
  .figure-title strong { color:#182a34 !important; }
  .figure-title > div:first-child { display:flex; gap:12px; align-items:flex-start; }
  .figure-title h3 { margin:0 0 5px; font-size:1.05rem; color:#182a34 !important; }
  .figure-title p { margin:0; max-width:780px; color:#42535f; font-size:.9rem; line-height:1.45; }
  .figure-number { display:grid; place-items:center; min-width:34px; height:34px; border-radius:50%; border:1px solid var(--border,#cbd5d9); font:700 .72rem var(--font-mono,monospace); }
  .figure-actions { display:flex; gap:8px; align-items:center; }
  .figure-card button, .figure-card select { border:1px solid #a8b9c4; border-radius:9px; background:#fff; padding:7px 10px; font:inherit; color:#172b36; }
  button { cursor:pointer; font-weight:700; }
  .svg-scroll { overflow:auto; border:1px solid #e4e8ea; border-radius:12px; background:white; }
  svg { display:block; width:100%; min-width:680px; height:auto; }
  .svg-heading { font-size:14px; font-weight:800; fill:#263238; }
  .svg-note { font-size:11px; fill:#344b59; }
  .sample-label { font-size:10px; fill:#253c4a; }
  .condition-mark { font-size:8px; font-weight:800; fill:#263238; }
  .feature-label { font-size:10px; fill:#203441; }
  .method-note, .empty-note { margin:10px 0 0; color:#3b5060; font-size:.82rem; line-height:1.45; }
  .metabologram-layout { display:grid; grid-template-columns:minmax(300px,1fr) minmax(300px,.95fr); gap:20px; align-items:start; }
  .metabologram-layout > svg { min-width:0; max-width:520px; margin:auto; }
  .center-label { font-size:11px; font-weight:800; fill:#263238; }
  .center-value { font-size:10px; fill:#263238; }
  .half-label { font-size:12px; font-weight:800; fill:#263238; }
  /* The plot has a white SVG background, whereas the rest of the portal can
     use a dark theme. Explicit ink + surface are required on every legend row:
     inheriting light text on a white legend caused unreadable labels. */
  .metabologram-keys { display:flex; flex-direction:column; gap:13px; min-width:0; font-size:.85rem;
    padding:12px; border:1px solid #cbd5dc; border-radius:12px;
    color:#1b2b34; background:#fff; color-scheme:light; }
  .color-key { display:flex; flex-wrap:wrap; align-items:center; gap:7px; margin:0; font-size:.83rem; color:#1b2b34; }
  .color-key i { display:inline-block; width:13px; height:13px; border:1px solid #829099; border-radius:3px; margin-left:8px; }
  .color-key .key-blue { background:#3567b2; margin-left:0; }
  .color-key .key-white { background:#f7f7f7; border:2px solid #566673; }
  .color-key .key-red { background:#cc3e3e; }
  .color-key small { display:block; width:100%; color:#415360; margin-top:5px; }
  .legend-section { border:1px solid #cbd5dc; border-radius:12px; min-width:0; overflow:hidden; color:#1b2b34; background:#fff; }
  .legend-section h4 { margin:0; color:#182731; font-size:.9rem; padding:11px 12px; background:#eaf0f3; border-bottom:1px solid #cbd5dc; }
  .legend-section h4 small { font-weight:400; color:#394e5b; }
  .legend-items { list-style:none; padding:4px 10px; margin:0; display:grid; color:#1b2b34; background:#fff; }
  .legend-entry { display:grid; grid-template-columns:14px minmax(0,1fr) 70px; align-items:center; gap:10px; padding:9px 3px; border-bottom:1px solid #dbe3e8; min-width:0; color:#1b2b34; background:#fff; }
  .legend-entry:last-child { border-bottom:0; }
  .legend-swatch { width:13px; height:13px; border-radius:3px; border:1px solid #839098; }
  .legend-name { display:flex; flex-direction:column; gap:1px; min-width:0; overflow-wrap:anywhere; line-height:1.3; }
  .legend-name strong { color:#182731; font-size:.88rem; font-weight:700; }
  .legend-name small { font-size:.73rem; color:#405463; font-family:var(--font-mono,monospace); }
  .legend-effect { color:#182731; text-align:right; font-family:var(--font-mono,monospace); font-variant-numeric:tabular-nums; font-size:.87rem; font-weight:750; }
  .legend-empty, .legend-more { display:block; padding:6px 12px; color:#405463; font-size:.82rem; }
  .off-map { margin:12px 0; padding:12px 14px; border:1px solid #c3d0d8; border-radius:10px; font-size:.85rem; line-height:1.5; color:#182a34; background:#fff; }
  .off-map strong { font-size:.9rem; }
  .off-map p { margin:4px 0 8px; color:#3f5260; }
  .off-map details { margin-top:6px; }
  .off-map summary { cursor:pointer; }
  .map-summary { display:flex; gap:16px; flex-wrap:wrap; margin-bottom:10px; font-size:.82rem; }
  .pathway-label { font-size:16px; font-weight:800; letter-spacing:.08em; fill:#364c59; }
  .metabolite-label { font-size:10px; font-weight:800; fill:#172126; pointer-events:none; }
  .enzyme-label { font-size:8px; font-weight:800; fill:#172126; pointer-events:none; }
  .legend-label { font-size:9px; fill:#48555b; }
  .interpretation-boundary { border-left:4px solid #176c83; padding:12px 14px; background:#f3f8f9; color:#203845; display:grid; gap:4px; font-size:.85rem; }
  .interpretation-boundary summary { margin-bottom:4px; }
  @media (max-width: 820px) {
    .figure-title, .visual-head { flex-direction:column; align-items:stretch; }
    .figure-actions { justify-content:space-between; flex-wrap:wrap; }
    .figure-switcher { width:100%; }
    .figure-switcher button { flex:1 1 auto; }
    .metabologram-layout { grid-template-columns:1fr; }
    svg { min-width:620px; }
    .metabologram-layout > svg { min-width:0; width:100%; }
  }
</style>
