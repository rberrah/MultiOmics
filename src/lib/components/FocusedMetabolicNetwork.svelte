<script>
  // @ts-nocheck
  import { focusMetabolicRegion } from '$lib/multiomics/metabolic-network.js';
  export let graph = null;
  export let language = 'fr';

  let selectedRegion = 'auto';
  let neighborhood = '1';
  const tr = (fr,en) => language === 'en' ? en : fr;
  const fmt = (n) => Number.isFinite(n) ? n.toFixed(2) : '—';
  function colorLimit(card, kind) {
    const vals = card.nodes.filter((node)=>node.kind===kind && Number.isFinite(node.measurement?.effect))
      .map((node)=>node.measurement);
    if (!vals.length || vals.every((row)=>row.effectScale==='log2')) return 3;
    return Math.max(1e-9,...vals.map((row)=>Math.abs(row.effect)));
  }
  function color(measure, limit=3) {
    if (measure?.status === 'duplicate_features') return '#f5dba8';
    if (!Number.isFinite(measure?.effect)) return '#e5e9eb';
    const value = Math.max(-1,Math.min(1,measure.effect / limit));
    const white=[248,248,248],target=value<0?[51,104,178]:[204,61,65];
    return 'rgb(' + white.map((c,i)=>Math.round(c + Math.abs(value)*(target[i]-c))).join(',') + ')';
  }
  function displayName(node) {
    if(node.kind==='gene') return node.name;
    return language==='en' ? node.labelEn : node.labelFr;
  }
  function shortName(node) {
    const name=displayName(node);
    return name.length>21 ? name.slice(0,20)+'…' : name;
  }
  function positionCard(card) {
    const metab=card.nodes.filter((n)=>n.kind==='metabolite');
    const genes=card.nodes.filter((n)=>n.kind==='gene');
    const columns=4, x0=110, dx=200, y0=65, dy=106;
    const positions={};
    for(const [index,node] of metab.entries())positions[node.id]={x:x0+(index%columns)*dx,y:y0+Math.floor(index/columns)*dy};
    const geney=y0+Math.ceil(metab.length/columns)*dy+10;
    for(const [index,node] of genes.entries())positions[node.id]={x:x0+(index%columns)*dx,y:geney+Math.floor(index/columns)*92};
    return {...card,positions,height:Math.max(160,geney+Math.ceil(genes.length/columns)*92+40),
      metLimit:colorLimit(card,'metabolite'),geneLimit:colorLimit(card,'gene')};
  }
  $: views = focusMetabolicRegion(graph,selectedRegion,neighborhood)
    .map(positionCard).filter((v)=>v.nodes.length);
  $: warnings = graph?.audit?.filter((a)=>!a.id || a.status==='name_only') || [];
  $: unknown = graph?.audit?.filter((a)=>!a.id) || [];
  $: exactMatches = graph?.audit?.filter((a)=>['verified_id','verified_crossref'].includes(a.status))?.length || 0;
  $: nameMatches = graph?.audit?.filter((a)=>a.status==='name_only')?.length || 0;
  $: declaredMatches = graph?.audit?.filter((a)=>a.status==='user_declared')?.length || 0;
</script>

<section class="network" data-testid="multiomics-focused-network">
  <header>
    <h3>{tr('Réseau métabolique centré sur vos données','Metabolic network focused on your data')}</h3>
    <p>{tr('Seules les voies liées à vos mesures apparaissent. Les points gris sont des voisins biologiques, pas des mesures inventées.',
      'Only pathways connected to your measurements are shown. Gray nodes are biological context, not invented measurements.')}</p>
  </header>
  <div class="network-summary">
    <span><strong>{graph?.coverage?.metabolites || 0}</strong> {tr('métabolites dans le dictionnaire','metabolites in catalog')}</span>
    <span><strong>{graph?.coverage?.genes || 0}</strong> {tr('gènes dans le réseau','genes in network')}</span>
    <span><strong>{graph?.coverage?.regions || 0}</strong> {tr('régions disponibles','available regions')}</span>
  </div>
  <div class="network-controls">
    <label>{tr('Région affichée','Shown region')}
      <select bind:value={selectedRegion} data-testid="multiomics-network-region">
        <option value="auto">{tr('Automatique · régions contenant des mesures','Automatic · regions with measurements')}</option>
        {#each graph?.regions || [] as region}
          <option value={region.id}>{language==='en' ? region.labelEn : region.labelFr} ({region.measuredCount})</option>
        {/each}
      </select>
    </label>
    <label>{tr('Autour des molécules mesurées','Around measured molecules')}
      <select bind:value={neighborhood} data-testid="multiomics-network-hops">
        <option value="0">{tr('Mesures seules','Measurements only')}</option>
        <option value="1">{tr('Voisins directs (conseillé)','Direct neighbors (recommended)')}</option>
        <option value="2">{tr('Voisinage élargi','Wider neighborhood')}</option>
      </select>
    </label>
  </div>
  <div class="network-key" aria-label={tr('Légende du réseau','Network legend')}>
    <span><i style="background:#3368b2"></i>{tr('Baisse','Decrease')}</span>
    <span><i style="background:#cc3d41"></i>{tr('Hausse','Increase')}</span>
    <span><i style="background:#e5e9eb"></i>{tr('Non mesuré / effet indisponible','Not measured / effect unavailable')}</span>
    <span><i style="background:#f5dba8"></i>{tr('Plusieurs valeurs non fusionnées','Multiple values not merged')}</span>
  </div>
  {#if !views.length}
    <p class="network-empty">{tr('Aucune correspondance exploitable dans cette région. Choisissez une autre voie ou vérifiez le dictionnaire ci-dessous.',
      'No matching measurement for this region. Select another pathway or review the dictionary below.')}</p>
  {/if}
  {#each views as card (card.id)}
    <article class="network-card" data-testid={'multiomics-region-'+card.id}>
      <div class="network-card-title">
        <strong>{language==='en' ? card.labelEn : card.labelFr}</strong>
        <span>{card.measuredCount} {tr('mesure(s) colorée(s)','colored measurement(s)')} · {card.nodes.length} {tr('éléments visibles','visible nodes')}</span>
      </div>
      <div class="network-scroll">
        <svg viewBox={'0 0 820 '+card.height} role="img" aria-label={tr('Voisinage métabolique : ','Metabolic neighborhood: ')+(language==='en'?card.labelEn:card.labelFr)}>
          <rect x="0" y="0" width="820" height={card.height} fill="#fff"/>
          {#each card.links as edge}
            {@const a=card.positions[edge.from]}
            {@const b=card.positions[edge.to]}
            {#if a && b}
              <line x1={a.x} y1={a.y} x2={b.x} y2={b.y}
                stroke={edge.relation==='gene_context'?'#899fab':'#b3bfc7'}
                stroke-width={edge.relation==='gene_context'?1.5:2.3}
                stroke-dasharray={edge.relation==='gene_context'?'4 4':undefined}/>
            {/if}
          {/each}
          {#each card.nodes as node}
            {@const p=card.positions[node.id]}
            {@const m=node.measurement}
            <g>
              {#if node.kind==='gene'}
                <rect x={p.x-65} y={p.y-20} width="130" height="40" rx="7"
                  stroke={m?.matchStatus==='name_only'?'#ab802f':'#4c6876'}
                  stroke-width={Number.isFinite(m?.effect)?2.5:1}
                  stroke-dasharray={m?.matchStatus==='name_only'?'4 2':undefined}
                  fill={color(m,card.geneLimit)}/>
                <text x={p.x} y={p.y+5} text-anchor="middle" class="gene-name">{shortName(node)}</text>
              {:else}
                <ellipse cx={p.x} cy={p.y} rx="79" ry="27"
                  stroke={m?.matchStatus==='name_only'?'#ab802f':'#4c6876'}
                  stroke-width={Number.isFinite(m?.effect)?2.5:1}
                  stroke-dasharray={m?.matchStatus==='name_only'?'4 2':undefined}
                  fill={color(m,card.metLimit)}/>
                <text x={p.x} y={p.y+5} text-anchor="middle" class="metabolite-name">{shortName(node)}</text>
              {/if}
              <title>{displayName(node)}{node.chebi ? ' · '+node.chebi : ''}{m?.features ? ' · '+m.features.join(', ') : ''}{m?.feature ? ' · '+m.feature : ''}{m?.matchStatus ? ' · '+m.matchStatus : ''} · {m?.effectScale==='log2'?'log2FC':tr('différence (échelle fournie)','difference (supplied scale)')}={fmt(m?.effect)} · q={fmt(m?.qValue)}</title>
            </g>
          {/each}
        </svg>
      </div>
    </article>
  {/each}
  <details class="identity-audit" data-testid="multiomics-network-dictionary">
    <summary>{tr('Vérifier les correspondances des identifiants','Inspect identifier matching')} — {exactMatches} {tr('ID vérifiés','verified IDs')}, {nameMatches} {tr('noms seuls','name-only')}, {declaredMatches} {tr('annotations déclarées','declared annotations')}, {unknown.length} {tr('non reconnus','unmatched')}</summary>
    <p>{tr('Le dictionnaire privilégie les identifiants ChEBI exacts. Un nom seul est une correspondance provisoire, signalée par une bordure pointillée. Les identifiants inconnus ne sont jamais rattachés à une molécule par simple ressemblance.',
      'The registry prioritizes exact ChEBI IDs. Name-only matches are provisional and have dashed borders. Unknown identifiers are never assigned by fuzzy matching.')}</p>
    {#if graph?.audit?.length}
      <div class="audit-table-wrap">
        <table><thead><tr><th>{tr('Entrée','Input')}</th><th>{tr('Correspondance','Matched as')}</th><th>{tr('Statut','Status')}</th></tr></thead>
          <tbody>
            {#each graph.audit as item}
              <tr><td><code>{item.input}</code></td>
                <td>{item.label || '—'}
                  {#if item.chebi}<a href={'https://www.ebi.ac.uk/chebi/'+item.chebi} target="_blank" rel="noopener noreferrer">{item.chebi} ↗</a>{/if}
                </td>
                <td>{item.status==='verified_id'?tr('Identifiant ChEBI exact','Exact ChEBI ID')
                  :item.status==='verified_crossref'?tr('Identifiant croisé vérifié','Verified cross-reference')
                  :item.status==='user_declared'?tr('Annotation importée · identité non vérifiée','Imported annotation · identity unverified')
                  :item.status==='declared_not_in_network'?tr('Annotation hors catalogue','Annotation outside catalog')
                  :item.status==='name_only'?tr('Nom seul · à confirmer','Name only · confirm')
                  :item.status==='ambiguous_name'?tr('Nom ambigu','Ambiguous name')
                  :item.status==='unmapped_identifier'?tr('Identifiant hors dictionnaire','Identifier outside dictionary')
                  :tr('Non reconnu','Not recognized')}</td></tr>
            {/each}
          </tbody>
        </table>
      </div>
    {/if}
  </details>
  <p class="network-limit">{tr(
    'Réseau de voisinage pédagogique, construit à partir de voies simplifiées. Les liens indiquent un contexte métabolique, pas des réactions complètes, des flux, une activité enzymatique ou une causalité. Une mesure d’expression génique ne prouve pas un flux métabolique.',
    'Educational neighborhood graph derived from simplified pathways. Links indicate metabolic context, not complete reactions, flux, enzyme activity or causation. Gene expression does not establish metabolic flux.'
  )}</p>
</section>

<style>
  .network { display:grid; gap:14px; margin:0 0 18px; color:#182b35; background:#fff; color-scheme:light; padding:14px; border:1px solid #c3d0d8; border-radius:12px; }
  .network header h3 { font-size:1.12rem; margin:0 0 5px; color:#182b35 !important; }
  .network header p { margin:0; color:#3e5260; font-size:.88rem; line-height:1.45; }
  .network-summary { display:flex; flex-wrap:wrap; gap:9px; }
  .network-summary span { padding:7px 10px; border:1px solid #cbd5dc; border-radius:8px; font-size:.79rem; color:#182b35; background:#fff; }
  .network-controls { display:flex; flex-wrap:wrap; gap:12px; }
  .network-controls label { display:flex; flex:1 1 230px; flex-direction:column; gap:5px; font-size:.85rem; font-weight:600; color:#182b35; }
  .network-controls select { width:100%; min-width:0; border:1px solid #aebcc5; border-radius:8px; padding:9px; background:#fff; color:#182b35; }
  .network-key { display:flex; flex-wrap:wrap; gap:8px 14px; color:#3e5260; font-size:.8rem; }
  .network-key span { display:inline-flex; align-items:center; gap:5px; }
  .network-key i { width:12px; height:12px; border-radius:3px; border:1px solid #8a9ca6; }
  .network-card { border:1px solid var(--border,#dce4e8); border-radius:12px; overflow:hidden; }
  .network-card-title { display:flex; align-items:baseline; justify-content:space-between; flex-wrap:wrap; gap:6px; padding:12px 14px; background:#f2f6f8; }
  .network-card-title strong { font-size:.95rem; color:#203440; }
  .network-card-title span { font-size:.78rem; color:#5d6c74; }
  .network-scroll { overflow:auto; background:#fff; }
  .network-scroll svg { display:block; width:100%; min-width:660px; height:auto; }
  .network-scroll text { font-family:inherit; pointer-events:none; fill:#1c2d35; font-weight:650; }
  .metabolite-name { font-size:12px; }
  .gene-name { font-size:13px; }
  .network-empty { padding:18px; border:1px dashed #c7d5db; border-radius:10px; }
  .identity-audit { border:1px solid #c8d4dc; border-radius:10px; padding:12px 14px; color:#182b35; background:#fff; }
  .identity-audit summary { cursor:pointer; font-weight:700; font-size:.87rem; }
  .identity-audit p { font-size:.83rem; line-height:1.5; }
  .audit-table-wrap { max-height:320px; overflow:auto; }
  table { border-collapse:collapse; width:100%; min-width:520px; font-size:.8rem; }
  th,td { text-align:left; padding:8px; border-bottom:1px solid #dde5e9; }
  td a { margin-left:8px; white-space:nowrap; font-size:.76rem; }
  code { font-family:var(--font-mono,monospace); font-size:.76rem; }
  .network-limit { font-size:.8rem; line-height:1.55; color:#3e5260; margin:0; }
</style>
