// @ts-nocheck
// DISPLAY ONLY. Human-readable, deliberately strict metabolite identity registry
// and expandable biochemical context. This module does not infer reactions,
// enzyme activities, isotope flux, pathway enrichment, or causality.
// ChEBI primary identifiers are manually reviewed; unknown/ambiguous IDs are
// explicitly reported, never fuzzy-matched. See metabolic-network.md.
const metaboliteRecords = [
  // id ; French label ; English label ; exact primary ChEBI ; exact-name aliases
  ['tryptophan','L-tryptophane','L-tryptophan','CHEBI:16828','L-tryptophan|tryptophan'],
  ['formylkyn','N-formyl-L-kynurénine','N-formyl-L-kynurenine','CHEBI:30249','N-formyl-L-kynurenine'],
  ['kynurenine','L-kynurénine','L-kynurenine','CHEBI:16946','kynurenine|L-kynurenine'],
  ['hydroxykyn','3-hydroxy-L-kynurénine','3-hydroxy-L-kynurenine','CHEBI:17380','3-hydroxy-L-kynurenine'],
  ['quinolinate','Acide quinolinique','Quinolinic acid','CHEBI:16675','quinolinic acid'],
  ['kynurenate','Acide kynurénique','Kynurenic acid','','kynurenic acid'],
  ['anthranilate','Acide anthranilique','Anthranilic acid','','anthranilic acid'],
  ['hydroxyanthranilate','Acide 3-hydroxyanthranilique','3-hydroxyanthranilic acid','','3-hydroxyanthranilic acid'],
  ['serotonin','Sérotonine','Serotonin','','serotonin|5-hydroxytryptamine'],
  ['hydroxytrp','5-hydroxytryptophane','5-hydroxytryptophan','','5-hydroxytryptophan'],
  ['melatonin','Mélatonine','Melatonin','','melatonin'],
  ['niacin','Acide nicotinique','Nicotinic acid','','nicotinic acid|niacin'],
  ['nad','NAD+','NAD+','','nad+|nicotinamide adenine dinucleotide'],
  ['nadp','NADP+','NADP+','','nadp+'],
  ['nmna','Acide nicotinique mononucléotide','Nicotinic acid mononucleotide','','nicotinic acid mononucleotide'],
  ['glutathione','Glutathion réduit','Reduced glutathione','','reduced glutathione|gsh'],
  ['gssg','Glutathion oxydé','Oxidized glutathione','','oxidized glutathione|gssg'],
  ['cysteine','L-cystéine','L-cysteine','','l-cysteine|cysteine'],
  ['cystathionine','Cystathionine','Cystathionine','','cystathionine'],
  ['methionine','L-méthionine','L-methionine','','l-methionine|methionine'],
  ['sam','S-adénosylméthionine','S-adenosylmethionine','','s-adenosylmethionine|sam'],
  ['sah','S-adénosylhomocystéine','S-adenosylhomocysteine','','s-adenosylhomocysteine|sah'],
  ['homocysteine','Homocystéine','Homocysteine','','homocysteine'],
  ['ornithine','L-ornithine','L-ornithine','','l-ornithine|ornithine'],
  ['citrulline','L-citrulline','L-citrulline','','l-citrulline|citrulline'],
  ['argininosuccinate','Argininosuccinate','Argininosuccinate','','argininosuccinate'],
  ['arginine','L-arginine','L-arginine','','l-arginine|arginine'],
  ['urea','Urée','Urea','','urea'],
  ['carbamoylp','Carbamoyl-phosphate','Carbamoyl phosphate','','carbamoyl phosphate'],
  ['ammonium','Ammonium','Ammonium','','ammonium'],
  ['imp','IMP','IMP','','inosine monophosphate|imp'],
  ['amp','AMP','AMP','','adenosine monophosphate|amp'],
  ['gmp','GMP','GMP','','guanosine monophosphate|gmp'],
  ['inosine','Inosine','Inosine','','inosine'],
  ['hypoxanthine','Hypoxanthine','Hypoxanthine','','hypoxanthine'],
  ['xanthine','Xanthine','Xanthine','','xanthine'],
  ['urate','Acide urique','Uric acid','','uric acid|urate'],
  ['atp','ATP','ATP','','adenosine triphosphate|atp'],
  ['adp','ADP','ADP','','adenosine diphosphate|adp'],
  ['adenosine','Adénosine','Adenosine','','adenosine'],
  ['palmitate','Palmitate','Palmitate','','palmitate|palmitic acid'],
  ['palmitoylcoa','Palmitoyl-CoA','Palmitoyl-CoA','','palmitoyl-coa'],
  ['palmitoylcarnitine','Palmitoylcarnitine','Palmitoylcarnitine','','palmitoylcarnitine'],
  ['acetylcarnitine','Acétylcarnitine','Acetylcarnitine','','acetylcarnitine'],
  ['carnitine','L-carnitine','L-carnitine','','l-carnitine|carnitine'],
  ['octanoylcoa','Octanoyl-CoA','Octanoyl-CoA','','octanoyl-coa'],
  ['hydroxybutyrate','3-hydroxybutyrate','3-hydroxybutyrate','','3-hydroxybutyrate|beta-hydroxybutyrate'],
  ['acetoacetate','Acétoacétate','Acetoacetate','','acetoacetate'],
  ['hmgcoa','HMG-CoA','HMG-CoA','','hmg-coa'],
  ['cholesterol','Cholestérol','Cholesterol','','cholesterol'],
  ['cholicacid','Acide cholique','Cholic acid','','cholic acid'],
  ['chenodeoxycholicacid','Acide chénodésoxycholique','Chenodeoxycholic acid','','chenodeoxycholic acid'],
  ['deoxycholicacid','Acide désoxycholique','Deoxycholic acid','','deoxycholic acid'],
  ['phosphatidylcholine','Phosphatidylcholine','Phosphatidylcholine','','phosphatidylcholine'],
  ['choline','Choline','Choline','CHEBI:15354','choline'],
  ['phosphocholine','Phosphocholine','Phosphocholine','CHEBI:18132','phosphocholine'],
  ['glycerol','Glycérol','Glycerol','','glycerol'],
  ['dihydroxyacetonep','DHAP','Dihydroxyacetone phosphate','','dihydroxyacetone phosphate|dhap'],
  ['serinehydroxymethyl','5,10-méthylène-THF','5,10-methylene-THF','','5,10-methylene-thf'],
  ['thf','Tétrahydrofolate','Tetrahydrofolate','','tetrahydrofolate|thf'],
  ['folate','Folate','Folate','','folate|folic acid'],
  ['formylthf','10-formyl-THF','10-formyl-THF','','10-formyl-thf'],
  ['dtmp','dTMP','dTMP','','thymidine monophosphate|dtmp'],
  ['dump','dUMP','dUMP','','deoxyuridine monophosphate|dump'],
  ['heme','Hème','Heme','','heme|haem'],
  ['biliverdin','Biliverdine','Biliverdin','','biliverdin'],
  ['bilirubin','Bilirubine','Bilirubin','','bilirubin'],
  ['protoporphyrin','Protoporphyrine IX','Protoporphyrin IX','','protoporphyrin ix'],
  ['succinylcoa','Succinyl-CoA','Succinyl-CoA','','succinyl-coa'],
  ['deltaala','Acide δ-aminolévulinique','5-aminolevulinic acid','','5-aminolevulinic acid|delta-ala'],
  ['glutamine','L-glutamine','L-glutamine','CHEBI:18050','l-glutamine|glutamine'],
  ['glutamate','L-glutamate (1−)','L-glutamate (1−)','CHEBI:29985','l-glutamate|glutamate'],
  ['aspartate','L-aspartate','L-aspartate','CHEBI:29993','l-aspartate|aspartate'],
  ['asparagine','L-asparagine','L-asparagine','CHEBI:17196','l-asparagine|asparagine'],
  ['lactate','Lactate (anion)','Lactate (anion)','CHEBI:24996','lactate'],
  ['succinate','Succinate (2−)','Succinate (2−)','CHEBI:30031','succinate|succinate(2-)'],
  ['gaba','GABA','GABA','CHEBI:30566','4-aminobutyrate|gamma-aminobutyric acid|gaba']
];
const regionSpecs = [
  // ID, French, English, ordered metabolites, related enzyme/transcript names
  ['glycolysis','Glycolyse','Glycolysis','glucose g6p f6p fbp g3p 3pg pep pyruvate lactate alanine','hk gpi pfk aldo gapdh eno pkm ldha ldhb'],
  ['tca','Cycle de Krebs','TCA cycle','pyruvate acetylcoa citrate akg succinylcoa succinate fumarate malate oaa 2hg','pdh cs idh ogdh suclg sdh fh mdh'],
  ['ppp','Pentoses phosphates','Pentose phosphate','g6p 6pg r5p prpp','g6pd pgd'],
  ['amino','Acides aminés','Amino acids','glutamine glutamate akg gaba serine glycine aspartate asparagine alanine 3pg','gls glud glul got phgdh shmt'],
  ['kynurenine','Tryptophane et kynurénine','Tryptophan and kynurenine','tryptophan formylkyn kynurenine hydroxykyn hydroxyanthranilate quinolinate kynurenate anthranilate nmna nad','IDO1 IDO2 TDO2 AFMID KMO KYNU AADAT HAAO QPRT'],
  ['serotonin','Sérotonine et mélatonine','Serotonin and melatonin','tryptophan hydroxytrp serotonin melatonin','TPH1 TPH2 DDC AANAT ASMT'],
  ['redox','Glutathion et stress oxydant','Glutathione and oxidative stress','glutamate cysteine glycine glutathione gssg nadp','GCLC GCLM GSS GPX1 GPX4 GSR G6PD'],
  ['urea','Cycle de l’urée','Urea cycle','glutamate ammonium carbamoylp citrulline ornithine argininosuccinate arginine urea aspartate','CPS1 OTC ASS1 ASL ARG1 NAGS'],
  ['purines','Purines et nucléotides','Purine and nucleotide metabolism','r5p prpp imp amp gmp inosine hypoxanthine xanthine urate adenosine atp adp','PPAT HPRT1 IMPDH1 IMPDH2 GMPS ADSS ADSL XDH ADA PNP'],
  ['lipids','Acides gras et β-oxydation','Fatty acid oxidation','palmitate palmitoylcoa carnitine palmitoylcarnitine octanoylcoa acetylcoa acetylcarnitine acetoacetate hydroxybutyrate','ACSL1 CPT1A CPT2 SLC25A20 ACADM HADHA HADHB ACAT1'],
  ['cholesterol','Cholestérol et acides biliaires','Cholesterol and bile acids','acetylcoa hmgcoa cholesterol cholicacid chenodeoxycholicacid deoxycholicacid','HMGCS1 HMGCR CYP7A1 CYP8B1 CYP27A1'],
  ['choline','Choline et membranes','Choline and membranes','choline phosphocholine phosphatidylcholine glycerol gly3p dihydroxyacetonep','CHKA PCYT1A CHPT1 PEMT GPD1'],
  ['onecarbon','Folates et métabolisme à un carbone','Folate and one-carbon metabolism','serine glycine folate thf serinehydroxymethyl formylthf dump dtmp methionine sam sah homocysteine cystathionine cysteine','SHMT1 SHMT2 MTHFD1 MTHFR TYMS MTR MAT1A AHCY CBS CTH'],
  ['heme','Métabolisme de l’hème','Heme metabolism','succinylcoa glycine deltaala protoporphyrin heme biliverdin bilirubin','ALAS1 ALAS2 FECH HMOX1 BLVRA']
];
// Curated, intentionally schematic connections; not stoichiometrically complete reactions.
const extraLinks = [
  'tryptophan>formylkyn','formylkyn>kynurenine','kynurenine>hydroxykyn',
  'hydroxykyn>hydroxyanthranilate','hydroxyanthranilate>quinolinate','quinolinate>nmna','nmna>nad',
  'kynurenine>kynurenate','kynurenine>anthranilate',
  'tryptophan>hydroxytrp','hydroxytrp>serotonin','serotonin>melatonin',
  'glutamate>glutathione','cysteine>glutathione','glycine>glutathione','glutathione>gssg',
  'methionine>sam','sam>sah','sah>homocysteine','homocysteine>cystathionine','cystathionine>cysteine',
  'glutamate>ammonium','ammonium>carbamoylp','carbamoylp>citrulline','ornithine>citrulline',
  'citrulline>argininosuccinate','aspartate>argininosuccinate','argininosuccinate>arginine',
  'arginine>ornithine','arginine>urea',
  'r5p>prpp','prpp>imp','imp>amp','imp>gmp','amp>adenosine','adenosine>inosine',
  'inosine>hypoxanthine','hypoxanthine>xanthine','xanthine>urate','amp>adp','adp>atp',
  'palmitate>palmitoylcoa','palmitoylcoa>palmitoylcarnitine','carnitine>palmitoylcarnitine',
  'palmitoylcarnitine>octanoylcoa','octanoylcoa>acetylcoa','acetylcoa>acetylcarnitine',
  'acetylcoa>acetoacetate','acetoacetate>hydroxybutyrate',
  'acetylcoa>hmgcoa','hmgcoa>cholesterol','cholesterol>cholicacid',
  'cholesterol>chenodeoxycholicacid','cholicacid>deoxycholicacid',
  'choline>phosphocholine','phosphocholine>phosphatidylcholine',
  'glycerol>gly3p','gly3p>dihydroxyacetonep',
  'serine>glycine','folate>thf','thf>serinehydroxymethyl','serine>serinehydroxymethyl',
  'serinehydroxymethyl>formylthf','dump>dtmp',
  'succinylcoa>deltaala','glycine>deltaala','deltaala>protoporphyrin',
  'protoporphyrin>heme','heme>biliverdin','biliverdin>bilirubin','glutamate>gaba',
  'akg>succinylcoa','succinylcoa>succinate'
];
const geneAnchor = {
  hk:'glucose',gpi:'g6p',pfk:'f6p',aldo:'fbp',gapdh:'g3p',
  eno:'3pg',pkm:'pep',ldha:'pyruvate',ldhb:'lactate',
  pdh:'pyruvate',cs:'acetylcoa',idh:'citrate',ogdh:'akg',
  suclg:'succinylcoa',sdh:'succinate',fh:'fumarate',mdh:'malate',
  got:'aspartate',gls:'glutamine',glud:'glutamate',glul:'glutamate',
  g6pd:'g6p',pgd:'6pg',phgdh:'3pg',shmt:'serine',
  ido1:'tryptophan',ido2:'tryptophan',tdo2:'tryptophan',afmid:'formylkyn',
  kmo:'kynurenine',kynu:'kynurenine',aadat:'kynurenine',haao:'hydroxyanthranilate',qprt:'quinolinate',
  tph1:'tryptophan',tph2:'tryptophan',ddc:'hydroxytrp',aanat:'serotonin',asmt:'serotonin',
  gclc:'glutamate',gclm:'glutamate',gss:'glutathione',gpx1:'glutathione',gpx4:'glutathione',gsr:'gssg',
  cps1:'ammonium',otc:'ornithine',ass1:'citrulline',asl:'argininosuccinate',arg1:'arginine',nags:'glutamate',
  ppat:'prpp',hprt1:'hypoxanthine',impdh1:'imp',impdh2:'imp',gmps:'gmp',adss:'imp',adsl:'imp',xdh:'xanthine',ada:'adenosine',pnp:'inosine',
  acsl1:'palmitate',cpt1a:'palmitoylcoa',cpt2:'palmitoylcarnitine',slc25a20:'palmitoylcarnitine',
  acadm:'octanoylcoa',hadha:'octanoylcoa',hadhb:'octanoylcoa',acat1:'acetylcoa',
  hmgcs1:'acetylcoa',hmgcr:'hmgcoa',cyp7a1:'cholesterol',cyp8b1:'cholicacid',cyp27a1:'cholesterol',
  chka:'choline',pcyt1a:'phosphocholine',chpt1:'phosphatidylcholine',pemt:'phosphatidylcholine',gpd1:'gly3p',
  shmt1:'serine',shmt2:'serine',mthfd1:'formylthf',mthfr:'serinehydroxymethyl',
  tyms:'dump',mtr:'homocysteine',mat1a:'methionine',ahcy:'sah',cbs:'homocysteine',cth:'cystathionine',
  alas1:'deltaala',alas2:'deltaala',fech:'protoporphyrin',hmox1:'heme',blvra:'biliverdin'
};
const key = (s) => String(s ?? '').trim().toLowerCase().normalize('NFKD')
  .replace(/[\u0300-\u036f]/g, '').replace(/[^a-z0-9]+/g, '');
const isIdentifier = (s) => /^(?:CHEBI:\d+|HMDB\d+|C\d{5}|CID:?\d+|[A-Z]{14}-[A-Z]{10}-[A-Z])$/i.test(String(s ?? '').trim());
const list = (s) => String(s || '').split(' ').filter(Boolean);
function toRegistry(central) {
  const chemicals = new Map();
  for (const n of central.metabolites || []) {
    chemicals.set(n.id, { id:n.id, labelFr:n.label, labelEn:n.label, chebi:null,
      aliases:n.aliases || [], kind:'metabolite' });
  }
  for (const [id,fr,en,chebi,aliases] of metaboliteRecords) {
    const old = chemicals.get(id);
    chemicals.set(id, { ...old, id, labelFr:fr, labelEn:en, chebi:chebi || old?.chebi || null,
      aliases:[...new Set([...aliases.split('|').filter(Boolean),fr,en])],
      kind:'metabolite' });
  }
  const genes = new Map((central.enzymes || []).map((e) => [e.id,
    { id:'gene:'+e.id, name:e.label, aliases:e.aliases || [e.label], kind:'gene' }]));
  const allGenes = new Set(regionSpecs.flatMap((r) => list(r[4])));
  for (const symbol of allGenes) {
    const id = symbol.toLowerCase();
    if (!genes.has(id)) genes.set(id,{ id:'gene:'+id, name:symbol, aliases:[symbol], kind:'gene' });
  }
  return { chemicals, genes };
}
function identityDictionary(chemicals) {
  const exact = new Map(), names = new Map();
  for (const mol of chemicals.values()) {
    if (mol.chebi) {
      const id = mol.chebi.toUpperCase();
      if (exact.has(id) && exact.get(id) !== mol.id) throw new Error('Duplicate ChEBI: '+id);
      exact.set(id,mol.id);
    }
    for (const alias of [mol.labelFr,mol.labelEn,...mol.aliases]) {
      if (isIdentifier(alias)) continue;
      const norm = key(alias);
      if (!norm) continue;
      if (!names.has(norm)) names.set(norm,new Set());
      names.get(norm).add(mol.id);
    }
  }
  return { exact,names };
}
function resolveIdentity(raw,dictionary,mappings,userAnnotations) {
  const original = String(raw || '').trim();
  const uppercase = original.toUpperCase();
  if (dictionary.exact.has(uppercase)) return { status:'verified_id', id:dictionary.exact.get(uppercase), provenance:'ChEBI' };
  const declared = userAnnotations?.get(original);
  if (declared) {
    if (dictionary.exact.has(declared)) return {
      status:'user_declared', id:dictionary.exact.get(declared), provenance:'unverified_user_crosswalk'
    };
    return { status:'declared_not_in_network', id:null, provenance:'unverified_user_crosswalk' };
  }
  if (isIdentifier(original)) {
    const mapped = mappings?.get(original);
    const status = mapped?.status;
    if (mapped?.resolved && ['canonical','resolved_external_id','resolved'].includes(status)
      && /^CHEBI:\d+$/i.test(String(mapped.resolved))
      && dictionary.exact.has(String(mapped.resolved).toUpperCase())) {
      return {status:'verified_crossref',id:dictionary.exact.get(String(mapped.resolved).toUpperCase()),provenance:status || 'external'};
    }
    return { status:'unmapped_identifier',id:null,provenance:'unresolved' };
  }
  const found = dictionary.names.get(key(original));
  if (found?.size === 1) return { status:'name_only',id:[...found][0],provenance:'curated_name' };
  if (found?.size > 1) return { status:'ambiguous_name',id:null,provenance:'curated_name' };
  return { status:'unrecognized',id:null,provenance:'unresolved' };
}
function measurement(row) {
  const usable = Number.isFinite(row.effect)
    && ['log2','as_supplied','transformed_unknown'].includes(row.effectScale);
  return { feature:row.feature, effect:usable ? row.effect : null,
    effectScale:row.effectScale || null,
    qValue:Number.isFinite(row.qValue) ? row.qValue : null,
    status:usable ? 'usable' : 'unavailable_effect_scale' };
}
/** Build an auditable network from explicit exact ChEBI IDs or unambiguous names only. */
export function buildFocusedMetabolicNetwork(result, central) {
  const { chemicals,genes } = toRegistry(central || {metabolites:[],enzymes:[]});
  const dictionary = identityDictionary(chemicals);
  const mappings = new Map((result?.identifierResolution?.metabolomics?.mappings || []).map((m)=>[m.original,m]));
  const userAnnotations = new Map((result?.userMetaboliteAnnotations || []).map((m)=>[m.feature,m.chebi]));
  const audit = [];
  const matched = new Map();
  for (const row of result?.layers?.metabolomics?.rows || []) {
    const identity = resolveIdentity(row.feature,dictionary,mappings,userAnnotations);
    const m = measurement(row);
    audit.push({input:row.feature,...identity,usable:m.status === 'usable',
      label:identity.id ? chemicals.get(identity.id)?.labelFr : null,
      chebi:identity.id ? chemicals.get(identity.id)?.chebi : null});
    if (!identity.id || m.status !== 'usable') continue;
    const prior = matched.get(identity.id);
    // Multiple distinct features for the same metabolite are NOT pooled or cherry-picked.
    if (!prior) matched.set(identity.id,{...m,matchStatus:identity.status,features:[row.feature]});
    else matched.set(identity.id,{...prior,effect:null,status:'duplicate_features',
      features:[...prior.features,row.feature]});
  }
  const geneMeasurements = new Map();
  const transcriptRows = result?.layers?.transcriptomics?.rows || [];
  for (const e of genes.values()) {
    const wanted = new Set([e.name,...e.aliases].map(key));
    const candidates = transcriptRows.filter((row)=>wanted.has(key(row.feature)) &&
      Number.isFinite(row.effect) && row.effectScale === 'log2');
    if (!candidates.length) continue;
    if (candidates.length !== 1) {
      geneMeasurements.set(e.id,{ status:'duplicate_features',effect:null,
        features:candidates.map((c)=>c.feature) });
    } else geneMeasurements.set(e.id,{...measurement(candidates[0]),matchStatus:'gene_symbol'});
  }
  const nodes = new Map();
  for (const mol of chemicals.values()) nodes.set(mol.id,{...mol,measurement:matched.get(mol.id) || null});
  for (const e of genes.values()) nodes.set(e.id,{...e,measurement:geneMeasurements.get(e.id) || null});
  const allLinks = [];
  const known = (a,b) => nodes.has(a) && nodes.has(b);
  for (const [a,b] of central?.edges || []) if (known(a,b))
    allLinks.push({from:a,to:b,relation:'schematic'});
  for (const pair of extraLinks) {
    const [a,b]=pair.split('>');
    if (known(a,b)) allLinks.push({from:a,to:b,relation:'schematic'});
  }
  for (const [symbol,anchor] of Object.entries(geneAnchor)) {
    const id='gene:'+symbol;
    if (known(id,anchor)) allLinks.push({from:id,to:anchor,relation:'gene_context'});
  }
  const regions = regionSpecs.map(([id,labelFr,labelEn,metabs,ens]) => {
    const ids = new Set([...list(metabs), ...list(ens).map((x)=>'gene:'+x.toLowerCase())]);
    const measured = [...ids].filter((n)=>nodes.get(n)?.measurement?.effect !== null &&
      nodes.get(n)?.measurement?.effect !== undefined);
    return {id,labelFr,labelEn,ids:[...ids].filter((n)=>nodes.has(n)),
      measuredIds:measured,measuredCount:measured.length};
  });
  return {
    nodes:[...nodes.values()], links:allLinks, regions, audit,
    coverage:{metabolites:chemicals.size,genes:genes.size,links:allLinks.length,regions:regions.length},
    caveat:'Only exact curated ChEBI IDs are guaranteed. Unambiguous names are explicitly marked name-only. Simplified pathway links are schematic rather than atom-balanced biochemical reactions or flux estimates.'
  };
}
/**
 * Display a small induced subgraph of regions containing observed molecular
 * effects. The hop count adjusts CONTEXT, not the statistical evidence.
 */
export function focusMetabolicRegion(network, regionId='auto', hops=1) {
  if (!network) return [];
  const byId=new Map(network.nodes.map((n)=>[n.id,n]));
  const depth=Math.max(0,Math.min(2,Number(hops)||0));
  const selected=regionId==='auto' ? (()=>{
    // Prefer specific measured regions over broad shared precursors (e.g. Trp
    // participates in both kynurenine and serotonin routes). Show a second
    // region only when it adds at least one previously unseen measured anchor.
    const covered=new Set();
    return network.regions.filter((r)=>r.measuredCount>0)
      .sort((a,b)=>b.measuredCount-a.measuredCount || a.id.localeCompare(b.id))
      .filter((region)=>{
        const adds=region.measuredIds.some((id)=>!covered.has(id));
        if (adds) for (const id of region.measuredIds) covered.add(id);
        return adds;
      });
  })() : network.regions.filter((r)=>r.id===regionId);
  return selected.map((region)=>{
    const allowed=new Set(region.ids);
    const base=new Set(region.measuredIds);
    const selectedIds=new Set(base);
    for(let step=0;step<depth;step++){
      const next=new Set(selectedIds);
      for(const edge of network.links){
        if (!allowed.has(edge.from)||!allowed.has(edge.to)) continue;
        if(selectedIds.has(edge.from))next.add(edge.to);
        if(selectedIds.has(edge.to))next.add(edge.from);
      }
      for(const id of next)selectedIds.add(id);
    }
    const ordered=region.ids.filter((id)=>selectedIds.has(id));
    const links=network.links.filter((e)=>selectedIds.has(e.from) && selectedIds.has(e.to) &&
      allowed.has(e.from) && allowed.has(e.to));
    return {...region,nodes:ordered.map((id)=>byId.get(id)),links};
  });
}
