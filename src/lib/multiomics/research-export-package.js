// @ts-nocheck
// Researcher-focused, offline ZIP dossier. No dependence on a remote service,
// publication API, or R backend. Browser statistics and R fits are identified
// separately. The ZIP contains raw study inputs and may be CONFIDENTIAL.
import {buildResearchMethodsManifest} from './research-report.js';
import {referenceReplayEntries} from './reference-replay.js';
import {RESEARCH_R_FIGURE_SOURCE} from './r-figure-source.js';
import {RESEARCH_RUN_ALL_R_SOURCE} from './r-run-all-source.js';
import {RESEARCH_R_MARKDOWN_SOURCE} from './r-markdown-source.js';
import {RESEARCH_R_MARKDOWN_RENDERER} from './r-markdown-renderer-source.js';
import {buildDownstreamExchange} from './downstream-exchange.js';
import {parseReferenceCsv} from './research-reference-csv.js';
import {RESEARCH_DOSSIER_VERIFY_PY,RESEARCH_DOSSIER_README} from './research-archive-verifier.js';

const encoder=new TextEncoder();
const CSV_COLUMNS=[
  ['feature','feature'],['effect','effect'],['effect_scale','effectScale'],
  ['ci95_low','ciLow'],['ci95_high','ciHigh'],
  ['exponentiated_effect','exponentiatedEffect'],
  ['exponentiated_ci95_low','exponentiatedCiLow'],
  ['exponentiated_ci95_high','exponentiatedCiHigh'],
  ['fold_ratio_if_log2','foldRatio'],['p_value','pValue'],
  ['q_value','qValue'],['n','n'],['n_reference','nReference'],
  ['n_comparison','nComparison'],['model','model']
];
const UTF8_FLAG=0x0800;
const MAX_ZIP_BYTES=130*1024*1024;
const LOCALES=['fr','en'];
/** The R backend currently exports method.top via top_frame(default 100).
 * It MUST NOT be passed off as complete reference R differential data.
 * This is separate from the complete browser exploratory result tables. */
export function referencePreviewAudit(result){
  const methods=result?.referenceBackend?.methods;
  if(!methods||typeof methods!=='object')return [];
  return Object.entries(methods).filter(([,m])=>m && typeof m==='object')
    .map(([key,m])=>{
      const table=m.referenceTable&&typeof m.referenceTable==='object'?m.referenceTable:m;
      const complete=table.referenceTableComplete===true &&
        table.referenceTableStatus==='complete';
      return {
        method_key:key,
        reported_name:typeof m.method==='string'?m.method:key,
        status:typeof m.status==='string'?m.status:'unknown',
        preview_rows:Array.isArray(m.top)?m.top.length:null,
        complete_reference_table_exported:complete,
        reference_table_rows:Number.isInteger(table.referenceTableRows)?
          table.referenceTableRows:null,
        reference_export_status:typeof table.referenceTableStatus==='string'?
          table.referenceTableStatus:'preview_only',
        local_archive_id:typeof table.referenceTableLocalArchiveId==='string' &&
          /^pmx_reference_[A-Za-z0-9._-]+$/.test(table.referenceTableLocalArchiveId)?
          table.referenceTableLocalArchiveId:null,
        local_archive_sha256:typeof table.referenceTableLocalArchiveSha256==='string' &&
          /^[a-f0-9]{64}$/.test(table.referenceTableLocalArchiveSha256)?
          table.referenceTableLocalArchiveSha256:null,
        paired_design:m.pairedDesignContract&&
          m.pairedDesignContract.pairing==='same_subject_two_conditions'
          ? {
              pairing:'same_subject_two_conditions',
              design:m.pairedDesignContract.design||null,
              independent_biological_subjects:m.pairedDesignContract.independent_biological_subjects??null,
              biological_pairs:m.pairedDesignContract.biological_pairs??null,
              model_samples:m.pairedDesignContract.model_samples??null,
              residual_degrees_of_freedom:m.pairedDesignContract.residual_degrees_of_freedom??null,
              contrast:m.pairedDesignContract.explicit_primary_contrast??null,
              contrast_was_prespecified:m.pairedDesignContract.contrast_was_prespecified===true,
              certified:false
            }:null,
        note:complete
          ?'Complete original-identifier R table is included in reference_r/, separate from browser exploration.'
          :'Not a complete fitted R feature table; top-frame previews cannot be cited as a full differential analysis.'
      };
    });
}
/** Extract only certified-complete-by-count R method CSV from a consented
 * local backend response, never mislabel browser tables as R output. */
export function fullReferenceMethodCsvs(result){
  const methods=result?.referenceBackend?.methods;
  if(!methods||typeof methods!=='object')return [];
  const collected=[];
  for(const [key,method] of Object.entries(methods)){
    if(!method||typeof method!=='object')continue;
    const table=method.referenceTable&&typeof method.referenceTable==='object'?
      method.referenceTable:method;
    if(table.referenceTableComplete===true && method.status!=='ok')
      throw Error('Reference '+key+' claims a complete fit with non-success status.');
    if(table.referenceTableComplete!==true){
      if(typeof table.referenceTableCsv==='string'&&table.referenceTableCsv.length)
        throw Error('Reference '+key+' has raw rows without verified completeness.');
      continue;
    }
    if(table.referenceTableStatus!=='complete'||typeof table.referenceTableCsv!=='string')
      throw Error('Reference '+key+' falsely claims a complete table.');
    if(!/^[A-Za-z0-9_]+$/.test(key))throw Error('Unsafe R method name in archive');
    const parsed=parseReferenceCsv(table.referenceTableCsv);
    if(!parsed.headers.includes('feature'))throw Error('Reference '+key+' has no original feature ID');
    if(!Number.isInteger(table.referenceTableRows)||table.referenceTableRows!==parsed.rows.length||
        table.referenceTableExportedRows!==parsed.rows.length)
      throw Error('Reference '+key+' row count mismatch');
    if(table.referenceTableBytes!==undefined &&
       table.referenceTableBytes!==encoder.encode(table.referenceTableCsv).length)
      throw Error('Reference '+key+' byte count mismatch');
    if(parsed.rows.some(row=>!row.feature)||new Set(parsed.rows.map(row=>row.feature)).size!==parsed.rows.length)
      throw Error('Reference '+key+' missing or duplicate fitted feature IDs');
    collected.push({name:'reference_r/'+key+'.csv',csv:table.referenceTableCsv,
      method_key:key,rows:parsed.rows.length,columns:parsed.headers,
      original_features:table.referenceTableOriginalFeatures??null,
      excluded_before_fit:table.referenceTableExcludedBeforeFit??null});
  }
  return collected;
}

function requireThat(value,message){if(!value)throw Error('Research dossier: '+message)}
function plain(value){return value==null||typeof value==='number'&&!Number.isFinite(value)?'':String(value)}
function decimal(value){return typeof value==='number'&&Number.isFinite(value)?value:''}
/** Preserve scientific identifiers in the JSON source. In the spreadsheet
 * CSV, avoid formula execution by explicitly prefixing unsafe text IDs.
 * This is disclosed in the manifest; it must never be a silent ID rewrite. */
function spreadsheetId(value){
  const s=plain(value);
  return /^[\s]*[=+@-]/.test(s)?"'"+s:s;
}
function cell(value){
  const s=plain(value);
  return /[",\r\n]/.test(s)?'"'+s.replaceAll('"','""')+'"':s;
}
/** Complete UTF-8 RFC 4180 compatible table, CRLF, Excel-safe identifiers. */
export function fullLayerCsv(rows){
  requireThat(Array.isArray(rows),'Results rows must be an array');
  const out=[CSV_COLUMNS.map(x=>x[0]).join(',')];
  for(const r of rows){
    requireThat(r&&typeof r.feature==='string'&&r.feature.length,
      'Missing measured feature identifier');
    for(const field of ['pValue','qValue']){
      if(r[field]!=null)
        requireThat(typeof r[field]==='number'&&Number.isFinite(r[field])&&
          r[field]>=0&&r[field]<=1,'Invalid '+field+' for '+r.feature);
    }
    out.push(CSV_COLUMNS.map(([heading,key])=>{
      const value=key==='feature'?spreadsheetId(r.feature):
        key==='effect'||key==='ciLow'||key==='ciHigh'||
        key==='exponentiatedEffect'||key==='exponentiatedCiLow'||
        key==='exponentiatedCiHigh'||key==='foldRatio'||
        key==='pValue'||key==='qValue'||key==='n'||
        key==='nReference'||key==='nComparison'?decimal(r[key]):r[key];
      return cell(value);
    }).join(','));
  }
  return '\uFEFF'+out.join('\r\n')+'\r\n';
}
const layerLabel={transcriptomics:['Transcriptomique','Transcriptomics'],
                  proteomics:['Protéomique','Proteomics'],
                  metabolomics:['Métabolomique','Metabolomics']};
function text(fr,en,lang){return lang==='en'?en:fr}
function print(v){if(Array.isArray(v))return v.join(', ')||'—';if(v===null||v===undefined||v==='')return '—';return String(v)}
/** A usable, conservative methods and interpretation document. */
export function methodsMarkdown(result,locale='fr'){
  requireThat(LOCALES.includes(locale),'Unknown language');
  requireThat(result?.engine&&result?.protocol&&result?.layers,
    'No completed omics analysis');
  const a=buildResearchMethodsManifest(result);
  const rPreviews=referencePreviewAudit(result);
  const fullR=fullReferenceMethodCsvs(result);
  const t=(fr,en)=>text(fr,en,locale);
  const h=[
    '# '+t('Analyse MultiOmics — méthodes et résultats','MultiOmics analysis — methods and results'),
    '',
    '> '+t('Dossier contenant des données potentiellement sensibles. Ne pas diffuser sans autorisation.',
      'This dossier may contain sensitive research data. Do not distribute without authorization.'),
    '',
    '## '+t('Résumé de l’étude','Study summary'),
    '',
    '| '+t('Champ','Field')+' | '+t('Valeur','Value')+' |',
    '| --- | --- |',
    '| '+t('Date du calcul','Analysis date')+' | '+print(result.generatedAt)+' |',
    '| '+t('Moteur réellement exécuté','Engine actually executed')+' | '+print(result.engine.name)+' |',
    '| '+t('Version du moteur','Engine version')+' | '+print(result.engine.version)+' |',
    '| '+t('Objectif déclaré','Declared objective')+' | '+print(a.objective)+' |',
    '| '+t('Type de plan','Study design')+' | '+print(a.experimental_design)+' |',
    '| '+t('Intention scientifique','Inference intent')+' | '+print(a.intended_inference)+' |',
    '| '+t('Sujets distincts','Distinct subjects')+' | '+print(a.independent_subjects)+' |',
    '| '+t('Groupes déclarés','Declared groups')+' | '+print(result.metadataSummary?.conditions)+' |',
    '| '+t('Statut de la revue','Review status')+' | '+print(a.scientific_assurance_status)+' |',
    '',
    '## '+t('Méthodes exécutées par type de données','Methods actually run for each data type'),
    '',
  ];
  for(const [layer,entry] of Object.entries(result.layers)){
    if(!layerLabel[layer])continue;
    const label=layerLabel[layer][locale==='en'?1:0];
    const audit=a.omics[layer]||{};
    const omitted=Number.isFinite(audit.omitted_features)&&audit.omitted_features>=0?
      audit.omitted_features:null;
    const omittedRatio=omitted!==null&&Number.isFinite(audit.features_input)&&
      audit.features_input>0?100*omitted/audit.features_input:null;
    const omittedDisplay=print(omitted)+(omittedRatio!==null?
      ' ('+omittedRatio.toFixed(1).replace('.',locale==='fr'?',':'.')+' %)':'');
    const nrows=Array.isArray(entry?.rows)?entry.rows.length:0;
    h.push('### '+label,'',
      '- '+t('Calcul présenté','Reported calculation')+': '+print(audit.reported_model),
      '- '+t('Données déclarées','Declared value type')+': '+print(audit.declared_value_type),
      '- '+t('Variables avant / après filtrage','Features before / after QC')+': '+print(audit.features_input)+' / '+print(audit.features_after_filter),
      '- '+t('Variables écartées par le QC navigateur','Features omitted by browser QC')+': '+omittedDisplay,
      '- '+t('Lignes du tableau exporté','Exported table rows')+': '+nrows,
      '- '+t('Proportion médiane manquante','Median missing fraction')+': '+print(audit.median_missing_fraction),
      '- '+t('Valeur maximale manquante par échantillon','Worst missing fraction per sample')+': '+print(audit.maximum_sample_missing_fraction),
      '- '+t('Moteur R de référence pour cette analyse','R reference method for this analysis')+': '+
           (audit.reference_executed?t('exécuté (sans certification)','executed (not certified)'):
             t('non exécuté','not executed')),
      '- '+t('Contrôle de multiplicité déclaré','Declared multiple-testing policy')+': '+print(audit.fdr_family),
      '',
      t('Le tableau CSV contient toutes les variables renvoyées par le moteur navigateur, sans troncature. Une exécution R de référence, lorsqu’elle existe, est conservée séparément dans analysis_complete.json ; les lignes navigateur ne doivent jamais être attribuées à DESeq2, limma ou mixOmics par défaut.',
        'The CSV contains every feature returned by the browser engine without truncation. When R reference computations exist, they are preserved separately in analysis_complete.json; browser result rows must not be attributed to DESeq2, limma or mixOmics by default.'),
      '');
  }


  // Internal outer-fold predictions are NOT a frozen independent external
  // validation. Report proper probability scores alongside ranking metrics.
  if(result.predictiveOutcome?.status==='ok' &&
     Number.isFinite(result.predictiveOutcome.metrics?.auc)){
    const prediction=result.predictiveOutcome.metrics;
    h.push('## '+t('Prédiction binaire — discrimination et calibration',
       'Binary prediction — discrimination and calibration'),'',
       '- AUC (nested CV): '+print(prediction.auc),
       '- Brier (nested CV): '+print(prediction.brier),
       '- log-loss (nested CV): '+print(prediction.logLoss),
       '',
       '> '+t('AUC élevée ≠ probabilités fiables. Le score de Brier mesure une erreur quadratique moyenne, pas à lui seul une validation de calibration. Comparer à une prédiction de référence, examiner courbes et pente/intercept de calibration, puis valider sur une cohorte externe sans réentraîner le modèle. Une recalibration peut dégrader les scores sur un autre centre.',
         'High AUC does not imply reliable probabilities. Brier is a mean squared probability error, not a stand-alone calibration certificate. Compare against a reference predictor, inspect calibration curves, slope and intercept, then validate on an independent cohort without retraining. Recalibration can worsen performance at another center.'),'');
  }
  if(a.omics?.transcriptomics?.null_rnaseq_calibration_review_required){
    h.push('## '+t('Limites de calibration RNA-seq','RNA-seq statistical calibration limits'), '',
      '> '+t('Des études de simulation sous nul global ont observé, dans certains scénarios, davantage de faux positifs que le seuil nominal de 5 % sous DESeq2. Cela ne prouve pas un défaut général de DESeq2, mais une q-value BH n’est pas à elle seule la preuve que la FDR est contrôlée pour ce plan. Vérifier les facteurs de taille, dispersions, contrastes, dépendances et hypothèses avec une revue statistique indépendante avant une conclusion confirmatoire.',
        'Global-null simulations showed excess false discoveries above the nominal 5% for DESeq2 in some tested settings. This is not evidence of a general DESeq2 defect, but BH q-values alone do not prove FDR calibration for this study design. Independently review size factors, dispersions, contrasts, dependence and assumptions before making confirmatory claims.'), '');
  }

  if(rPreviews.length){
    h.push('## '+t('Résultats de référence R — exhaustivité vérifiée par méthode',
      'R reference results — per-method completeness'),'',
      t('Les CSV de reference_r/ sont distincts des CSV exploratoires de tableaux/. Les identifiants R originaux sont conservés. Les fichiers R sont destinés à R/Python : un identifiant malveillant peut être interprété comme une formule par un tableur. N’ouvrez pas ces CSV non fiables directement dans Excel.',
        'The reference_r/ CSVs are distinct from exploratory tableaux/ CSVs. Exact original R identifiers are preserved. R CSVs are intended for R/Python; malicious identifiers can execute as spreadsheet formulas, so do not open untrusted R CSVs directly in Excel.'),
      '');
    for(const r of rPreviews)
      h.push('- '+r.method_key+': '+r.status+
        ', '+t('aperçu','preview')+'='+print(r.preview_rows)+
        ', '+t('export complet','complete export')+'='+
        (r.complete_reference_table_exported?t('oui','yes'):t('non','no'))+
        ', '+t('lignes','rows')+'='+print(r.reference_table_rows)+
        ', '+t('état','status')+'='+r.reference_export_status+
        (r.local_archive_id && r.reference_export_status==='too_large_archived_locally'?
          ', '+t('sauvegarde R locale (HORS ZIP)','local R-only archive (NOT IN ZIP)')+
          '='+r.local_archive_id+', SHA-256='+print(r.local_archive_sha256) :
          ''));
    h.push('');
    for(const r of rPreviews.filter(x=>x.paired_design)){
      const design=r.paired_design;
      h.push('**'+r.reported_name+' — '+t('plan apparié réel','actual paired design')+'**',
        '- '+t('Modèle','Model')+': \`'+print(design.design)+'\`',
        '- '+t('Unités biologiques indépendantes','Independent biological units')+': '+print(design.independent_biological_subjects),
        '- '+t('Paires','Pairs')+': '+print(design.biological_pairs),
        '- '+t('Prélèvements ajustés','Fitted samples')+': '+print(design.model_samples),
        '- '+t('Degrés de liberté résiduels','Residual degrees of freedom')+': '+print(design.residual_degrees_of_freedom),
        '- '+t('Contraste (comparaison vs référence)','Contrast (comparison vs reference)')+': '+print(design.contrast),
        '- '+t('Contraste déclaré avant calcul','Contrast declared before fitting')+': '+(design.contrast_was_prespecified?'yes':'no'),
        '- '+t('Certification scientifique automatique','Automatic scientific certification')+': **non / no**','');
    }
    if(rPreviews.some(r=>r.reference_export_status==='too_large_archived_locally'))
      h.push('> '+t('Des résultats R dépassant la limite JSON ont été conservés UNIQUEMENT sur le disque du backend local dans le dossier défini par PMX_MULTIOMICS_RESULTS_DIR. Ils ne sont PAS inclus dans ce ZIP. Récupérer le fichier fitted_reference.csv dans le sous-dossier ID indiqué et vérifier le SHA-256.',
        'Oversized R results were retained ONLY in the R backend operator\'s local PMX_MULTIOMICS_RESULTS_DIR. They are NOT included in this ZIP. Retrieve fitted_reference.csv from the recorded folder ID and verify SHA-256.'),'');
    if(rPreviews.some(r=>r.status==='ok'&&r.preview_rows!==null&&!r.complete_reference_table_exported))
      h.push('> '+t('Certaines analyses R ne disposent que d’un aperçu (ou dépassent la taille de transfert). Conserver leur CSV complet depuis le moteur R local avant toute publication.',
        'Some R analyses have only a preview (or exceed the transfer size limit). Preserve their full CSV from the local R engine before publication.'),'');
  }
  h.push('## '+t('Comprendre les colonnes exportées','Understanding exported columns'),'',
    '| '+t('Colonne','Column')+' | '+t('Définition','Definition')+' |',
    '| --- | --- |',
    '| \`feature\` | '+t('Identifiant de la mesure (valeur exacte dans le JSON).','Measured feature (exact identifier in JSON).')+' |',
    '| \`effect\` / \`effect_scale\` | '+t('Taille d’effet et échelle associée ; le signe dépend du contraste déclaré.','Effect size and associated scale; the sign depends on the declared contrast.')+' |',
    '| \`ci95_low\` / \`ci95_high\` | '+t('Bornes disponibles de l’intervalle, sinon vides.','Interval bounds if available, otherwise blank.')+' |',
    '| \`p_value\` | '+t('p-value du test indiqué, pas la probabilité que l’hypothèse soit vraie.','Reported test p-value, not the probability the hypothesis is true.')+' |',
    '| \`q_value\` | '+t('p-value ajustée pour tests multiples, sans prouver une découverte biologique.','Multiple-testing-adjusted p-value, not proof of a biological discovery.')+' |',
    '| \`n\` / \`n_reference\` / \`n_comparison\` | '+t('Observations contribuant au calcul : ne pas confondre avec les sujets indépendants.','Contributing observations: do not confuse with independent subjects.')+' |',
    '| \`model\` | '+t('Nom de la méthode correspondant à la ligne.','Method name corresponding to the row.')+' |',
    '',
    t('Les CSV sont directement importables dans R ou Python. Vérifier les hypothèses et les limites avant de reprendre les estimations dans un article.',
      'CSV tables can be imported in R or Python. Check assumptions and limitations before citing estimates in a paper.'),
    '');
  h.push('## '+t('Interprétation prudente','Interpretation boundaries'),'',
    '- '+t('Une p-value ou une q-value n’établit ni causalité ni validation biologique. Les comparaisons sont limitées au plan déclaré.',
      'A p-value or q-value does not demonstrate causality or biological validation. Comparisons are limited to the declared study design.'),
    '- '+t('Un résultat exploratoire ne devient pas confirmatoire parce qu’un moteur R est disponible.',
      'Exploratory evidence is not confirmatory solely because R is available.'),
    '- '+t('Les variables dont le nom peut exécuter une formule dans un tableur sont préfixées par une apostrophe dans le CSV ; les identifiants exacts sont préservés dans analysis_complete.json.',
      'Feature names that could execute spreadsheet formulas are prefixed by an apostrophe in CSV; exact names are preserved in analysis_complete.json.'),
    '- '+t('Les valeurs absentes sont des cellules vides, jamais des zéros inventés.',
      'Missing estimates are blank cells, never fabricated zeros.'),
    '- '+t('Les bases de données externes et leurs versions ne peuvent être reproduites que si leurs identifiants et versions sont capturés dans analysis_complete.json.',
      'External database results are only independently reproducible if their identifiers and releases are captured in analysis_complete.json.'),
    '- '+t('Une analyse utilisant un serveur R nécessite également ses versions de packages et son environnement pour un rejeu indépendant.',
      'Analyses using the R backend also require its package versions and environment for independent reruns.'),
    '',
    '## '+t('Organisation du dossier','Folder contents'),'',
    '- `tableaux/*.csv`: '+t('tableaux R natifs originaux lorsque la couche affiche les chiffres R, sinon tableaux exploratoires JS clairement identifiés','original native R full tables when a layer uses R values, otherwise clearly identified exploratory browser tables'),
    '- `methods/methods.md`: '+t('ce texte modifiable','this editable text'),
    '- `methods/methods_audit.json`: '+t('méthodes et QC sans matrices brutes','method/QC summary without raw matrices'),
    '- `provenance/analysis_complete.json`: '+t('résultat complet, y compris traces éventuelles du backend R','complete analysis including any R backend trace'),
    '- `source/projet_reproductible.json`: '+t('fichiers originaux et paramètres pour réimport puis recalcul ; ne contient pas de résultat considéré nouveau','original input files and settings for manual rerun; no archived result reused as a new analysis'),
    '- `provenance/manifest.json`: '+t('empreintes SHA-256 des autres fichiers','SHA-256 digests of the other files'),
    '- `verification/verifier_dossier.py`: '+t('script autonome Python 3 de vérification sans internet','offline Python 3 script for independent integrity checks'),
    '',
    t('Le dossier ne constitue pas une certification statistique ni clinique.',
      'This dossier is not a statistical or clinical certification.'),
    ''
  );
  return h.join('\n');
}
let table=null;
function crc32(data){
  if(!table){table=new Uint32Array(256);for(let n=0;n<256;n++){let c=n;for(let k=0;k<8;k++)c=(c&1)?0xedb88320^(c>>>1):c>>>1;table[n]=c>>>0}}
  let crc=0xffffffff;for(const v of data)crc=table[(crc^v)&255]^(crc>>>8);
  return (crc^0xffffffff)>>>0;
}
function u16(view,offset,x){view.setUint16(offset,x,true)}
function u32(view,offset,x){view.setUint32(offset,x>>>0,true)}
/** Stable, standards-compliant uncompressed ZIP with UTF-8 file names;
 * deliberately small-scope to avoid browser supply-chain dependencies. */
export function zipEntries(entries){
  requireThat(Array.isArray(entries)&&entries.length>0&&entries.length<=60,
    'Invalid number of research files');
  const names=new Set(),local=[],central=[];
  let offset=0,sourceBytes=0;
  for(const entry of entries){
    const name=entry?.name;
    requireThat(typeof name==='string'&&/^[A-Za-z0-9._/-]+$/.test(name)&&
      !name.startsWith('/')&&!name.split('/').includes('..')&&!names.has(name),
      'Unsafe or duplicate ZIP member path');
    names.add(name);
    const filename=encoder.encode(name);
    const data=entry.data instanceof Uint8Array?entry.data:encoder.encode(String(entry.data));
    requireThat(filename.length<=65535,'ZIP file name too long');
    sourceBytes+=data.length;
    requireThat(sourceBytes<=MAX_ZIP_BYTES,'Too much research data for browser ZIP; use individual exports');
    const crc=crc32(data),head=new Uint8Array(30+filename.length),v=new DataView(head.buffer);
    u32(v,0,0x04034b50);u16(v,4,20);u16(v,6,UTF8_FLAG);
    u16(v,8,0);u16(v,10,0);u16(v,12,0x0021); // 1980-01-01, ZIP DOS epoch
    u32(v,14,crc);u32(v,18,data.length);u32(v,22,data.length);
    u16(v,26,filename.length);u16(v,28,0);head.set(filename,30);
    local.push(head,data);
    const center=new Uint8Array(46+filename.length),c=new DataView(center.buffer);
    u32(c,0,0x02014b50);u16(c,4,20);u16(c,6,20);
    u16(c,8,UTF8_FLAG);u16(c,10,0);u16(c,12,0);u16(c,14,0x0021);
    u32(c,16,crc);u32(c,20,data.length);u32(c,24,data.length);
    u16(c,28,filename.length);u16(c,30,0);u16(c,32,0);
    u16(c,34,0);u16(c,36,0);u32(c,38,0);u32(c,42,offset);
    center.set(filename,46);central.push(center);
    offset+=head.length+data.length;
    requireThat(offset<0xffffffff,'ZIP size > 4 GiB not supported');
  }
  const centralBytes=central.reduce((sum,x)=>sum+x.length,0);
  const trailer=new Uint8Array(22),e=new DataView(trailer.buffer);
  u32(e,0,0x06054b50);u16(e,4,0);u16(e,6,0);
  u16(e,8,entries.length);u16(e,10,entries.length);
  u32(e,12,centralBytes);u32(e,16,offset);u16(e,20,0);
  return new Blob([...local,...central,trailer],{type:'application/zip'});
}
async function sha256(bytes){
  const result=await globalThis.crypto.subtle.digest('SHA-256',bytes);
  return [...new Uint8Array(result)].map(x=>x.toString(16).padStart(2,'0')).join('');
}
/** PortableProject is constructed from the EXACT original source File objects,
 * keeping any sensitive rows exclusively in the user's downloaded archive. */
/** @param {{result:any,portableProject:any,referenceReplay?:any,figureSnapshots?:Array<{name:string,data:string}>,locale?:string}} options */
export async function buildResearchBundle({result,portableProject,referenceReplay=null,figureSnapshots=[],locale='fr'}){
  requireThat(result?.engine&&result?.layers&&result?.protocol,'Missing completed analysis');
  requireThat(portableProject?.format==='pmx-multiomics-portable-project'&&
    portableProject?.schemaVersion===1&&portableProject?.inputs?.metadata,
    'Original study files required (no provenance-free package)');
  requireThat(!portableProject.archivedResult,
    'Bundle source project must not present archived results as a fresh analysis');
  const members=[];
  const add=(name,data)=>members.push({name,data});
  const count={};
  let num=0;
  for(const [layer,entry] of Object.entries(result.layers)){
    requireThat(Object.hasOwn(layerLabel,layer),'Unsupported data layer '+layer);
    requireThat(Array.isArray(entry?.rows),'Missing results table for '+layer);
    count[layer]=entry.rows.length;
    // A layer displayed with native R statistics must export the ORIGINAL R
    // CSV as its primary table. A reserialization would discard native
    // columns such as DESeq2 baseMean/lfcSE/stat and break exact provenance.
    if(entry.numericalProvenance==='exact_native_r_reference_table'){
      const methodKey=entry.primaryReferenceMethodKey;
      requireThat(typeof methodKey==='string'&&/^[A-Za-z0-9_]+$/.test(methodKey),
        'Native R primary layer lacks a safe method key');
      const native=result.referenceBackend?.methods?.[methodKey];
      const table=native?.referenceTable&&typeof native.referenceTable==='object'?
        native.referenceTable:native;
      requireThat(native?.status==='ok'&&table?.referenceTableComplete===true &&
        table.referenceTableStatus==='complete'&&typeof table.referenceTableCsv==='string',
        'Native R primary layer lacks its complete reference CSV');
      const parsed=parseReferenceCsv(table.referenceTableCsv);
      requireThat(parsed.rows.length===entry.rows.length &&
        parsed.rows.every((r,i)=>r.feature===entry.rows[i].feature),
        'Native R primary table differs from displayed feature order');
      add('tableaux/'+layer+'.csv',table.referenceTableCsv);
    }else{
      add('tableaux/'+layer+'.csv',fullLayerCsv(entry.rows));
    }
    num++;
  }
  requireThat(num>0,'No analyzable layer');
  const fullR=fullReferenceMethodCsvs(result);
  for(const reference of fullR)add(reference.name,reference.csv);
  // These are the ACTUAL R backend sources, JSON request and returned values.
  // Browser figures are only immutable snapshots: R never drew them.
  const replay=await referenceReplayEntries(referenceReplay);
  for(const entry of replay.entries)add(entry.name,entry.data);
  requireThat(Array.isArray(figureSnapshots)&&figureSnapshots.length<=12,'Invalid number of figure snapshots');
  for(const item of figureSnapshots){
    requireThat(/^figure_[0-9]{2}\.svg$/.test(item?.name)&&
      typeof item?.data==='string'&&item.data.length<2_000_000&&
      item.data.includes('<svg'), 'Invalid browser SVG snapshot');
    add('figures_browser/'+item.name,item.data);
  }
  add('methods/methods.md',methodsMarkdown(result,locale));
  add('verification/verifier_dossier.py',RESEARCH_DOSSIER_VERIFY_PY);
  add('verification/LISEZ-MOI.txt',RESEARCH_DOSSIER_README);
  add('methods/methods_audit.json',JSON.stringify(buildResearchMethodsManifest(result),null,2)+'\n');
  add('provenance/analysis_complete.json',JSON.stringify(result,null,2)+'\n');
  add('source/projet_reproductible.json',JSON.stringify(portableProject,null,2)+'\n');
  // R figure recipes use the SAME archived plot-ready values, rather than
  // separately refitting or silently altering the browser statistical output.
  // This is R-generated, scientifically similar plotting; not pixel identity.
  add('COMMENCER_ICI_FR_EN.md',[
    '# MultiOmics — Reproduire les résultats / Reproduce the results',
    '',
    'FR — Aucun changement du code n’est nécessaire :',
    '1. Conservez ce ZIP dans un dossier sécurisé : il contient vos données scientifiques.',
    '2. Extrayez entièrement le ZIP. Installez R puis le package R jsonlite.',
    '3. Dans ce dossier extrait, lancez : Rscript run_all.R',
    '   Pour le compte rendu pédagogique HTML commenté : Rscript render_notebook.R (installe rmarkdown uniquement si absent).',
    '4. Consultez research_reproduction_status.json, figures_r/output/ et les CSV des valeurs tracées.',
    '5. Si le calcul de référence R était disponible, vérifiez PASS_R_REFERENCE_ONLY.',
    '   Sinon, les chiffres exploratoires JS ne sont PAS présentés comme reproduits en R.',
    '6. Les figures générées par R sont scientifiquement similaires mais peuvent différer graphiquement des SVG interactifs.',
    '',
    'EN — No code edits required:',
    '1. Keep the archive secure: it can contain confidential research data.',
    '2. Extract all files. Install R and the R package jsonlite.',
    '3. Open a terminal inside the extracted folder and run: Rscript run_all.R',
    '   For a fully commented HTML scientific notebook: Rscript render_notebook.R (installs rmarkdown only if missing).',
    '4. Read research_reproduction_status.json, figures_r/output/ and the plotted-values CSVs.',
    '5. PASS_R_REFERENCE_ONLY means the native R reference result was numerically replayed.',
    '   Browser-only exploratory numbers are never misrepresented as a native R method.',
    '6. R figure equivalents use the same plot-ready values, not identical SVG pixel geometry.',
    '',
    'Pour le rapport HTML : RStudio inclut généralement Pandoc. Sous Linux, la compilation de rmarkdown peut exiger libuv1-dev, libssl-dev, libcurl4-openssl-dev et libxml2-dev.',
    'To render the HTML report: RStudio usually includes Pandoc. Linux installations may require libuv1-dev, libssl-dev, libcurl4-openssl-dev and libxml2-dev.',
    'Ubuntu/Debian (system administrator): sudo apt-get install libuv1-dev libssl-dev libcurl4-openssl-dev libxml2-dev pandoc',
    'Windows: open PowerShell in the extracted folder and run Rscript.exe run_all.R.',
    'If Rscript is not on PATH, run it from the R installation bin directory.',
    'For strict numerical native-R replay, restore the ORIGINAL R and Bioconductor package versions.',
    'A new/latest installation is NOT evidence of matching the original statistical results.',
    'Report a FAIL or NOT_CHECKED; do not call an unverified replay reproducible.',
    '',
  ].join('\n'));
  add('ANALYSE_REPRODUCTIBLE.Rmd',RESEARCH_R_MARKDOWN_SOURCE);
  add('render_notebook.R',RESEARCH_R_MARKDOWN_RENDERER);
  add('run_all.R',RESEARCH_RUN_ALL_R_SOURCE);
  add('figures_r/generate_figures.R',RESEARCH_R_FIGURE_SOURCE);
  add('figures_r/README.md',[
    '# Recreate the scientific display figures using R',
    '',
    'From the extracted ZIP root run:',
    '',
    'Rscript figures_r/generate_figures.R provenance/analysis_complete.json figures_r/output',
    '',
    'Requires R + jsonlite and an SVG-capable R build (PDF is always generated).',
    'The program uses original heatmap row z-scores, annotations, and reported',
    'model effect values from analysis_complete.json without statistical refitting.',
    'It renders a heatmap for each available layer, a circular effect plot and',
    'a central-carbon network when those display datasets exist.',
    'Plots are scientifically similar to browser/Svelte views, NOT pixel-identical.',
    'It does not reproduce arbitrary interactive tabs or local selections.',
    'No biological significance or pathway activity is computed by plotting.',
    'Inspect output/figure_manifest.json: every saved PDF/SVG identifies its input.',
    'The original browser SVG snapshots are in figures_browser/, not produced by R.',
    'Where an R reference fit ran, the exact native numerical replay is separately',
    'available under reproduction/ (requires its own source/environment attestation).',
    ''
  ].join('\n'));
  // Exchange is an optional derivative of the EXACT uploaded matrices. Do not
  // block any existing mono-omic ZIP export when their layout differs.
  let downstreamExchange={status:'not_available',reason:null};
  try {
    const exchange=await buildDownstreamExchange(portableProject);
    for(const entry of exchange.entries)add(entry.name,entry.data);
    downstreamExchange={status:'included',schema_version:'1.0.0',matrices:exchange.manifest.matrices.map(m=>m.omic)};
  } catch(error) {
    downstreamExchange={status:'not_available',reason:error instanceof Error?error.message:'Unknown incompatibility'};
  }
  const files={};
  for(const {name,data} of members){
    const bytes=encoder.encode(data);
    files[name]={sha256:await sha256(bytes),bytes:bytes.length};
  }
  const manifest={
    schema:'pmx-multiomics-research-dossier/v1',
    created_at:result.generatedAt||null,
    analysis_engine:result.engine.name,
    analysis_engine_version:result.engine.version,
    ui_language:locale,
    input_file_sha256:Object.fromEntries(
      Object.entries(portableProject.inputs).map(([key,f])=>[key,f.sha256])),
    browser_layer_rows:count,
    primary_table_sources:Object.fromEntries(
      Object.entries(result.layers).map(([layer,entry])=>[
        layer,entry.numericalProvenance==='exact_native_r_reference_table'
          ? {engine:'native_R',method:entry.primaryReferenceMethodKey}
          : {engine:'browser_exploratory',method:entry.inferenceMethod||null}
      ])),
    r_source_replay:replay.status,
    browser_svg_snapshots:{count:figureSnapshots.length,source:'original_rendered_Svelte_SVG',r_replayable:false},
    exact_r_numeric_and_figure_equivalence_certified:false,
    downstream_exchange:downstreamExchange,
    reference_method_previews:referencePreviewAudit(result),
    r_reference_full_feature_tables_exported:fullR.length>0,
    r_reference_full_tables:fullR.map(({name,method_key,rows,columns,original_features,excluded_before_fit})=>
      ({name,method_key,rows,columns,original_features,excluded_before_fit})),
    note_on_r_reference:'Only reference_r/*.csv with independently checked row counts are full R method tables. Missing or oversized R methods remain previews and are NOT publication-complete.',

    files,
    inference_scope:'browser exploratory tables; any R backend output is separately recorded in full JSON',
    clinical_validation:false,
    confirmatory_certification:false,
    real_backend_proved_by_this_package:false,
    original_inputs_included:true,
    sensitive_research_data:true,
    provenance_notes:[
      'File SHA-256 values refer to the uncompressed UTF-8 member bytes.',
      'Result CSV spreadsheet protection prefixes potentially executable feature identifiers with a single apostrophe; raw identifiers are unchanged in analysis_complete.json.',
      'All layer features are exported, not only the HTML preview top 50.',
      'Run the included offline Python verifier to check the original files and complete CSV counts; this does not certify statistical validity.',
      'ZIP source project contains exact original uploaded CSV/TSV and parameters, but not a fresh inference when re-imported.',
      'For R reference methods, also archive version-locked backend and package library.',
      'For R reference differential results, obtain full fitted tables separately; the backend currently returns top-frame previews of up to 100 rows, not exhaustive reference tables.',
      'Never publish the bundle publicly without reviewing subjects, assays and patient confidentiality.'
    ]
  };
  add('provenance/manifest.json',JSON.stringify(manifest,null,2)+'\n');
  return {
    filename:'multiomics_dossier_recherche_v1.zip',
    blob:zipEntries(members),
    manifest
  };
}
