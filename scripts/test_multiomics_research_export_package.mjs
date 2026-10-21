import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
import {parseReferenceCsv} from '../src/lib/multiomics/research-reference-csv.js';
import {File} from 'node:buffer';
import {createPortableProject} from '../src/lib/multiomics/portable-project.js';
import {buildResearchBundle,fullLayerCsv,methodsMarkdown,zipEntries} from '../src/lib/multiomics/research-export-package.js';

if(!globalThis.File)globalThis.File=File;
if(!globalThis.crypto)globalThis.crypto=crypto.webcrypto;
const subject='PRIVATE_SUBJECT_NEVER_LOG';
const rows=[
 {feature:'GENE,1',effect:-0.21,effectScale:'log2',pValue:0.02,qValue:0.07,approxPValue:0.021,approxQValue:0.071,n:8,nReference:4,nComparison:4,model:'Welch, robust'},
 {feature:'=HYPERLINK("http://malicious.test")',effect:1.25,pValue:null,qValue:null,n:8,model:'reference'},
 {feature:'G"ene\nLine',effect:0,qValue:1,pValue:0,model:'normal'},
 ...Array.from({length:53},(_,i)=>({feature:'GENE_'+i,effect:i-5,
   pValue:(i+1)/60,qValue:(i+2)/60,model:'test'}))
];
const csv=fullLayerCsv(rows);
assert.ok(csv.startsWith('\uFEFFfeature,effect,effect_scale'));
assert.match(csv,/approx_p_value_browser_uncalibrated,approx_q_value_browser_uncalibrated/);
assert.equal(csv.split('\r\n').filter(Boolean).length,57); // header + 56 rows (quoted LF does not split CRLF)
assert.ok(csv.includes('"GENE,1"'));
assert.ok(csv.includes('G""ene'+String.fromCharCode(10)+'Line')); // doubled quote in feature
assert.ok(csv.includes("'=HYPERLINK"));
assert.ok(csv.includes(',-0.21,'));
assert.ok(csv.includes(',0,'));
assert.ok(!csv.includes('undefined'));
assert.throws(()=>fullLayerCsv([{feature:'G1',pValue:1.2}]),/Invalid pValue/);
assert.throws(()=>fullLayerCsv([{feature:'G1',approxQValue:-0.1}]),/Invalid approxQValue/);
assert.throws(()=>fullLayerCsv([{feature:'',pValue:0.1}]),/Missing measured feature/);

const result={
 generatedAt:'2026-10-13T14:02:00Z',
 engine:{name:'PMx deterministic browser',version:'2.1',analysisMode:'single_omic'},
 protocol:{objective:'groups',analysisIntent:'exploratory',designType:'independent'},
 metadataSummary:{subjects:8,conditions:['control','treated']},
 inputManifest:{sourceSha256:'a'.repeat(64)},
 layers:{transcriptomics:{rows,inferenceMethod:'Welch exploratory screening',
  qc:{featuresBefore:60,featuresAfter:56,medianMissingFraction:0,
    preprocessingAudit:{declaredValueType:'normalized'}}}},
 scientificAssurance:{status:'independent_review_required'},
 privateSubject:subject
};
const files={metadata:new File(['subject_id,sample_id,condition\n'+subject+',S01,control\n'],
 'private_metadata.csv',{type:'text/csv'}),
 transcriptomics:new File(['feature_id,S01\nGENE1,10\n'],'rna.csv',{type:'text/csv'})};
const portable=await createPortableProject({files,settings:{objective:'groups',analysisIntent:'exploratory'},result:null});
const md=methodsMarkdown(result,'fr');
assert.match(md,/Transcriptomique/);
assert.match(md,/56/);
assert.match(md,/Variables écartées par le QC navigateur: 4 \(6,7 %\)/);
assert.match(md,/non exécuté/);
assert.match(md,/toutes les variables/);
assert.throws(()=>methodsMarkdown(result,'invalid'),/Unknown language/);
const rawRnaResult={...result,layers:{...result.layers,
  transcriptomics:{...result.layers.transcriptomics,
    qc:{...result.layers.transcriptomics.qc,
      preprocessingAudit:{declaredValueType:'raw_counts'}}}}};
const rawRnaMethods=methodsMarkdown(rawRnaResult,'fr');
assert.match(rawRnaMethods,/Limites de calibration RNA-seq/);
assert.match(rawRnaMethods,/faux positifs/);
assert.doesNotMatch(md,/Limites de calibration RNA-seq/);
const rawRnaManifest=JSON.stringify((await buildResearchBundle({
  result:rawRnaResult,
  portableProject:await createPortableProject({
    files,settings:{objective:'groups'},result:null
  })
})).manifest);
assert.match(rawRnaManifest,/confirmatory_certification/);

const bundle=await buildResearchBundle({result,portableProject:portable,locale:'en'});
assert.equal(bundle.filename,'multiomics_dossier_recherche_v1.zip');
assert.ok(bundle.blob.size>1000);
assert.equal(bundle.manifest.browser_layer_rows.transcriptomics,56);
assert.equal(bundle.manifest.confirmatory_certification,false);
assert.equal(bundle.manifest.clinical_validation,false);
assert.equal(bundle.manifest.sensitive_research_data,true);
assert.equal(bundle.manifest.files['tableaux/transcriptomics.csv'].sha256.length,64);
// Check ZIP structures, byte-for-byte SHA-256, no external libraries.
const b=Buffer.from(await bundle.blob.arrayBuffer());
let offset=0,seen={},names=[];
while(offset<b.length && b.readUInt32LE(offset)===0x04034b50){
 const method=b.readUInt16LE(offset+8),crc=b.readUInt32LE(offset+14);
 const len=b.readUInt32LE(offset+18),nlen=b.readUInt16LE(offset+26);
 const extra=b.readUInt16LE(offset+28),name=b.subarray(offset+30,offset+30+nlen).toString();
 const file=b.subarray(offset+30+nlen+extra,offset+30+nlen+extra+len);
 assert.equal(method,0);
 assert.ok(crc>=0);
 seen[name]=file.toString('utf8');names.push(name);
 const hash=crypto.createHash('sha256').update(file).digest('hex');
 if(bundle.manifest.files[name])assert.equal(hash,bundle.manifest.files[name].sha256);
 offset+=30+nlen+extra+len;
}
assert.equal(b.readUInt32LE(offset),0x02014b50);
assert.ok(names.includes('methods/methods.md'));
assert.ok(names.includes('provenance/manifest.json'));
assert.equal(JSON.parse(seen['source/projet_reproductible.json']).archivedResult,null);
assert.ok(seen['source/projet_reproductible.json'].includes(subject)); // PRIVATE
assert.ok(seen['provenance/analysis_complete.json'].includes(subject)); // PRIVATE
assert.ok(!seen['methods/methods_audit.json'].includes(subject)); // PUBLIC-SAFE
assert.ok(!seen['methods/methods.md'].includes(subject)); // PUBLIC-SAFE
assert.ok(!seen['tableaux/transcriptomics.csv'].includes(subject));
assert.ok(seen['tableaux/transcriptomics.csv'].includes('GENE_52'));
assert.ok(seen['tableaux/transcriptomics.csv'].includes("'=HYPERLINK"));
assert.equal(JSON.parse(seen['provenance/analysis_complete.json']).layers.transcriptomics.rows[1].feature,
  '=HYPERLINK("http://malicious.test")');
assert.throws(()=>zipEntries([{name:'../bad.txt',data:'x'}]),/Unsafe or duplicate/);
assert.throws(()=>zipEntries([{name:'same.txt',data:'a'},{name:'same.txt',data:'b'}]),/Unsafe or duplicate/);
await assert.rejects(buildResearchBundle({result,portableProject:{format:'unknown'}}),/Original study files/);
const html=methodsMarkdown(result,'en');
assert.match(html,/Methods actually run/);
assert.match(html,/Features omitted by browser QC: 4 \(6\.7 %\)/);
assert.match(html,/not executed/);

const binaryReport={...result,protocol:{...result.protocol,objective:'outcome',outcomeType:'binary'},
  predictiveOutcome:{status:'ok',metrics:{auc:.9855,brier:.2793,logLoss:.895}}};
const binaryMethods=methodsMarkdown(binaryReport,'en');
assert.match(binaryMethods,/Binary prediction — discrimination and calibration/);
assert.match(binaryMethods,/AUC \(nested CV\): 0\.9855/);
assert.match(binaryMethods,/Brier \(nested CV\): 0\.2793/);
assert.match(binaryMethods,/log-loss \(nested CV\): 0\.895/);
assert.match(binaryMethods,/High AUC does not imply reliable probabilities/);
assert.match(methodsMarkdown(binaryReport,'fr'),/AUC élevée ≠ probabilités fiables/);
assert.doesNotMatch(html,/Binary prediction — discrimination and calibration/);



const rReference={...result,referenceBackend:{
  methods:{transcriptomics_differential_deseq2:{
    method:'DESeq2',status:'ok',top:Array.from({length:100},(_,i)=>(
      {feature:'REF'+i,pvalue:0.1,adjustedP:0.3}))
  }},
  runtime:{RVersion:'4.5.2',packageVersions:{DESeq2:'1.50.2'}}
}};
const referenceMarkdown=methodsMarkdown(rReference,'fr');
assert.match(referenceMarkdown,/exhaustivité vérifiée par méthode/);
assert.match(referenceMarkdown,/export complet/);
assert.match(referenceMarkdown,/DESeq2/);
const refArchive=await buildResearchBundle({result:rReference,portableProject:portable});
assert.equal(refArchive.manifest.r_reference_full_feature_tables_exported,false);
assert.equal(refArchive.manifest.reference_method_previews.length,1);
assert.equal(refArchive.manifest.reference_method_previews[0].preview_rows,100);
assert.equal(refArchive.manifest.reference_method_previews[0].complete_reference_table_exported,false);

const specialCsv = 'feature,pvalue,padj\r\n"gene,""quoted""\nline",0.015,NA\r\nGENE2,0,0\r\n';
const specialParsed=parseReferenceCsv(specialCsv);
assert.equal(specialParsed.rows.length,2);
assert.equal(specialParsed.rows[0].feature,'gene,"quoted"\nline');
assert.equal(specialParsed.rows[0].padj,'NA');
assert.equal(specialParsed.rows[1].pvalue,'0');
assert.throws(()=>parseReferenceCsv('feature,pvalue\n"bad\n'),/Unclosed/);
assert.throws(()=>parseReferenceCsv('feature,pvalue\nG1,0.1,extra\n'),/Ragged/);
assert.throws(()=>parseReferenceCsv('feature,feature\nG1,x\n'),/duplicate/);
const realCsv = ['feature,pvalue,padj', ...Array.from({length:165},(_,i)=>
  ['GENE_'+String(i+1),0.01,0.04].join(','))].join('\n')+'\n';
const fullMethod={...result,referenceBackend:{methods:{transcriptomics_differential_deseq2:{
  method:'DESeq2',status:'ok',top:Array.from({length:100},(_,i)=>({feature:'GENE_'+i})),
  referenceTableComplete:true,referenceTableStatus:'complete',
  referenceTableRows:165,referenceTableExportedRows:165,
  referenceTableOriginalFeatures:180,referenceTableExcludedBeforeFit:15,
  pairedDesignContract:{
    pairing:'same_subject_two_conditions',design:'~ subject_id + condition',
    independent_biological_subjects:4,biological_pairs:4,model_samples:8,
    residual_degrees_of_freedom:3,
    explicit_primary_contrast:['condition','treated','control'],
    contrast_was_prespecified:true,certified:false
  },
  referenceTableCsv:realCsv,referenceTableBytes:Buffer.byteLength(realCsv)
}},runtime:{RVersion:'4.5.2'}}};
const fullBundle=await buildResearchBundle({result:fullMethod,portableProject:portable});
assert.equal(fullBundle.manifest.r_reference_full_feature_tables_exported,true);
assert.equal(fullBundle.manifest.r_reference_full_tables[0].rows,165);
assert.match(methodsMarkdown(fullMethod,'en'),/reference_r/);
assert.match(methodsMarkdown(fullMethod,'fr'),/Unités biologiques indépendantes.*4/);
assert.match(methodsMarkdown(fullMethod,'en'),/Fitted samples.*8/);
assert.equal(fullBundle.manifest.reference_method_previews[0].paired_design.biological_pairs,4);
assert.equal(fullBundle.manifest.reference_method_previews[0].paired_design.certified,false);
const nativePrimary={
  ...fullMethod,
  layers:{
    transcriptomics:{
      ...fullMethod.layers.transcriptomics,
      rows:Array.from({length:165},(_,i)=>({
        feature:'GENE_'+String(i+1),effect:0,pValue:0.01,qValue:0.04,
        numericalProvenance:'exact_native_r_reference_table'
      })),
      numericalProvenance:'exact_native_r_reference_table',
      primaryReferenceMethodKey:'transcriptomics_differential_deseq2'
    }
  }
};
const nativePrimaryMethods=methodsMarkdown(nativePrimary,'en');
assert.match(nativePrimaryMethods,/primary CSV for this layer is the verified complete native-R table/);
assert.doesNotMatch(nativePrimaryMethods,/every feature returned by the browser engine without truncation/);
assert.match(methodsMarkdown(result,'en'),/every feature returned by the browser engine without truncation/);
const failedFit=structuredClone(fullMethod);
failedFit.referenceBackend.methods.transcriptomics_differential_deseq2.status='blocked';
await assert.rejects(buildResearchBundle({result:failedFit,portableProject:portable}),
  /non-success status/);
const truncated=structuredClone(fullMethod);
truncated.referenceBackend.methods.transcriptomics_differential_deseq2.referenceTableCsv=realCsv.slice(0,-16);
await assert.rejects(buildResearchBundle({result:truncated,portableProject:portable}),/row count mismatch|byte count mismatch|Ragged scientific CSV/);

// A too-large output exists ONLY on the local R operator's machine.
// Its folder ID and SHA-256 must be disclosed, but the browser ZIP MUST
// never assert that reference_r/ contains that omitted table.
const localOnly={...result,referenceBackend:{methods:{transcriptomics_differential_deseq2:{
  method:'DESeq2',status:'ok',top:Array.from({length:100},(_,i)=>({feature:'G'+i})),
  referenceTableStatus:'too_large_archived_locally',referenceTableComplete:false,
  referenceTableRows:30000,referenceTableExportedRows:0,referenceTableCsv:null,
  referenceTableLocalArchiveId:'pmx_reference_abcdef1234',
  referenceTableLocalArchiveSha256:'a'.repeat(64)
}}}};
const localZip=await buildResearchBundle({result:localOnly,portableProject:portable});
assert.equal(localZip.manifest.r_reference_full_feature_tables_exported,false);
assert.equal(localZip.manifest.reference_method_previews[0].local_archive_id,
  'pmx_reference_abcdef1234');
assert.equal(localZip.manifest.reference_method_previews[0].local_archive_sha256,
  'a'.repeat(64));
assert.match(methodsMarkdown(localOnly,'fr'),/HORS ZIP/);
assert.match(methodsMarkdown(localOnly,'en'),/NOT included in this ZIP/);
const noArchive=structuredClone(localOnly);
const nt=noArchive.referenceBackend.methods.transcriptomics_differential_deseq2;
nt.referenceTableStatus='too_large';
nt.referenceTableLocalArchiveId=null;
nt.referenceTableLocalArchiveSha256=null;
const noArchiveZip=await buildResearchBundle({result:noArchive,portableProject:portable});
assert.equal(noArchiveZip.manifest.reference_method_previews[0].local_archive_id,null);

assert.ok(names.includes('verification/verifier_dossier.py'));
assert.ok(names.includes('verification/LISEZ-MOI.txt'));
const scratch=fs.mkdtempSync(path.join(os.tmpdir(),'pmx-research-zip-'));
try {
  const original=path.join(scratch,'original.zip');
  const script=path.join(scratch,'verifier_dossier.py');
  fs.writeFileSync(original,b);
  fs.writeFileSync(script,seen['verification/verifier_dossier.py']);
  const verified=spawnSync('python3',[script,original],{encoding:'utf8'});
  assert.equal(verified.status,0,'External Python verifier rejected original archive: '+verified.stderr);
  assert.match(verified.stdout,/INTEGRITY PASS/);
  const index=b.indexOf(Buffer.from(subject));
  assert.ok(index>=0,'Synthetic private marker not found in ZIP data for corruption test');
  const corrupted=Buffer.from(b);
  corrupted[index]^=0x01;
  const tampered=path.join(scratch,'tampered.zip');
  fs.writeFileSync(tampered,corrupted);
  const rejected=spawnSync('python3',[script,tampered],{encoding:'utf8'});
  assert.notEqual(rejected.status,0,'Tampered archive unexpectedly verified');
  assert.match(rejected.stderr,/INVALID DOSSIER|DOSSIER INVALIDE/);
} finally {
  fs.rmSync(scratch,{recursive:true,force:true});
}

console.log('MULTIOMICS RESEARCH ZIP: PASS, full 56-feature CSV, genuine ZIP structure, SHA-256, privacy boundaries, roundtrippable sources');
