// @ts-nocheck
// Source-location provenance is an analytical constraint, NOT a proxy model.
// No diagnosis, hepatic function, dose response or tissue equivalence inferred.
const VALID_BULK=new Set(['bulk','single_cell','single_nucleus','sorted_cells','cell_line','unknown']);
const VALID_COMPONENT_STATUS=new Set(['measured','estimated','not_measured','unknown']);
const VALID_RELATIVE_DOSE=new Set(['pre_dose','post_dose','at_dose','not_applicable','unknown']);
const KNOWN_LIVER=new Set(['liver','hepatic_tissue','liver_biopsy','liver_microsomes','liver_s9','hepatocytes']);
const KNOWN_BLOOD=new Set(['whole_blood','blood','pbmc','plasma','serum','peripheral_blood','extracellular_vesicles']);
function trim(v){return String(v??'').trim();}
function key(v){return trim(v).toLowerCase().replace(/[\s-]+/g,'_');}
export const SOURCE_LOCATION_COLUMNS=Object.freeze([
  'specimen_type','tissue_ontology_id','cell_type_or_bulk',
  'cell_composition_status','cell_composition_method','sample_matrix',
  'assay_platform','source_preprocessing','pharmacology_target_organ',
  'relative_to_dose','source_protocol_id'
]);
/** Validate only what the investigator explicitly provided.
 * Blank remains unknown, including bulk/single-cell status and relative dose time.
 */
export function validateSourceLocation(row){
  const result={};
  for(const col of SOURCE_LOCATION_COLUMNS)result[col]=trim(row[col]);
  if(result.cell_type_or_bulk && !VALID_BULK.has(key(result.cell_type_or_bulk)))
    throw Error('cell_type_or_bulk requires bulk/single_cell/single_nucleus/sorted_cells/cell_line/unknown');
  if(result.cell_composition_status && !VALID_COMPONENT_STATUS.has(key(result.cell_composition_status)))
    throw Error('cell_composition_status requires measured/estimated/not_measured/unknown');
  if(result.relative_to_dose && !VALID_RELATIVE_DOSE.has(key(result.relative_to_dose)))
    throw Error('relative_to_dose must be pre_dose/post_dose/at_dose/not_applicable/unknown');
  if(result.tissue_ontology_id && !/^(?:UBERON:\d{7}|CL:\d{7}|EFO:\d{7})$/.test(result.tissue_ontology_id))
    throw Error('Source tissue ontology must be an actual UBERON/CL/EFO identifier or blank; do not fabricate');
  if(result.cell_composition_method && !['measured','estimated'].includes(key(result.cell_composition_status)))
    throw Error('Cell composition method cannot claim an undeclared measurement/estimation');
  return result;
}
/** Aggregate, no individual sample IDs, no assumed biological proxy. */
export function summarizeSourceLocation(samples){
  const categories={
    n_assays:samples.length,
    n_original_tissue_declared:samples.filter(r=>trim(r.tissue)).length,
    n_specimen_type_declared:samples.filter(r=>trim(r.specimen_type)).length,
    n_cellularity_declared:samples.filter(r=>trim(r.cell_type_or_bulk)).length,
    n_source_ontology_declared:samples.filter(r=>trim(r.tissue_ontology_id)).length,
    n_cell_composition_measured_or_estimated:samples.filter(r=>['measured','estimated'].includes(key(r.cell_composition_status))).length,
    n_explicit_predose:samples.filter(r=>key(r.relative_to_dose)==='pre_dose').length,
    n_explicit_postdose:samples.filter(r=>key(r.relative_to_dose)==='post_dose').length,
    n_explicit_target_organ:samples.filter(r=>trim(r.pharmacology_target_organ)).length,
    n_cross_organ_proxy_not_validated:0,
  };
  const sourceScope=(s)=>{
    const t=key(s.tissue),p=key(s.specimen_type),m=key(s.sample_matrix);
    if(KNOWN_BLOOD.has(t)||KNOWN_BLOOD.has(p)||KNOWN_BLOOD.has(m))return 'circulating';
    if(KNOWN_LIVER.has(t)||KNOWN_LIVER.has(p))return 'hepatic';
    return 'other_or_unclassified';
  };
  for(const sample of samples){
    const source=sourceScope(sample),target=key(sample.pharmacology_target_organ);
    if(source==='circulating' && KNOWN_LIVER.has(target))
      categories.n_cross_organ_proxy_not_validated++;
  }
  return {
    ...categories,
    status: categories.n_original_tissue_declared===samples.length
      && categories.n_specimen_type_declared===samples.length
      && categories.n_cellularity_declared===samples.length
      ?'source_specimen_and_cellular_origin_declared_not_functionally_validated'
      :'incomplete_source_provenance_no_organ_proxy_claim',
    biological_source_to_target_organ_translation_validated:false,
    direct_enzyme_or_transporter_activity_inferred:false,
    predose_exposure_validated:false,
    caveat:'Whole blood, PBMC, serum, plasma, extracellular vesicles and liver are different compartments. Bulk composition or hepatic CYP functional equivalence cannot be inferred from an expression label.',
  };
}
