#!/usr/bin/env python3
"""Synthetic-only demonstration. Not original human Omics or PK/PD rows."""
import json
from scripts.variant_omic_evidence_sidecar import reconcile_variant_omic_evidence

SOURCE="https://pmc.ncbi.nlm.nih.gov/articles/PMC7071893/"
def readout(name,kind,omic,signal):
    return {"assay_id":name,"readout":kind,"omic":omic,
            "source_qualified_relative_to_reference":signal,
            "original_assay_matched_reference_validated":True,
            "assay_method":"SYNTHETIC_REFERENCE_QUALIFIED",
            "original_feature_id":"SYNTH_"+name,
            "original_assay_source_url":SOURCE,
            "biological_unit_id":"SYNTHETIC_MUTATION_CONSTRUCT",
            "tissue":"synthetic_cell_context","collection_landmark":"before_drug"}

def case(consequence,observations,scope="same_original_specimen"):
    return {"gene":"NUDT15",
            "genome_assay_status":"original_variant_called",
            "variant_consequence":consequence,
            "source_variant_id":"HYPOTHETICAL_VARIANT_NOT_GENOTYPED",
            "reference_transcript":"NM_018283.4",
            "genetic_evidence_url":SOURCE,
            "original_pairing_scope":scope,
            "observations":observations}

def run():
    rn=readout("total_rna","total_RNA","transcriptomics","reference_like")
    r_low=readout("protein_low","protein_total","proteomics","decreased")
    p_high=readout("protein_hi","protein_total","proteomics","reference_like")
    fun=readout("enzyme_low","measured_substrate_activity","functional_activity","decreased")
    sensitive=readout("drug_low","drug_cell_tolerance","cellular_PD","decreased")
    splice=readout("canonical","canonical_splice_RNA","transcriptomics","decreased")
    return {
        "source_status":"SYNTHETIC_QUALITATIVE_EXAMPLES_NOT_ORIGINAL_SAMPLES",
        "rna_normal_protein_low_no_decay_proof":reconcile_variant_omic_evidence(
            case("missense",[rn,r_low])),
        "rna_and_protein_normal_function_low":reconcile_variant_omic_evidence(
            case("missense",[rn,p_high,fun])),
        "gene_total_RNA_normal_canonical_splice_low":reconcile_variant_omic_evidence(
            case("splice_acceptor",[rn,splice])),
        "nudt15_variant_matched_construct_stability_vs_cell_PD":reconcile_variant_omic_evidence(
            case("missense",[p_high,sensitive],
                 scope="same_construct_variant_independent_assays")),
    }
if __name__=="__main__":
    print(json.dumps(run(),indent=2,sort_keys=True))
