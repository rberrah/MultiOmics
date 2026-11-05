"""MultiOmics research sidecar: original biofluid/tissue/assay-site provenance.

The biological organ of origin is not inferred from the detected molecular
feature's gene symbol, a plasma metabolite's pathway, or a paper title.
This module classifies the *actual collection and assay preparation*;
it does not estimate whether a feature improves clinical PK.
Do not change frozen Research-v1 or DOWNSTREAM_EXCHANGE_V1.
"""
from __future__ import annotations
from collections.abc import Mapping,Sequence

SCHEMA="omics_specimen_context_v1"
FRACTIONS={
 "bulk_liver_tissue":("liver","bulk_tissue"),
 "liver_microsomes":("liver","cell_fraction"),
 "primary_hepatocytes":("liver","isolated_native_cells"),
 "small_intestinal_tissue":("small_intestine","bulk_tissue"),
 "kidney_tissue":("kidney","bulk_tissue"),
 "whole_blood":("blood","whole_blood"),
 "PBMC":("blood","isolated_blood_cells"),
 "white_blood_cells":("blood","isolated_blood_cells"),
 "plasma":("blood","acellular_biofluid"),
 "serum":("blood","acellular_biofluid"),
 "urine":("urinary_tract","acellular_biofluid"),
 "cultured_human_cell_line":("laboratory_culture","engineered_or_immortalized_cells"),
 "engineered_yeast":("laboratory_culture","engineered_or_immortalized_cells"),
 "engineered_human_cells":("laboratory_culture","engineered_or_immortalized_cells"),
 "cell_culture_medium":("laboratory_culture","acellular_culture_medium"),
}
LAYERS={"genomics","transcriptomics","proteomics","metabolomics",
        "functional_activity","cellular_PD","drug_concentration"}
METHODS={
 "genomics":{"SNP_array","WGS","WES","targeted_genotyping","DMS_engineered_variant"},
 "transcriptomics":{"bulk_RNAseq_gene","RNAseq_exon_junction","microarray_gene_probe",
                    "qPCR_total","qPCR_isoform","miRNAseq","small_RNAseq"},
 "proteomics":{"SWATH_MS","DIA_MS","QTAP_targeted_peptide","VAMPseq_abundance",
               "ELISA","targeted_MS","untargeted_plasma_proteome"},
 "metabolomics":{"GC_TOF_untargeted","LC_MS_targeted","LC_MS_untargeted","lipidomics"},
 "functional_activity":{"CYP_shared_substrate_probe","CYP_selective_substrate_probe",
                        "Clickseq_activity","protein_turnover_pulse_chase"},
 "cellular_PD":{"drug_cell_survival","GI50_growth_inhibition"},
 "drug_concentration":{"LC_MS_plasma","LC_MS_whole_blood","LC_MS_urine"},
}
LANDMARKS={"untreated_predose","previous_dose_trough","post_first_dose",
           "post_second_dose","postdrug_response","ex_vivo_baseline",
           "ex_vivo_drug_incubation","unknown"}
PROVENANCES={"native_human_specimen","human_donor_derived_preparation",
             "engineered_variant_construct","cancer_cell_line","aggregated_reference_tissue"}
KNOWN_ORIGINS={"blood","liver","kidney","small_intestine","urinary_tract",
               "laboratory_culture"}
class SpecimenContextError(ValueError):
    pass


def qualify_specimen_records(records, *, original_matched_biological_unit_verified,
                             original_collection_time_verified):
    """Fail closed on assay-site ambiguity, don't infer tissues from biomarkers.

    Records may be from different sites if they're clearly declared as
    such and from the SAME ORIGINAL DONOR. This method does not assert
    same aliquot: that requires the same sample key, collection and prep.
    Engineered variants can be paired on the *construct* only, not donor.
    """
    if not isinstance(records,(list,tuple)) or not records or len(records)>50:
        raise SpecimenContextError("At least one original assay-site record required")
    if original_matched_biological_unit_verified is not True:
        raise SpecimenContextError("Cannot join unrelated donor and tumor lineage samples")
    if original_collection_time_verified is not True:
        raise SpecimenContextError("Original drug/pre-dose timing and processing must be qualified")
    seen=set();contexts=[]
    for r in records:
        if not isinstance(r,Mapping):
            raise SpecimenContextError("Original provenance record missing")
        assay=r.get("assay_id")
        if not isinstance(assay,str) or not assay or assay in seen:
            raise SpecimenContextError("Missing or duplicated original assay ID")
        seen.add(assay)
        fraction=r.get("collected_fraction")
        if fraction not in FRACTIONS:
            raise SpecimenContextError("Unknown original tissue, isolated cells or biofluid compartment")
        organ,format=FRACTIONS[fraction]
        if r.get("anatomical_origin")!=organ:
            raise SpecimenContextError("Claimed organ disagrees with actual collected fraction")
        layer=r.get("omic")
        if layer not in LAYERS or r.get("assay_method") not in METHODS[layer]:
            raise SpecimenContextError("Method must match the material actually measured")
        if r.get("collection_landmark") not in LANDMARKS:
            raise SpecimenContextError("Original treatment/pre-dose time not qualified")
        status=r.get("biological_provenance")
        if status not in PROVENANCES:
            raise SpecimenContextError("Original donor-derived versus engineered culture status absent")
        if status=="engineered_variant_construct" and organ!="laboratory_culture":
            raise SpecimenContextError("Engineered construct is not a directly genotyped human tissue")
        if status=="cancer_cell_line" and organ!="laboratory_culture":
            raise SpecimenContextError("Cancer cell-line measure is not a liver or plasma biopsy")
        if status=="aggregated_reference_tissue" and format!="bulk_tissue" and fraction!="whole_blood":
            raise SpecimenContextError("Tissue atlas median is not a purified patient donor cell")
        if (layer=="transcriptomics" and format in
                {"acellular_biofluid","acellular_culture_medium"} and
                r.get("cell_free_RNA_or_EV_source_qualified") is not True):
            raise SpecimenContextError("Bulk tissue RNA cannot be sourced from plasma without qualified cfRNA/EV assay")
        if (layer=="proteomics" and format=="acellular_biofluid"
                and r.get("direct_tissue_protein_claim") is True):
            raise SpecimenContextError("Plasma secretome peptides are not directly collected intracellular liver CYP protein")
        if (layer=="metabolomics" and
                r.get("asserts_single_origin_enzyme") is True):
            raise SpecimenContextError("Circulating metabolite alone cannot identify unique liver or kidney enzyme")
        sample=r.get("original_source_sample_key_scope")
        if sample not in {"same_original_specimen","same_donor_different_fraction",
                          "same_engineered_variant_different_screens",
                          "source_cohort_only_unpaired"}:
            raise SpecimenContextError("Source specimen pairing must be explicitly scoped")
        if (sample=="same_engineered_variant_different_screens"
                and status!="engineered_variant_construct"):
            raise SpecimenContextError("Same engineered variant is not germline donor pairing")
        if sample=="source_cohort_only_unpaired":
            raise SpecimenContextError("Cohort-wide study labels cannot establish same original donor or construct")
        if sample=="same_original_specimen" and status in {
                "aggregated_reference_tissue","engineered_variant_construct"}:
            raise SpecimenContextError("Reference medians or engineered protein screens are not one patient tissue aliquot")
        method=r["assay_method"]
        contexts.append({
            "assay_id":assay,
            "collected_fraction":fraction,
            "anatomical_origin":organ,
            "sample_preparation_class":format,
            "biological_provenance":status,
            "omic":layer,
            "assay_method":method,
            "collection_landmark":r["collection_landmark"],
            "original_source_sample_key_scope":sample,
            "gene_or_metabolite_identity_not_equal_anatomical_origin":True,
        })
    # Longitudinal and across-study feature labels can describe distinct
    # source tissues, but only actual independent donor proof upgrades a join.
    return {
        "schema_version":SCHEMA,
        "status":"QUALIFIED_MEASURED_SITE_AND_ASSAY_PROVENANCE_NOT_TARGET_ORGAN_FUNCTION",
        "n_original_assay_contexts":len(contexts),
        "contexts":contexts,
        "does_not_identify_tissue_of_origin_from_plasma_metabolite":True,
        "does_not_turn_blood_RNA_into_liver_gene_expression":True,
        "does_not_upgrade_variant_construct_to_patient_genotype":True,
        "in_vivo_PK_omic_superiority_inferred":False,
        "original_numeric_patient_omics_or_PK_released":False,
    }


def reconcile_variant_with_original_specimens(case, specimen_records, *,
                                             original_matched_biological_unit_verified=True,
                                             original_collection_time_verified=True):
    """Qualify complete assay↔physical specimen mapping before hypotheses.

    This is an *optional* v1 companion to the existing frozen
    DOWNSTREAM_EXCHANGE_V1, NOT a modification to its data schema.
    """
    from scripts.variant_omic_evidence_sidecar import (
        reconcile_variant_omic_evidence, VariantOmicSidecarError,
    )
    source=qualify_specimen_records(
        specimen_records,
        original_matched_biological_unit_verified=original_matched_biological_unit_verified,
        original_collection_time_verified=original_collection_time_verified)
    if not isinstance(case,Mapping) or not isinstance(case.get("observations"),(tuple,list)):
        raise SpecimenContextError("Original variant/assay observations absent")
    a={r.get("assay_id"):r for r in case["observations"] if isinstance(r,Mapping)}
    b={r["assay_id"]:r for r in specimen_records}
    if len(a)!=len(case["observations"]) or set(a)!=set(b):
        raise SpecimenContextError("Every molecular source assay must be mapped to one physical sample")
    for k in a:
        if a[k].get("omic")!=b[k].get("omic"):
            raise SpecimenContextError("Original sample preparation mismatched across assay and biological evidence")
    original=reconcile_variant_omic_evidence(case)
    if original.get("status")!="QUALIFIED_ORIGINAL_ASSAY_DESCRIPTIVE_HYPOTHESES_NOT_CAUSAL_PROOF":
        raise SpecimenContextError("Source assay context cannot rescue an abstained biological evidence join")
    actual_pairing=original.get("source_pairing_scope")
    if actual_pairing=="same_construct_variant_independent_assays":
        if any(x["original_source_sample_key_scope"]!="same_engineered_variant_different_screens"
               for x in source["contexts"]):
            raise SpecimenContextError("Cell variants require construct-level source pairing, never same patient")
    elif actual_pairing=="same_original_specimen":
        if any(x["original_source_sample_key_scope"] not in {
                "same_original_specimen","same_donor_different_fraction"}
               for x in source["contexts"]):
            raise SpecimenContextError("Patient tissues cannot be paired by construct or cohort label")
        # Different organs/samples require the original donor match in the
        # specimen contract. The legacy biological engine declines
        # comparing those fractions as a within-tissue discordance.
    else:
        raise SpecimenContextError("Unknown biological source pairing type")
    return {
        **original,
        "separate_specimen_sidecar_schema":SCHEMA,
        "actual_original_assay_specimens":source["contexts"],
        "specimen_provenance_qualified_without_inferring_function":True,
        "not_a_clinical_PK_or_Bayesian_function_prior":True,
    }
