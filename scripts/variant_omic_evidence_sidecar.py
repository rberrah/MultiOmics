"""Optional MultiOmics research sidecar: genetic variant ↔ molecular discordance.

No modification to the Research v1 scientific core, DESeq2/MOFA, or
DOWNSTREAM_EXCHANGE_V1 schema. This report is a *qualitative* independent
annotation, never a prediction of real patient PK or a mechanism proof.

Assumes externally source-verified original assay calls, a qualified
comparison to an appropriate reference (not comparing raw RNA counts to
proteomic intensity) and matching biological subject/tissue/time.

Clinical or pharmacological causal conclusions require allele phasing,
isogenic perturbation, validated proteoform/turnover/enzymatic assays,
patient PK/PD validation and externally qualified evidence.
"""
from __future__ import annotations
from collections.abc import Mapping,Sequence
from urllib.parse import urlparse

class VariantOmicSidecarError(ValueError):
    pass

CONSEQUENCES={
    "missense","synonymous","stop_gained","frameshift","splice_donor",
    "splice_acceptor","splice_region","intron","regulatory","unknown",
}
GENOTYPE={"original_variant_called","original_no_variant_detected_in_qualified_panel",
          "untyped","partial","engineered_protein_variant_construct"}
OMICS={"transcriptomics","proteomics","metabolomics","functional_activity","cellular_PD"}
CALLED={"reference_like","decreased","increased","absent","unknown"}
READOUTS={
    "total_RNA":"transcriptomics",
    "canonical_splice_RNA":"transcriptomics",
    "aberrant_splice_RNA":"transcriptomics",
    "protein_total":"proteomics",
    "protein_full_length":"proteomics",
    "protein_variant_specific":"proteomics",
    "endogenous_metabolite":"metabolomics",
    "measured_substrate_activity":"functional_activity",
    "drug_cell_tolerance":"cellular_PD",
}
FUNCTION={"measured_substrate_activity"}
OBSERVED_SIGNAL_PROOF="OBSERVED_ORIGINAL_SOURCE_ASSAY_CONTRAST_NOT_BIOCHEMICAL_CAUSAL_PROOF"


def _short(value,label):
    if not isinstance(value,str) or not value.strip() or len(value)>512:
        raise VariantOmicSidecarError(label+" requires explicit string")
    return value.strip()


def _https(value,label):
    v=_short(value,label)
    if urlparse(v).scheme!="https" or not urlparse(v).hostname:
        raise VariantOmicSidecarError(label+" needs source HTTPs evidence")
    return v


def reconcile_variant_omic_evidence(case):
    """Validate original qualitative readouts and output *testable* mechanisms.

    Each readout is already normalized/classified by the source assay to
    an independent control. This engine does not calculate those classes.
    The common biological subject, organ/tissue and landmark are required
    to distinguish within-subject signals from unrelated cohort proxies.
    """
    if not isinstance(case,Mapping):
        raise VariantOmicSidecarError("Original case metadata required")
    gene=_short(case.get("gene"),"source HGNC gene").upper()
    status=case.get("genome_assay_status")
    if status not in GENOTYPE:
        raise VariantOmicSidecarError("Genotype coverage must distinguish no variants from no assay")
    consequence=case.get("variant_consequence","unknown")
    if consequence not in CONSEQUENCES:
        raise VariantOmicSidecarError("Exact annotated genetic consequence required")
    if status in {"original_variant_called","engineered_protein_variant_construct"}:
        _short(case.get("source_variant_id"),"original source variant or construct ID")
        if status=="original_variant_called":
            _short(case.get("reference_transcript"),"source versioned transcript")
        else:
            _short(case.get("reference_protein_accession"),"original engineered protein reference accession")
        _https(case.get("genetic_evidence_url"),"original variant / construct source")
    else:
        if consequence!="unknown" or case.get("source_variant_id"):
            raise VariantOmicSidecarError("No named functional variant when untyped or only negative assayed sites")
    evidence=case.get("curated_variant_function")
    if evidence is not None:
        if not isinstance(evidence,Mapping) or evidence.get("kind") not in {
            "published_no_function","published_reduced_function",
            "published_normal_function","not_established"}:
            raise VariantOmicSidecarError("Published allele-function assertion requires externally curated proof")
        _https(evidence.get("source_url"),"published variant function")
    pairing=case.get("original_pairing_scope")
    if pairing not in {"same_original_specimen","same_construct_variant_independent_assays"}:
        raise VariantOmicSidecarError("Original pairing scope must distinguish specimen from independent construct assays")
    if (pairing=="same_construct_variant_independent_assays" and
            status not in {"original_variant_called","engineered_protein_variant_construct"}):
        raise VariantOmicSidecarError("Independent constructs require source-qualified identical variant identity")
    if pairing=="same_original_specimen" and status=="engineered_protein_variant_construct":
        raise VariantOmicSidecarError("Engineered DMS construct cannot be a directly genotyped human biological specimen")
    observations=case.get("observations")
    if (not isinstance(observations,(list,tuple)) or len(observations)<1
            or len(observations)>100):
        raise VariantOmicSidecarError("Nonempty original per-assay evidence descriptors required")
    seen=set()
    relevant={}
    contexts=set()
    conflicting=[]
    for entry in observations:
        if not isinstance(entry,Mapping):
            raise VariantOmicSidecarError("Original molecular assay record must be a mapping")
        aid=_short(entry.get("assay_id"),"original assay ID")
        if aid in seen:
            raise VariantOmicSidecarError("Duplicate original assay ID")
        seen.add(aid)
        target=entry.get("readout")
        if target not in READOUTS or entry.get("omic")!=READOUTS[target]:
            raise VariantOmicSidecarError("Assay method must match what the molecular readout truly measures")
        call=entry.get("source_qualified_relative_to_reference")
        if call not in CALLED:
            raise VariantOmicSidecarError("Measured assay effect must be source-qualified; no raw counts ratio")
        if entry.get("original_assay_matched_reference_validated") is not True:
            raise VariantOmicSidecarError("Reported qualitative effect lacks a prequalified comparator")
        _short(entry.get("assay_method"),"assay technology")
        _short(entry.get("original_feature_id"),"original feature")
        _https(entry.get("original_assay_source_url"),"original molecular assay source")
        if pairing=="same_construct_variant_independent_assays":
            # Two independent systems may assay *exactly the same variant*
            # in different tissues, hosts and post-drug time windows.
            # Do NOT pretend they are one patient's original specimen.
            matched_construct=_short(
                entry.get("source_variant_construct_id"),
                "independently source-verified construct variant identity")
            if matched_construct!=case["source_variant_id"]:
                raise VariantOmicSidecarError(
                    "Cross-assay variant identity mismatch: ABSTAIN instead of arbitrary gene-level join")
        contexts.add((
            _short(entry.get("biological_unit_id"),"biological subject"),
            _short(entry.get("tissue"),"measured tissue"),
            _short(entry.get("collection_landmark"),"pre/post-dose timing"),
        ))
        if target in relevant:
            conflicting.append(target)
        relevant[target]=call
    if len(contexts)>1 and pairing=="same_original_specimen":
        return {
            "status":"ABSTAIN_UNMATCHED_SUBJECT_TISSUE_OR_SAMPLE_LANDMARK",
            "gene":gene,
            "source_qualified_observation_count":len(observations),
            "hypotheses_to_test":[],
            "reason":"ASSAYS_FROM_DIFFERENT_BIOLOGICAL_CONTEXTS_CANNOT_ESTABLISH_WITHIN_SUBJECT_DISCORDANCE",
            "no_real_patient_PK_PD_function_or_causal_effect_claim":True,
        }
    if conflicting:
        raise VariantOmicSidecarError("Duplicated readout needs explicit replicate aggregation and QC")
    notes=[];required=[];qualifiers=[]
    rna=relevant.get("total_RNA")
    canon=relevant.get("canonical_splice_RNA")
    aberrant=relevant.get("aberrant_splice_RNA")
    protein=relevant.get("protein_full_length",relevant.get("protein_total"))
    activity=relevant.get("measured_substrate_activity")
    cellular_tolerance=relevant.get("drug_cell_tolerance")
    if (rna in {"reference_like","increased"} and protein in {"decreased","absent"}):
        notes.append("HYPOTHESIS_POSTTRANSCRIPTIONAL_TRANSLATION_OR_PROTEIN_INSTABILITY")
        required+=["pulse_chase_protein_turnover_or_proteasome_rescue","variant_specific_protein_mapping"]
        qualifiers.append("NORMAL_RNA_WITH_LOW_PROTEIN_DOES_NOT_PROVE_PROTEIN_DECAY")
    if (rna in {"reference_like","increased"} and protein in {"reference_like","increased"}
            and activity in {"decreased","absent"}):
        notes.append("HYPOTHESIS_CATALYTIC_OR_LOCALIZATION_OR_SUBSTRATE_BINDING_DEFECT")
        required+=["independent_substrate_specific_function","protein_localization_and_proteoform"]
        qualifiers.append("NORMAL_TOTAL_PROTEIN_CAN_BE_NONFUNCTIONAL")
    if (rna in {"reference_like","increased"} and canon in {"decreased","absent"}):
        notes.append("HYPOTHESIS_ABERRANT_SPLICING_DESPITE_TOTAL_RNA")
        required+=["splice_junction_RNA_validation","translated_isoform_specific_peptide"]
        qualifiers.append("TOTAL_RNA_NOT_EQUIVALENT_TO_FULL_LENGTH_TRANSCRIPT")
    if aberrant in {"increased"} and canon in {"decreased","absent"}:
        notes.append("HYPOTHESIS_CANONICAL_TO_ABERRANT_SPLICE_SHIFT")
        required+=["long_read_or_RTPCR_exon_junction","functional_activity"]
    if (protein in {"reference_like","increased"} and activity in {"decreased","absent"}
            and consequence=="missense"):
        if "HYPOTHESIS_CATALYTIC_OR_LOCALIZATION_OR_SUBSTRATE_BINDING_DEFECT" not in notes:
            notes.append("HYPOTHESIS_PRESERVED_ABUNDANCE_LOST_MISSENSE_FUNCTION")
            required+=["orthogonal_substrate_specific_function","allele_phase_and_proteoform"]
    if (protein in {"reference_like","increased"} and
            cellular_tolerance in {"decreased","absent"}):
        notes.append("HYPOTHESIS_PRESERVED_PROTEIN_BUT_REDUCED_CELL_DRUG_TOLERANCE")
        required+=["orthogonal_substrate_specific_enzyme_activity",
                   "cell_viability_and_drug_uptake_controls"]
        qualifiers.append("CELL_PD_IS_NOT_ENZYME_ACTIVITY_OR_CLINICAL_TOXICITY")
    if (protein in {"decreased","absent"} and
            cellular_tolerance in {"reference_like","increased"}):
        notes.append("HYPOTHESIS_LOW_PROTEIN_BUT_PRESERVED_CELL_DRUG_TOLERANCE")
        required+=["orthogonal_variant_specific_hydrolase_activity",
                   "independent_viability_and_drug_trafficking_controls"]
        qualifiers.append("LOW_PROTEIN_DOES_NOT_ESTABLISH_LOST_CELL_DRUG_PROTECTION")
    if consequence in {"stop_gained","frameshift"} and rna in {"reference_like","increased"}:
        notes.append("HYPOTHESIS_NMD_ESCAPE_OR_ALTERNATIVE_ISOFORM_POSSIBLE")
        required+=["original_transcript_specific_stop_position","allelic_aware_RNA_isoform"]
        qualifiers.append("CANONICAL_NMD_NOT_PROVEN_OR_EXCLUDED_FROM_RNA_TOTAL")
    if consequence in {"splice_donor","splice_acceptor","splice_region"} and canon is None:
        notes.append("UNRESOLVED_SPLICE_CONSEQUENCE_REQUIRES_ISOFORM_ASSAY")
        required+=["RNA_splice_junction","long_read_isoform"]
    if status=="original_no_variant_detected_in_qualified_panel":
        qualifiers.append("NEGATIVE_SOURCE_PANEL_DOES_NOT_EXCLUDE_OTHER_VARIANTS_OR_ABNORMAL_FUNCTION")
    if status in {"untyped","partial"}:
        qualifiers.append("NO_COMPLETE_GERMLINE_GENOTYPE_NO_VARIANT_FUNCTION_ASSIGNMENT")
    if evidence and evidence.get("kind") in {"published_no_function","published_reduced_function"}:
        qualifiers.append("PUBLISHED_ALLELE_FUNCTION_IS_NOT_INDEPENDENT_PATIENT_ACTIVITY_MEASUREMENT")
    if activity in {"decreased","absent"} and protein not in {"decreased","absent"}:
        required.append("assay_substrate_specificity_and_batch_controls")
    if cellular_tolerance is not None:
        qualifiers.append("CELL_TOLERANCE_RESPONSE_NOT_PATIENT_TOXICITY_OR_DRUG_CLEARANCE")
    if pairing=="same_construct_variant_independent_assays":
        qualifiers.append("INDEPENDENT_SCREENS_SHARE_VARIANT_NOT_ORIGINAL_BIOLOGICAL_SPECIMEN")
    if "endogenous_metabolite" in relevant:
        qualifiers.append("ENDOGENOUS_METABOLITE_LEVEL_NOT_UNIQUE_ENZYME_OR_TRANSPORT_ACTIVITY")
    if not notes:
        notes.append("NO_UNIQUE_VARIANT_MECHANISM_ESTABLISHED_FROM_AVAILABLE_READOUTS")
    return {
        "schema_version":"variant_omic_evidence_sidecar_v1",
        "status":"QUALIFIED_ORIGINAL_ASSAY_DESCRIPTIVE_HYPOTHESES_NOT_CAUSAL_PROOF",
        "gene":gene,"genomic_consequence":consequence,
        "genotype_coverage_state":status,
        "original_person_germline_genotype_measured":
            status=="original_variant_called",
        "source_variant_is_engineered_construct_not_human_genotype":
            status=="engineered_protein_variant_construct",
        "original_per_assay_signals_measured_not_inferred":dict(sorted(relevant.items())),
        "assays_share_source_biological_unit_tissue_and_time":(pairing=="same_original_specimen" and len(contexts)==1),
        "source_pairing_scope":pairing,
        "same_original_variant_identity_across_independent_assays_verified":
            pairing=="same_construct_variant_independent_assays",
        "do_not_use_postdrug_cell_PD_as_predose_pk_covariate":
            pairing=="same_construct_variant_independent_assays",
        "same_construct_variant_does_not_equal_same_patient":True,
        "hypotheses_to_test":sorted(set(notes)),
        "additional_independent_evidence_needed":sorted(set(required)),
        "methodological_interpretation_caveats":sorted(set(qualifiers)),
        "published_variant_effect_if_curated":evidence["kind"] if evidence else None,
        "source_readout_qualifications":OBSERVED_SIGNAL_PROOF,
        "causal_variant_mechanism_proven_by_this_report":False,
        "protein_decay_observed_directly":False,
        "in_vivo_drug_PK_or_PD_improvement_estimated":False,
        "can_export_as_optional_separate_annotation_not_exchange_V1_schema_change":True,
        "does_not_modify_original_omics_or_per_patient_measurements":True,
    }
