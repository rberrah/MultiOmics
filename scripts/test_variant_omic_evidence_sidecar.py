"""Source-method-aware observation→hypothesis tests; *synthetic* only."""
import unittest
from scripts.variant_omic_evidence_sidecar import (
    reconcile_variant_omic_evidence,VariantOmicSidecarError,
)
URL="https://pmc.ncbi.nlm.nih.gov/articles/PMC7071893/"


def obs(readout,call,unit="synthetic_subject_1",tissue="liver",landmark="predose"):
    typ=("transcriptomics" if "RNA" in readout else
         "proteomics" if "protein" in readout else
         "functional_activity" if "activity" in readout else
         "cellular_PD" if readout=="drug_cell_tolerance" else "metabolomics")
    return {"assay_id":readout,"omic":typ,"readout":readout,
            "source_qualified_relative_to_reference":call,
            "original_assay_matched_reference_validated":True,
            "assay_method":"original_source_qualified_assay_method",
            "original_feature_id":"source_"+readout,
            "original_assay_source_url":URL,
            "biological_unit_id":unit,"tissue":tissue,"collection_landmark":landmark}


def case(variant="missense",observations=None):
    return {"gene":"NUDT15","genome_assay_status":"original_variant_called",
            "variant_consequence":variant,
            "original_pairing_scope":"same_original_specimen",
            "source_variant_id":"rs116855232",
            "reference_transcript":"NM_018283.4",
            "genetic_evidence_url":URL,
            "curated_variant_function":{
                "kind":"published_reduced_function","source_url":URL},
            "observations":observations or [obs("total_RNA","reference_like")]}


class VariantMolecularSidecar(unittest.TestCase):
    def test_normal_rna_reduced_protein_is_instability_hypothesis_not_decay_proof(self):
        r=reconcile_variant_omic_evidence(case(observations=[
            obs("total_RNA","reference_like"),
            obs("protein_total","decreased"),
        ]))
        self.assertIn("HYPOTHESIS_POSTTRANSCRIPTIONAL_TRANSLATION_OR_PROTEIN_INSTABILITY",
                      r["hypotheses_to_test"])
        self.assertIn("NORMAL_RNA_WITH_LOW_PROTEIN_DOES_NOT_PROVE_PROTEIN_DECAY",
                      r["methodological_interpretation_caveats"])
        self.assertFalse(r["protein_decay_observed_directly"])
        self.assertFalse(r["causal_variant_mechanism_proven_by_this_report"])
        self.assertFalse(r["in_vivo_drug_PK_or_PD_improvement_estimated"])

    def test_normal_rna_and_protein_but_low_function_requires_catalytic_evidence(self):
        r=reconcile_variant_omic_evidence(case(observations=[
            obs("total_RNA","reference_like"),
            obs("protein_total","reference_like"),
            obs("measured_substrate_activity","decreased"),
        ]))
        self.assertIn("HYPOTHESIS_CATALYTIC_OR_LOCALIZATION_OR_SUBSTRATE_BINDING_DEFECT",
                      r["hypotheses_to_test"])
        self.assertNotIn("HYPOTHESIS_POSTTRANSCRIPTIONAL_TRANSLATION_OR_PROTEIN_INSTABILITY",
                         r["hypotheses_to_test"])
        self.assertIn("protein_localization_and_proteoform",
                      r["additional_independent_evidence_needed"])

    def test_splice_total_rna_normal_canonical_splice_low(self):
        r=reconcile_variant_omic_evidence(case("splice_acceptor",observations=[
            obs("total_RNA","reference_like"),
            obs("canonical_splice_RNA","decreased"),
            obs("aberrant_splice_RNA","increased"),
        ]))
        self.assertIn("HYPOTHESIS_ABERRANT_SPLICING_DESPITE_TOTAL_RNA",
                      r["hypotheses_to_test"])
        self.assertIn("HYPOTHESIS_CANONICAL_TO_ABERRANT_SPLICE_SHIFT",
                      r["hypotheses_to_test"])
        self.assertTrue(r["can_export_as_optional_separate_annotation_not_exchange_V1_schema_change"])

    def test_NUDT15_construct_protein_stability_not_same_as_cell_PD_or_patient(self):
        c=case(observations=[
            obs("protein_total","reference_like",tissue="HEK293T"),
            obs("drug_cell_tolerance","decreased",tissue="HEK293T"),
        ])
        c["original_pairing_scope"]="same_construct_variant_independent_assays"
        r=reconcile_variant_omic_evidence(c)
        self.assertIn("HYPOTHESIS_PRESERVED_PROTEIN_BUT_REDUCED_CELL_DRUG_TOLERANCE",
                      r["hypotheses_to_test"])
        self.assertIn("INDEPENDENT_SCREENS_SHARE_VARIANT_NOT_ORIGINAL_BIOLOGICAL_SPECIMEN",
                      r["methodological_interpretation_caveats"])
        self.assertFalse(r["assays_share_source_biological_unit_tissue_and_time"])
        self.assertFalse(r["in_vivo_drug_PK_or_PD_improvement_estimated"])
        self.assertNotIn("HYPOTHESIS_POSTTRANSCRIPTIONAL_TRANSLATION_OR_PROTEIN_INSTABILITY",
                         r["hypotheses_to_test"])

    def test_unverified_original_pairing_cannot_make_mechanism_hypothesis(self):
        x=case(observations=[obs("protein_total","decreased")])
        x.pop("original_pairing_scope")
        with self.assertRaisesRegex(VariantOmicSidecarError,"pairing scope"):
            reconcile_variant_omic_evidence(x)

    def test_stop_with_normal_total_RNA_not_automatic_normal_or_decay(self):
        r=reconcile_variant_omic_evidence(case("stop_gained",observations=[
            obs("total_RNA","reference_like"),
            obs("protein_full_length","absent"),
        ]))
        self.assertIn("HYPOTHESIS_NMD_ESCAPE_OR_ALTERNATIVE_ISOFORM_POSSIBLE",
                      r["hypotheses_to_test"])
        self.assertFalse(r["protein_decay_observed_directly"])

    def test_negative_targeted_PGx_is_not_normal_gene_function(self):
        c={"gene":"NUDT15",
           "genome_assay_status":"original_no_variant_detected_in_qualified_panel",
           "variant_consequence":"unknown",
           "original_pairing_scope":"same_original_specimen",
           "observations":[obs("protein_total","reference_like")]}
        r=reconcile_variant_omic_evidence(c)
        self.assertIn("NEGATIVE_SOURCE_PANEL_DOES_NOT_EXCLUDE_OTHER_VARIANTS_OR_ABNORMAL_FUNCTION",
                      r["methodological_interpretation_caveats"])
        with self.assertRaisesRegex(VariantOmicSidecarError,"No named functional variant"):
            reconcile_variant_omic_evidence({**c,"variant_consequence":"missense"})

    def test_no_variant_call_with_observed_low_function_is_not_a_fake_mutation(self):
        c={"gene":"NUDT15","genome_assay_status":"untyped",
           "variant_consequence":"unknown",
           "original_pairing_scope":"same_original_specimen",
           "observations":[obs("protein_total","reference_like"),
                           obs("measured_substrate_activity","decreased")]}
        r=reconcile_variant_omic_evidence(c)
        self.assertEqual(r["genomic_consequence"],"unknown")
        self.assertIn("NO_COMPLETE_GERMLINE_GENOTYPE_NO_VARIANT_FUNCTION_ASSIGNMENT",
                      r["methodological_interpretation_caveats"])

    def test_same_gene_not_same_subject_or_dose_timing_abstains(self):
        r=reconcile_variant_omic_evidence(case(observations=[
            obs("total_RNA","reference_like",unit="donor_1"),
            obs("protein_total","decreased",unit="donor_2"),
        ]))
        self.assertEqual(r["status"],"ABSTAIN_UNMATCHED_SUBJECT_TISSUE_OR_SAMPLE_LANDMARK")
        self.assertEqual(r["hypotheses_to_test"],[])
        r2=reconcile_variant_omic_evidence(case(observations=[
            obs("total_RNA","reference_like",landmark="predose"),
            obs("protein_total","decreased",landmark="postdose"),
        ]))
        self.assertEqual(r2["status"],"ABSTAIN_UNMATCHED_SUBJECT_TISSUE_OR_SAMPLE_LANDMARK")

    def test_invalid_assay_readout_and_raw_unscaled_values_fail(self):
        bad=obs("protein_total","decreased")
        bad["omic"]="transcriptomics"
        with self.assertRaisesRegex(VariantOmicSidecarError,"must match"):
            reconcile_variant_omic_evidence(case(observations=[bad]))
        bad2=obs("total_RNA","reference_like")
        bad2["original_assay_matched_reference_validated"]=False
        with self.assertRaisesRegex(VariantOmicSidecarError,"comparator"):
            reconcile_variant_omic_evidence(case(observations=[bad2]))
        with self.assertRaisesRegex(VariantOmicSidecarError,"Duplicate"):
            reconcile_variant_omic_evidence(case(observations=[
                obs("total_RNA","reference_like"),obs("total_RNA","decreased")]))

    def test_different_original_metabolites_not_automatically_enzyme_activity(self):
        r=reconcile_variant_omic_evidence(case(observations=[
            obs("endogenous_metabolite","increased"),
        ]))
        self.assertIn("ENDOGENOUS_METABOLITE_LEVEL_NOT_UNIQUE_ENZYME_OR_TRANSPORT_ACTIVITY",
                      r["methodological_interpretation_caveats"])
        self.assertIn("NO_UNIQUE_VARIANT_MECHANISM_ESTABLISHED_FROM_AVAILABLE_READOUTS",
                      r["hypotheses_to_test"])


if __name__=="__main__":
    unittest.main()
