"""Synthetic-only tissue/sampling-site and assay-origin scientific safety checks."""
import unittest
from scripts.specimen_context_sidecar import (
    qualify_specimen_records,reconcile_variant_with_original_specimens,SpecimenContextError,
)
from scripts.test_variant_omic_evidence_sidecar import obs,case


def sample(name="total_RNA",layer="transcriptomics",fraction="bulk_liver_tissue",
           origin="liver",method="bulk_RNAseq_gene",
           landmark="untreated_predose",study="native_human_specimen",
           scope="same_original_specimen",**kwargs):
    return {
        "assay_id":name,"omic":layer,"collected_fraction":fraction,
        "anatomical_origin":origin,"assay_method":method,
        "collection_landmark":landmark,
        "biological_provenance":study,
        "original_source_sample_key_scope":scope,**kwargs,
    }


def checked(rows):
    return qualify_specimen_records(rows,
        original_matched_biological_unit_verified=True,
        original_collection_time_verified=True)


class SpecimenSourceTests(unittest.TestCase):
    def test_hlc_direct_liver_RNA_vs_microsomal_activity_distinct_fractions(self):
        x=checked([
            sample("rna",fraction="bulk_liver_tissue",method="microarray_gene_probe"),
            sample("activity",layer="functional_activity",
                   fraction="liver_microsomes",method="CYP_shared_substrate_probe",
                   study="human_donor_derived_preparation",scope="same_donor_different_fraction"),
        ])
        self.assertEqual([z["collected_fraction"] for z in x["contexts"]],
                         ["bulk_liver_tissue","liver_microsomes"])
        self.assertTrue(x["does_not_turn_blood_RNA_into_liver_gene_expression"])

    def test_plasma_metabolomics_rotroff_predose_and_explicit_postdose(self):
        x=checked([
            sample("rotroff_A",layer="metabolomics",fraction="plasma",
                   origin="blood",method="GC_TOF_untargeted"),
            sample("rotroff_B",layer="metabolomics",fraction="plasma",
                   origin="blood",method="GC_TOF_untargeted",
                   landmark="post_first_dose",scope="same_donor_different_fraction"),
        ])
        self.assertEqual(x["contexts"][0]["anatomical_origin"],"blood")
        self.assertEqual(x["contexts"][1]["collection_landmark"],"post_first_dose")
        self.assertTrue(x["does_not_identify_tissue_of_origin_from_plasma_metabolite"])

    def test_blood_immuncell_RNA_never_called_liver_CYP_transcript(self):
        a=checked([sample("RNA",fraction="PBMC",origin="blood",
                          method="bulk_RNAseq_gene")])
        self.assertEqual(a["contexts"][0]["anatomical_origin"],"blood")
        with self.assertRaisesRegex(SpecimenContextError,"disagrees"):
            checked([sample("RNA",fraction="PBMC",origin="liver")])
        with self.assertRaisesRegex(SpecimenContextError,"cfRNA"):
            checked([sample("RNA",fraction="plasma",origin="blood")])
        checked([sample("cfRNA",fraction="plasma",origin="blood",
                        cell_free_RNA_or_EV_source_qualified=True)])

    def test_plasma_protein_not_intracellular_cyp_and_metabolite_not_specific_enzyme(self):
        with self.assertRaisesRegex(SpecimenContextError,"intracellular"):
            checked([sample("protein",layer="proteomics",fraction="plasma",
                            origin="blood",method="untargeted_plasma_proteome",
                            direct_tissue_protein_claim=True)])
        with self.assertRaisesRegex(SpecimenContextError,"unique"):
            checked([sample("met",layer="metabolomics",fraction="plasma",
                            origin="blood",method="LC_MS_targeted",
                            asserts_single_origin_enzyme=True)])

    def test_original_cell_line_variant_is_not_same_human_biospecimen(self):
        rows=[
            sample("protein_total","proteomics","engineered_human_cells",
                   "laboratory_culture","VAMPseq_abundance",
                   "ex_vivo_baseline","engineered_variant_construct",
                   "same_engineered_variant_different_screens"),
            sample("drug_cell_tolerance","cellular_PD","engineered_human_cells",
                   "laboratory_culture","drug_cell_survival",
                   "postdrug_response","engineered_variant_construct",
                   "same_engineered_variant_different_screens"),
        ]
        c=case(observations=[
            obs("protein_total","reference_like",landmark="before_drug"),
            obs("drug_cell_tolerance","decreased",landmark="post_drug"),
        ])
        c["genome_assay_status"]="engineered_protein_variant_construct"
        c["reference_protein_accession"]="UniProt:Q9NV35-1"
        c.pop("reference_transcript")
        c["original_pairing_scope"]="same_construct_variant_independent_assays"
        r=reconcile_variant_with_original_specimens(c,rows)
        self.assertEqual(r["actual_original_assay_specimens"][1]["collection_landmark"],"postdrug_response")
        self.assertFalse(r["in_vivo_drug_PK_or_PD_improvement_estimated"])
        self.assertTrue(r["specimen_provenance_qualified_without_inferring_function"])
        with self.assertRaisesRegex(SpecimenContextError,"Cell variants"):
            rows[0]["original_source_sample_key_scope"]="same_original_specimen"
            reconcile_variant_with_original_specimens(c,rows)

    def test_missing_record_unverified_cohort_and_ambiguous_assay_abstain(self):
        with self.assertRaisesRegex(SpecimenContextError,"Every molecular"):
            reconcile_variant_with_original_specimens(
                case(observations=[obs("total_RNA","reference_like")]),
                [sample("wrong")])
        with self.assertRaisesRegex(SpecimenContextError,"Cohort-wide"):
            checked([sample(scope="source_cohort_only_unpaired")])
        with self.assertRaisesRegex(SpecimenContextError,"must be qualified"):
            qualify_specimen_records([sample()],
                original_matched_biological_unit_verified=True,
                original_collection_time_verified=False)
        with self.assertRaisesRegex(SpecimenContextError,"Method"):
            checked([sample(method="QTAP_targeted_peptide")])

    def test_healthy_plasma_multiomic_not_claim_same_fraction_as_white_cells(self):
        a=checked([
            sample("WBC",layer="transcriptomics",fraction="white_blood_cells",
                   origin="blood",method="bulk_RNAseq_gene",
                   scope="same_donor_different_fraction"),
            sample("plasma_proteome",layer="proteomics",fraction="plasma",
                   origin="blood",method="untargeted_plasma_proteome",
                   scope="same_donor_different_fraction"),
        ])
        self.assertNotEqual(a["contexts"][0]["sample_preparation_class"],
                            a["contexts"][1]["sample_preparation_class"])


if __name__=="__main__":
    unittest.main()
