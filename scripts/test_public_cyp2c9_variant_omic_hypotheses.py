"""All artificial public-source schema tests, never original variant scores."""
import csv
import io
import unittest

from scripts.replay_public_cyp2c9_variant_omic_hypotheses import (
    OriginalCYP2C9Unavailable,_parse_original,replay_from_csv,
)


def artificial_data():
    out=io.StringIO()
    fields=("variant","class","position","start","end","activity_class",
            "activity_score","abundance_class","abundance_score")
    w=csv.DictWriter(out,fieldnames=fields)
    w.writeheader()
    for pos in range(1,41):
        for alt in ("C","D"):
            abnormal=alt=="C"
            w.writerow({
                "variant":f"A{pos}{alt}","class":"missense","position":pos,
                "start":"A","end":alt,
                "activity_class":"decreased" if abnormal else "wt-like",
                "activity_score":0.1 if abnormal else 1.0,
                "abundance_class":"wt-like" if abnormal else "decreased",
                "abundance_score":1.0 if abnormal else 0.1,
            })
    return out.getvalue().encode()


class CYP2C9VariantEvidenceTests(unittest.TestCase):
    def test_original_source_qualifier_fails_for_synthetic_checksum(self):
        with self.assertRaisesRegex(OriginalCYP2C9Unavailable,"SHA1"):
            _parse_original(artificial_data(),strict=True)

    def test_original_variant_building_blocks_are_not_original_human_subjects(self):
        a=replay_from_csv(artificial_data(),strict=False)
        self.assertEqual(a["n_original_source_paired_missense_constructs"],80)
        self.assertEqual(a["n_original_decreased_activity_but_wtlike_abundance"],40)
        self.assertEqual(a["n_automatically_proven_functional_or_decay_causal_mechanisms"],0)
        self.assertTrue(a["yeast_activity_and_human_cultured_cell_protein_not_same_original_biospecimen"])
        self.assertTrue(a["no_original_patient_germline_PGx_or_drug_PK"])
        self.assertEqual(a["original_evidence_hypothesis_class_counts"][
            "HYPOTHESIS_PRESERVED_ABUNDANCE_LOST_MISSENSE_FUNCTION"],40)
        self.assertNotIn("A1C",str(a))

    def test_exact_source_variant_and_both_assays_required(self):
        d=artificial_data()
        last=d.strip().splitlines()[-1]
        with self.assertRaisesRegex(OriginalCYP2C9Unavailable,"Duplicated"):
            _parse_original(d+last+b"\n",strict=False)
        with self.assertRaisesRegex(OriginalCYP2C9Unavailable,"headers"):
            _parse_original(d.replace(b"activity_class",b"other_class"),strict=False)
        with self.assertRaisesRegex(OriginalCYP2C9Unavailable,"incongruent"):
            _parse_original(d.replace(b"A1C",b"R1C"),strict=False)
        with self.assertRaises(OriginalCYP2C9Unavailable):
            _parse_original(d.replace(b"0.1",b"not_numeric"),strict=False)


if __name__=="__main__":
    unittest.main()
