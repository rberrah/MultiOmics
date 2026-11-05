"""All tests below are synthetic and carry NO genuine MaveDB variant values."""
import csv
import io
import unittest

from scripts.replay_public_nudt15_variant_omic_hypotheses import (
    URN,NUOriginalSourceError,
    _source_valid_amino_acid_substitution,
    replay_from_bytes, original_source_parse,
)


def fake_source():
    stream=io.StringIO()
    writer=csv.DictWriter(stream,fieldnames=[
        "accession","hgvs_nt","hgvs_splice","hgvs_pro",
        "score","abundance_score","sensitivity_score"])
    writer.writeheader()
    for pos in range(1,41):
        for variant in ("Arg","Gly"):
            # Two independently assayed phenotypes indexed by the
            # exact same original protein construct. One intentionally
            # stable-but-sensitive / one low-protein-but-tolerant.
            stable=variant=="Arg"
            writer.writerow({
                "accession":f"{URN}#{pos}_{variant}",
                "hgvs_nt":"","hgvs_splice":"",
                "hgvs_pro":f"p.(Ala{pos}{variant})",
                "score":0.45,
                "abundance_score":.9 if stable else .05,
                "sensitivity_score":.05 if stable else .9,
            })
    return stream.getvalue().encode()


class NURealStructureSyntheticTests(unittest.TestCase):
    def test_complete_length_164_and_unknown_variants(self):
        self.assertTrue(_source_valid_amino_acid_substitution("p.(Ala164Cys)"))
        self.assertTrue(_source_valid_amino_acid_substitution("p.Arg139Cys"))
        for value in ("p.Arg165Cys","p.(Arg12Arg)","p.=",
                      "p.Arg12Ter","c.12A>G","none"):
            self.assertFalse(_source_valid_amino_acid_substitution(value))

    def test_original_variant_paired_independent_assays_yield_two_directions(self):
        result=replay_from_bytes(fake_source(),strict_original=False)
        self.assertEqual(result["n_paired_original_missense_constructs"],80)
        self.assertEqual(result["n_preserved_protein_abundance_but_low_cell_drug_tolerance"],40)
        self.assertEqual(result["n_low_protein_abundance_but_preserved_cell_drug_tolerance"],40)
        self.assertEqual(result["n_causal_mechanism_claims"],0)
        self.assertTrue(result["source_same_construct_not_same_original_patient"])
        self.assertTrue(result["no_genetic_PGx_diplotypes_from_engineered_DMS_constructs"])
        self.assertTrue(result["postdrug_cell_PD_not_eligible_as_predose_PK_feature"])
        self.assertNotIn("Ala1Arg",str(result))
        counts=result["hypothesis_counts_aggregated_only"]
        self.assertEqual(counts["HYPOTHESIS_PRESERVED_PROTEIN_BUT_REDUCED_CELL_DRUG_TOLERANCE"],40)
        self.assertEqual(counts["HYPOTHESIS_LOW_PROTEIN_BUT_PRESERVED_CELL_DRUG_TOLERANCE"],40)

    def test_source_sha256_not_fabricated_by_toy_fixture(self):
        with self.assertRaisesRegex(NUOriginalSourceError,"SHA256"):
            original_source_parse(fake_source(),strict=True)
        rows,source=original_source_parse(fake_source(),strict=False)
        self.assertEqual(source["source_rows"],80)
        self.assertEqual(len(rows),80)

    def test_variant_id_duplicate_or_corrupted_source_fails_closed(self):
        raw=fake_source()
        last=raw.strip().splitlines()[-1]
        with self.assertRaisesRegex(NUOriginalSourceError,"duplicated"):
            original_source_parse(raw+last+b"\n",strict=False)
        with self.assertRaisesRegex(NUOriginalSourceError,"fields"):
            original_source_parse(raw.replace(b"sensitivity_score",b"other"),
                                  strict=False)


if __name__=="__main__":
    unittest.main()
