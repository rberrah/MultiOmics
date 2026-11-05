"""Pure synthetic ISA biological-source metadata tests, no real donor data."""
import unittest
from scripts.audit_eatrisplus_v4_original_isa_material import (
    audit_original_text,SourceISAMaterialAbstention,
)

class ISAMetadataTests(unittest.TestCase):
    def test_metadata_mentions_are_not_per_assay_proof(self):
        synthetic=b"""@prefix obo: <http://example.org/> .\n
        ex:sample a ex:White_blood_cell .\n
        ex:plasma a ex:plasma .\n
        ex:method a ex:mRNA sequencing, proteomics, mass spectrometry,
        acylcarnitine metabolomics and lipidomics.\n
        ex:assay a ex:PBMC .\n"""
        r=audit_original_text(synthetic,strict=False)
        d=r["aggregate_material_and_platform_term_mentions"]
        self.assertGreater(d["explicit_plasma"],0)
        self.assertGreater(d["explicit_white_blood_cell"],0)
        self.assertGreater(d["transcriptomics"],0)
        self.assertTrue(r["literal_keyword_occurrence_not_original_assay_sample_link_proof"])
        self.assertTrue(r["individual_source_person_IDs_or_RDF_triples_never_printed"])
        self.assertNotIn("ex:sample",str(r))
        with self.assertRaisesRegex(SourceISAMaterialAbstention,"checksum"):
            audit_original_text(synthetic,strict=True)

    def test_empty_and_bad_original_source_abstain(self):
        with self.assertRaises(SourceISAMaterialAbstention):
            audit_original_text(b"",strict=False)
        with self.assertRaises(SourceISAMaterialAbstention):
            audit_original_text(b"\xff"*200,strict=False)

if __name__=="__main__":
    unittest.main()
