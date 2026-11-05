"""Synthetic HDF5 MuData schema tests; no source EATRIS person values."""
import tempfile
import unittest
from pathlib import Path
from scripts.audit_eatrisplus_v4_original_h5mu import (
    audit_h5mu,EATRISSourceAbstention,_modality_origin_guess,
)


def _fixture(path, malformed=None):
    import h5py
    with h5py.File(path,"w") as f:
        f.create_group("obs").create_dataset(
            "_index",data=[f"SYNTHETIC_P{i}".encode() for i in range(8)])
        f.create_group("var").create_dataset("_index",data=[b"unused"])
        mods=f.create_group("mod")
        for mode,people in (
            ("mRNA",range(7)),("plasma_proteomics",range(2,8)),
            ("targeted_metabolomics",range(3,8)),
        ):
            m=mods.create_group(mode)
            if malformed=="alien_id" and mode=="plasma_proteomics":
                ids=[b"SYNTHETIC_NOT_IN_THE_SOURCE"]+[f"SYNTHETIC_P{i}".encode() for i in range(3,8)]
            elif malformed=="duplicate" and mode=="mRNA":
                ids=[b"SYNTHETIC_P0"]*7
            else:ids=[f"SYNTHETIC_P{i}".encode() for i in people]
            m.create_group("obs").create_dataset("_index",data=ids)
            m.create_group("var").create_dataset(
                "_index",data=[f"VALID_FEATURE_{j}".encode() for j in range(5)])
            if mode=="mRNA":
                # sparse HDF5 CSR form, with shape attribute
                x=m.create_group("X")
                x.attrs["shape"]=[len(ids),5]
                x.create_dataset("indptr",data=[0]*(len(ids)+1))
                x.create_dataset("data",data=[],dtype="f")
                x.create_dataset("indices",data=[],dtype="i")
            elif malformed=="missing_X" and mode=="plasma_proteomics":
                continue
            else:
                m.create_dataset("X",shape=(len(ids),5),dtype="f",
                                 fillvalue=2.0)


class RealSourceStructureMock(unittest.TestCase):
    def test_mudata_v4_true_original_sample_key_matching_is_permutation_invariant(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"fake.h5mu"
            _fixture(p)
            r=audit_h5mu(p,verify_original_checksum=False)
            self.assertEqual(r["status"],
                "ORIGINAL_ZENODO_V4_MUDATA_IDENTITY_METADATA_AUDIT_NOT_CLINICAL_PK")
            self.assertEqual(r["n_original_global_observation_keys"],8)
            self.assertEqual(r["n_omics_assay_modalities"],3)
            self.assertEqual(r["n_exact_original_keys_shared_across_all_modalities"],4)
            self.assertEqual(r["pairwise_exact_original_index_overlap"]
                             ["mRNA + plasma_proteomics"]
                             ["n_exact_original_identical_sample_keys"],5)
            self.assertFalse(r["original_human_drug_PK_endpoint_present"])
            self.assertFalse(r["joint_PGx_plus_omics_superiority_can_be_tested"])
            self.assertNotIn("SYNTHETIC_P",str(r))
            self.assertEqual(r["per_original_modality"]["mRNA"]
                             ["source_matrix_representation"],"HDF5_sparse")

    def test_checksum_and_source_identity_fail_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"fake.h5mu"
            _fixture(p)
            with self.assertRaisesRegex(EATRISSourceAbstention,"md5"):
                audit_h5mu(p,verify_original_checksum=True)
            for bug,expected in (
                ("alien_id","absent from true MuData global index"),
                ("duplicate","duplicated"),
                ("missing_X","matrix missing"),
            ):
                _fixture(p,malformed=bug)
                with self.subTest(bug=bug),self.assertRaisesRegex(
                    EATRISSourceAbstention,expected):
                    audit_h5mu(p,verify_original_checksum=False)

    def test_original_method_name_only_is_NOT_actual_organ_measurement(self):
        self.assertIn("REQUIRES_SOURCE_METHOD_MATCH",
                      _modality_origin_guess("mRNA"))
        self.assertIn("REQUIRES_SOURCE_METHOD_MATCH",
                      _modality_origin_guess("plasma_proteomics"))
        self.assertEqual(_modality_origin_guess("arbitrary_key"),
                         "UNRESOLVED_ORIGIN_DO_NOT_INFER_TISSUE_OR_OMICS_FROM_ASSAY_KEY")

if __name__=="__main__":
    unittest.main()
