"""Synthetic mock public-MuData quantitative cross-fraction prediction tests."""
import tempfile
import unittest
from pathlib import Path

from scripts.benchmark_eatrisplus_v4_cross_fraction import (
    benchmark_h5mu,EATRISBenchmarkAbstention,
)


def make_source(path,*,n=42,shuffle=False,unmatched=False):
    import h5py
    import numpy as np
    g=np.random.default_rng(204)
    z=g.normal(size=(n,4))
    key=[f"SYNTHETIC_PERSON_{i:03d}" for i in range(n)]
    RNA=np.concatenate((z@np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,1]])+
                   .05*g.normal(size=(n,3)),g.normal(size=(n,17))),axis=1)
    A=np.concatenate((z[:,:2]+.1*g.normal(size=(n,2)),g.normal(size=(n,8))),axis=1)
    M=np.concatenate((z[:,2:]+.1*g.normal(size=(n,2)),g.normal(size=(n,8))),axis=1)
    Y=np.concatenate((z@np.array([[1]*15,[.5]*15,[-.5]*15,[.25]*15])+
                    .1*g.normal(size=(n,15)),g.normal(size=(n,10))),axis=1)
    with h5py.File(path,"w") as h:
        h.create_group("obs").create_dataset("_index",data=[s.encode() for s in key])
        h.create_group("var").create_dataset("_index",data=[b"x"])
        mods=h.create_group("mod")
        for name,X in (
            ("mRNA",RNA),("acylcarnitines",A),
            ("aminoacids",M),("proteins",Y),
        ):
            sample=list(range(n))
            if name=="mRNA":sample=sample[:-2]
            if unmatched and name=="proteins":sample=sample[2:]
            if shuffle and name=="proteins":sample=sample[::-1]
            md=mods.create_group(name)
            md.create_group("obs").create_dataset(
                "_index",data=[key[i].encode() for i in sample])
            md.create_group("var").create_dataset(
                "_index",data=[f"{name}_{i}".encode() for i in range(X.shape[1])])
            md.create_dataset("X",data=X[sample])
    return key


class SyntheticPublicEATRISCrossCompartment(unittest.TestCase):
    def test_original_exact_donor_crossview_and_fixed_omics_ablation(self):
        with tempfile.TemporaryDirectory() as work:
            p=Path(work)/"fake.h5mu"
            make_source(p)
            result=benchmark_h5mu(p,strict_original=False,repeats=2)
            self.assertEqual(result["n_exact_same_original_subjects_in_all_four_model_input_output_assays"],40)
            self.assertEqual(result["source_MuData_modality_count"],4)
            self.assertEqual(len(result["heldout_standardized_plasma_protein_profile_MAE_by_arm"]),5)
            vals=result["heldout_standardized_plasma_protein_profile_MAE_by_arm"]
            # An intentionally finite/small synthetic cohort does NOT
            # guarantee additional omics improves generalization.
            # Regression must return finite comparable errors and
            # explicitly allow a negative/superfluous added-omic result.
            self.assertTrue(all(0.0 <= score < 20.0 for score in vals.values()))
            self.assertEqual(result,benchmark_h5mu(p,strict_original=False,repeats=2))
            self.assertFalse(result["original_germline_CYP3A_PGx_publicly_confirmed"])
            self.assertTrue(result["measured_endpoint_is_not_clinical_drug_PK_or_PD"])
            self.assertNotIn("SYNTHETIC_PERSON",str(result))
            ci=result["exploratory_original_biological_person_paired_contrasts"]
            self.assertTrue(all(not x["vs_train_mean_only"]
                                ["bootstrap_refits_models_or_adjusts_model_choice"]
                                for x in ci.values()))

    def test_original_index_reordering_only_never_changes_donor_prediction(self):
        with tempfile.TemporaryDirectory() as work:
            p=Path(work)/"fake.h5mu"
            make_source(p,shuffle=False)
            a=benchmark_h5mu(p,strict_original=False,repeats=1)
            make_source(p,shuffle=True)
            b=benchmark_h5mu(p,strict_original=False,repeats=1)
            self.assertEqual(a,b)

    def test_incompatible_source_no_claim_and_no_full_original_PK(self):
        with tempfile.TemporaryDirectory() as work:
            p=Path(work)/"fake.h5mu"
            make_source(p)
            with self.assertRaisesRegex(Exception,"md5"):
                benchmark_h5mu(p,strict_original=True,repeats=1)
            with self.assertRaisesRegex(EATRISBenchmarkAbstention,
                                        "Only fixed"):
                benchmark_h5mu(p,strict_original=False,repeats=8)


if __name__=="__main__":
    unittest.main()
