"""Synthetic-only controls for followup original n90 cross-compartment study."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scripts.benchmark_eatrisplus_v4_orthogonal_capacity_shrinkage import (
    _factorial_design,benchmark_original,ALPHAS,CAPACITIES,ARMS,
)
from scripts.test_eatrisplus_v4_cross_fraction_benchmark import make_source

class OrthogonalEATRISSyntheticContracts(unittest.TestCase):
    def test_exact_3_and_9_components_and_no_test_subject_PCA_refit(self):
        rng=np.random.default_rng(79)
        rna=rng.normal(size=(40,25))
        met=rng.normal(size=(40,20))
        train=np.arange(30);test=np.arange(30,40)
        first=_factorial_design(rna,met,train,test)
        self.assertEqual(set(first),{3,9})
        for cap in CAPACITIES:
            for arm in ARMS:
                x,y=first[cap][arm]
                self.assertEqual(x.shape,(30,cap))
                self.assertEqual(y.shape,(10,cap))
        # Perturb ONLY heldout molecular values: train PCA output
        # must not change at any setting, model may predict differently.
        changed=rna.copy();changed[test]+=rng.normal(size=changed[test].shape)*500
        other=_factorial_design(changed,met,train,test)
        for cap in CAPACITIES:
            for arm in ARMS:
                np.testing.assert_allclose(first[cap][arm][0],other[cap][arm][0],
                                           atol=1e-8,rtol=1e-8)

    def test_original_donor_shuffled_source_index_has_same_result(self):
        with tempfile.TemporaryDirectory() as work:
            path=Path(work)/"synthetic.h5mu"
            make_source(path,n=42)
            a=benchmark_original(path,verify_source=False,repeats=1)
            self.assertEqual(a["original_matched_unique_persons"],40)
            self.assertEqual(len(a["heldout_protein_profile_MAE_by_factorial_arm"]),13)
            self.assertTrue(a["model_setting_chosen_after_previous_original_experiments"])
            self.assertFalse(a["source_original_MD5_verified"])
            self.assertNotIn("SYNTHETIC_PERSON",str(a))
            self.assertEqual(a,benchmark_original(path,verify_source=False,repeats=1))
            make_source(path,n=42,shuffle=True)
            b=benchmark_original(path,verify_source=False,repeats=1)
            self.assertEqual(a,b)

    def test_non_original_data_cannot_claim_original_replay(self):
        with tempfile.TemporaryDirectory() as work:
            path=Path(work)/"synthetic.h5mu"
            make_source(path,n=42)
            with self.assertRaisesRegex(Exception,"md5"):
                benchmark_original(path,verify_source=True,repeats=1)
            with self.assertRaises(Exception):
                benchmark_original(path,verify_source=False,repeats=4)

if __name__=="__main__":
    unittest.main()
