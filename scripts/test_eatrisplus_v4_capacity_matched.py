"""Synthetic fixed-capacity nested CV original physical-source controls."""
import tempfile
import unittest
from pathlib import Path
from scripts.test_eatrisplus_v4_cross_fraction_benchmark import make_source
from scripts.benchmark_eatrisplus_v4_capacity_matched import (
    benchmark_h5mu_matched,EATRISBenchmarkAbstention,
)

class MatchedCapacityScientificTests(unittest.TestCase):
    def test_original_subject_outer_and_nested_ridge_capacity_three(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"mock_v4.h5mu"
            make_source(path,n=42)
            x=benchmark_h5mu_matched(path,strict_original=False,repeats=1)
            self.assertEqual(x["n_exact_matched_original_blood_and_plasma_donors"],40)
            self.assertTrue(x["all_omic_arms_have_exactly_three_total_train_only_PCs"])
            self.assertEqual(x["inner_biological_subject_split"],
                             "3 fold solely within each outer training sample")
            self.assertEqual(len(x["heldout_plasma_protein_target_standardized_MAE"]),6)
            for v in x["heldout_plasma_protein_target_standardized_MAE"].values():
                self.assertTrue(0 < v < 10)
            for a in x["inner_selected_alpha_counts_by_arm"].values():
                self.assertEqual(sum(a.values()),5)
                self.assertTrue(set(a).issubset({"1","10","100","1000","10000","100000","1000000"}))
            self.assertTrue(x["no_omic_winner_selected_on_outer_heldout_errors"])
            self.assertFalse(x["original_drug_exposure_or_germline_PGx_present"])
            self.assertNotIn("SYNTHETIC_PERSON",str(x))
            self.assertEqual(x,benchmark_h5mu_matched(path,strict_original=False,repeats=1))

    def test_original_subject_join_cannot_be_order_dependent(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"mock_v4.h5mu"
            make_source(path,n=42,shuffle=False)
            a=benchmark_h5mu_matched(path,strict_original=False,repeats=1)
            make_source(path,n=42,shuffle=True)
            b=benchmark_h5mu_matched(path,strict_original=False,repeats=1)
            self.assertEqual(a,b)

    def test_invalid_checksum_or_unprespecified_CV_fails(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"mock_v4.h5mu"
            make_source(path)
            with self.assertRaisesRegex(Exception,"md5"):
                benchmark_h5mu_matched(path,strict_original=True,repeats=1)
            with self.assertRaisesRegex(EATRISBenchmarkAbstention,"predeclared"):
                benchmark_h5mu_matched(path,strict_original=False,repeats=5)

if __name__=="__main__":
    unittest.main()
