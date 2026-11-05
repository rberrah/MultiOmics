"""Synthetic anti-leakage controls for two original EATRIS protein panel sensitivity."""
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from scripts.benchmark_eatrisplus_v4_protein_panels import (
    benchmark_proteome_panels, EATRISBenchmarkAbstention, TARGET_PANELS,
)
from scripts.test_eatrisplus_v4_cross_fraction_benchmark import make_source


def make_two_panel_source(path, *, reorder=False, altered_ipp=False):
    make_source(path, n=42)
    rng = np.random.default_rng(20261105)
    with h5py.File(path, "r+") as h:
        y = h["mod"]["proteins"]["X"][()]
        source_obs = [
            k.decode("utf-8") if isinstance(k, bytes) else str(k)
            for k in h["mod"]["proteins"]["obs"]["_index"][()]
        ]
        # These are different *measured targets*, not a duplicate of proteins.
        ipp = np.hstack([
            y[:, :15] * .30 + rng.normal(size=(42, 15)) * .9,
            rng.normal(size=(42, 30)),
        ])
        if altered_ipp:
            ipp = rng.normal(size=ipp.shape)
        rows = [i for i in range(42) if i not in (5, 41)]
        if reorder:
            rows = rows[::-1]
        m = h["mod"].create_group("proteomics-ipp")
        m.create_group("obs").create_dataset("_index", data=[
            source_obs[i].encode("utf-8") for i in rows
        ])
        m.create_group("var").create_dataset("_index", data=[
            f"SYNTHETIC_IPP_{i}".encode() for i in range(ipp.shape[1])
        ])
        m.create_dataset("X", data=ipp[rows])


class MatchedProteinPanelSyntheticContracts(unittest.TestCase):
    def test_same_person_panels_and_reproducible_splits(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "synthetic.h5mu"
            make_two_panel_source(p)
            a = benchmark_proteome_panels(p, verify_source=False, repeats=1)
            self.assertEqual(
                a["n_original_independent_matched_people_in_all_five_assays"], 39)
            self.assertEqual(set(a["panel_results"]), set(TARGET_PANELS))
            self.assertEqual(a["fixed_total_PC_budgets"], [3, 9])
            self.assertTrue(a["identical_donor_split_and_input_representation_all_panels"])
            self.assertFalse(a["original_germline_full_PGx_or_human_drug_PK_measured"])
            self.assertFalse(a["original_h5mu_MD5_verified"])
            self.assertNotIn("SYNTHETIC_PERSON", str(a))
            self.assertNotIn("SYNTHETIC_IPP_", str(a))
            self.assertEqual(a, benchmark_proteome_panels(
                p, verify_source=False, repeats=1))
            make_two_panel_source(p, reorder=True)
            self.assertEqual(a, benchmark_proteome_panels(
                p, verify_source=False, repeats=1))

    def test_target_panel_replacement_does_not_change_other_panel(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "synthetic.h5mu"
            make_two_panel_source(p)
            original = benchmark_proteome_panels(p, verify_source=False, repeats=1)
            make_two_panel_source(p, altered_ipp=True)
            changed = benchmark_proteome_panels(p, verify_source=False, repeats=1)
            self.assertEqual(
                original["panel_results"]["proteins"],
                changed["panel_results"]["proteins"])
            self.assertNotEqual(
                original["panel_results"]["proteomics-ipp"],
                changed["panel_results"]["proteomics-ipp"])

    def test_unverified_source_cannot_be_presented_as_original(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "synthetic.h5mu"
            make_two_panel_source(p)
            with self.assertRaisesRegex(Exception, "md5"):
                benchmark_proteome_panels(p, verify_source=True, repeats=1)
            with self.assertRaisesRegex(EATRISBenchmarkAbstention, "1 to 3"):
                benchmark_proteome_panels(p, verify_source=False, repeats=9)
            with h5py.File(p, "r+") as h:
                del h["mod"]["proteomics-ipp"]
            with self.assertRaisesRegex(EATRISBenchmarkAbstention, "not present"):
                benchmark_proteome_panels(p, verify_source=False, repeats=1)


if __name__ == "__main__":
    unittest.main()
