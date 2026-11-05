"""Offline EATRIS source provenance tests: no public subject-level data."""
import io
import tempfile
import unittest
from unittest.mock import patch
from scripts.audit_eatris_plus_v4_original_h5mu import (
    EXPECTED_MD5,EXPECTED_NAME,EXPECTED_RECORD_ID,EATRISSourceError,
    _source_metadata,audit_h5mu,
)

class EATRISSourceContract(unittest.TestCase):
    def test_metadata_pin_source_record_and_md5(self):
        fake={"id":EXPECTED_RECORD_ID,
              "files":[{"key":EXPECTED_NAME,"checksum":"md5:"+EXPECTED_MD5,
                        "size":141_700_000,
                        "links":{"self":"https://zenodo.org/api/records/19136013/files/test/content"}}]}
        class Response:
            def __enter__(self):return self
            def __exit__(self,*args):pass
            def read(self,*args):return b""
        with patch("scripts.audit_eatris_plus_v4_original_h5mu.urlopen",
                   return_value=io.BytesIO(__import__("json").dumps(fake).encode())):
            url,size=_source_metadata()
            self.assertIn("zenodo.org",url)
            self.assertEqual(size,141_700_000)
        fake["files"][0]["checksum"]="md5:"+"0"*32
        with patch("scripts.audit_eatris_plus_v4_original_h5mu.urlopen",
                   return_value=io.BytesIO(__import__("json").dumps(fake).encode())):
            with self.assertRaisesRegex(EATRISSourceError,"hash"):
                _source_metadata()
        fake["files"][0]["checksum"]="md5:"+EXPECTED_MD5
        fake["id"]=19136012
        with patch("scripts.audit_eatris_plus_v4_original_h5mu.urlopen",
                   return_value=io.BytesIO(__import__("json").dumps(fake).encode())):
            with self.assertRaisesRegex(EATRISSourceError,"record ID"):
                _source_metadata()

    def test_hdf5_modality_schema_from_synthetic_structure(self):
        try:import h5py
        except ImportError:self.skipTest("h5py optional source-reader extra")
        with tempfile.NamedTemporaryFile(suffix=".h5mu") as f:
            with h5py.File(f.name,"w") as h:
                x=h.create_group("mod").create_group("RNA")
                x.attrs["encoding-type"]="anndata"
                x.create_dataset("X",data=[[1.,2.],[3.,4.]])
                ob=x.create_group("obs")
                ob.create_dataset("_index",data=[b"SENSITIVE_SUBJECT_1",b"SENSITIVE_SUBJECT_2"])
            s=audit_h5mu(f.name)
        self.assertEqual(s["n_original_modality_groups"],1)
        self.assertTrue(s["original_individual_donor_labels_or_values_emitted"] is False)
        self.assertNotIn("SENSITIVE_SUBJECT_1",str(s))
        self.assertEqual(s["original_modality_group_schema"]["RNA"]["HDF5_group_schema"]
                         ["X"]["shape"],[2,2])

    def test_missing_modality_abstains(self):
        try:import h5py
        except ImportError:self.skipTest("h5py optional source-reader extra")
        with tempfile.NamedTemporaryFile(suffix=".h5mu") as f:
            with h5py.File(f.name,"w") as h:h.create_group("obs")
            with self.assertRaisesRegex(EATRISSourceError,"mod group"):
                audit_h5mu(f.name)
if __name__=="__main__":unittest.main()
