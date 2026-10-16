#!/usr/bin/env python3
"""Offline positive and adversarial tests for exchange v1; no third-party packages."""
import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_multiomics_exchange import validate


def fixture(root):
    def write(path, rows):
        dest = root / path
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows(rows)
        return {"path": path, "sha256": hashlib.sha256(dest.read_bytes()).hexdigest()}

    files = {
        "samples": write("samples.csv", [
            ["assay_id", "sample_id", "subject_id", "omic", "tissue", "collection_time"],
            ["RNA1", "S1", "P1", "transcriptomics", "blood", "0"],
            ["RNA2", "S2", "P2", "transcriptomics", "blood", "1"]
        ]),
        "features": write("features.csv", [
            ["omic", "feature_id", "id_namespace"],
            ["transcriptomics", "ENSG001", "Ensembl"],
            ["transcriptomics", "ENSG002", "Ensembl"]
        ])
    }
    mat = write("matrices/transcriptomics.csv", [
        ["feature_id", "RNA1", "RNA2"],
        ["ENSG001", "3", "NA"],
        ["ENSG002", "4", "9"]
    ])
    manifest = {
        "schema_version": "1.0.0",
        "analysis_id": "test-fixture",
        "study": {"organism": "human", "biological_unit": "participant", "data_access": "public"},
        "files": files,
        "matrices": [{"omic": "transcriptomics", "orientation": "feature_rows_sample_columns",
                      "value_type": "raw_counts", "missing_value": "NA", "file": mat}],
        "software": {"producer": "test", "commit": "0000000", "r_version": None,
                     "python_version": None, "dependency_lockfiles": []},
        "validation": {"status": "exploratory", "methods": []}
    }
    (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return manifest


def assert_invalid(root, fragment):
    problems = validate(root)
    assert any(fragment in problem for problem in problems), (fragment, problems)


def run():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        obj = fixture(root)
        assert validate(root) == [], validate(root)
        (root / "matrices/transcriptomics.csv").write_text(
            "feature_id,RNA1,RNA2\nENSG001,5,NA\nENSG002,4,9\n", encoding="utf-8")
        assert_invalid(root, "SHA-256 mismatch")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        obj = fixture(root)
        obj["matrices"][0]["file"]["path"] = "../outside.csv"
        (root / "manifest.json").write_text(json.dumps(obj), encoding="utf-8")
        assert_invalid(root, "Unsafe relative path")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        obj = fixture(root)
        obj["matrices"][0]["orientation"] = "sample_rows_feature_columns"
        (root / "manifest.json").write_text(json.dumps(obj), encoding="utf-8")
        assert_invalid(root, "Invalid matrix orientation")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        obj = fixture(root)
        p = root / "matrices/transcriptomics.csv"
        p.write_text("feature_id,RNA1,RNA3\nENSG001,3,NA\nENSG002,4,9\n")
        obj["matrices"][0]["file"]["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
        (root / "manifest.json").write_text(json.dumps(obj), encoding="utf-8")
        assert_invalid(root, "Matrix/sample assay_id mismatch")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        obj = fixture(root)
        p = root / "matrices/transcriptomics.csv"
        p.write_text("feature_id,RNA1,RNA2\nENSG001,NaN,NA\nENSG002,4,9\n")
        obj["matrices"][0]["file"]["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
        (root / "manifest.json").write_text(json.dumps(obj), encoding="utf-8")
        assert_invalid(root, "Non-finite numeric value")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        obj = fixture(root)
        p = root / "samples.csv"
        p.write_text("assay_id,sample_id,subject_id,omic,tissue,collection_time\nRNA1,S1,P1,transcriptomics,blood,0\nRNA2,S1,P2,transcriptomics,blood,1\n")
        obj["files"]["samples"]["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
        (root / "manifest.json").write_text(json.dumps(obj), encoding="utf-8")
        assert_invalid(root, "One sample_id assigned to multiple biological subjects")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        obj = fixture(root)
        p = root / "samples.csv"
        p.write_text("assay_id,sample_id,subject_id,omic,tissue,collection_time\nRNA1,S1,P1,transcriptomics,blood,0\nRNA2,S1,P1,transcriptomics,liver,0\n")
        obj["files"]["samples"]["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
        (root / "manifest.json").write_text(json.dumps(obj), encoding="utf-8")
        assert_invalid(root, "Conflicting tissue for shared specimen")
    print("PASS: 1 valid package, 6 adversarial cases")


if __name__ == "__main__":
    run()
