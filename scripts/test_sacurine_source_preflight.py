#!/usr/bin/env python3
"""No-network adversarial unit tests for the Sacurine source verification."""
from io import BytesIO
from pathlib import Path
import tempfile
import zipfile

from sacurine_source_preflight import identify_role, inspect_archive, inspect_metadata

assert identify_role("QC1_001.mzXML") == "pooled_qc"
assert identify_role("Blanc04.mzXML") == "blank"
assert identify_role("HU_neg_048.mzXML") == "biological"

def fixture(path, altered=None):
    names = (
        [f"HU_neg_{i:03d}.mzXML" for i in range(1, 21)] +
        [f"QC1_{i:03d}.mzXML" for i in range(1, 11)] +
        [f"Blanc{i:02d}.mzXML" for i in range(1, 7)]
    )
    if altered is not None:
        names[-1] = altered
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name in names:
            z.writestr(name, b"<mzXML />")
    return names

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    archive = root / "example.zip"
    fixture(archive)
    rows, others = inspect_archive(archive)
    assert len(rows) == 36 and not others
    assert sum(x["role_provisional"] == "blank" for x in rows) == 6
    meta = root / "source.tsv"
    meta.write_text("filename\tclass\n" + "".join(
        f"{r['filename']}\t{r['role_provisional']}\n" for r in rows
    ))
    audit = inspect_metadata(meta, rows)
    assert audit["exact_sample_name_matches"]["filename"] == 36
    assert not audit["official_roles_confirmed"]  # preflight is not conclusive
    for unsafe in ("../not-allowed.mzXML", "/absolute.mzXML",
                   "unknown-file.mzXML", "hu_neg_001.mzXML"):
        fixture(archive, unsafe)
        try:
            inspect_archive(archive)
        except ValueError:
            pass
        else:
            raise AssertionError("Unsafe filename, role or collision accepted: " + unsafe)
    fixture(archive, "HU_neg_999.mzXML")  # count balance fails
    try:
        inspect_archive(archive)
    except ValueError:
        pass
    else:
        raise AssertionError("Absent final blank was accepted")
print("PASS: 36-source ZIP, strict sample roles, traversal rejection, duplicated source and no invented meta certification")
