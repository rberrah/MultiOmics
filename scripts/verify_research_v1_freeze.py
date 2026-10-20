#!/usr/bin/env python3
"""Fail-closed Research v1 release consistency checks."""
from __future__ import annotations
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
RELEASE_VERSION = "1.0.0"
POST_V1_PATTERN = r"^1\.1\.[0-2](?:-dev\.\d+)?$"
DATE = "2026-10-19"

def load_json(path: str):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))

def main():
    pkg = load_json("package.json")
    lock = load_json("package-lock.json")
    manifest = load_json("research-v1-manifest.json")
    runtime = load_json("multiomics-engine/environment/reference-runtime.json")
    inventory = load_json("multiomics-engine/environment/reference-runtime-inventory.json")
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    freeze = (ROOT / "docs/multiomics/RESEARCH_V1_FREEZE_20261019.md").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")

    current_version = pkg["version"]
    assert lock["version"] == current_version
    assert lock["packages"][""]["version"] == current_version
    assert current_version == RELEASE_VERSION or re.match(POST_V1_PATTERN, current_version)
    assert manifest["release"]["softwareVersion"] == RELEASE_VERSION
    assert manifest["release"]["freezeDate"] == DATE
    assert manifest["release"]["intendedTag"] == "research-v1.0.0"
    assert re.search(rf"(?m)^version:\s*[\"']?{re.escape(current_version)}[\"']?\s*$", citation)
    expected_citation_date = DATE if current_version == RELEASE_VERSION else load_json("release-v1.1-manifest.json")["release"]["releaseDate"]
    assert re.search(rf"(?m)^date-released:\s*[\"']?{re.escape(expected_citation_date)}[\"']?\s*$", citation)

    assert runtime["rVersion"] == "4.5.2"
    assert runtime["bioconductorVersion"] == "3.22"
    assert runtime["strictLockReady"] is True
    assert "sha256:e34177896c1579f66d48b756e4229b58237499183a089da084af72d7aba1d40a" in runtime["container"]
    assert runtime["fullInstalledLibrary"]["packageCount"] == 171
    assert runtime["fullInstalledLibrary"]["sha256SortedPackageVersionLines"] == "02798608a29db003dba1a5f0d5116eb52f8c2a21e8d19401be2102db079c84db"
    assert inventory["schema"] == "multiomics.r-full-package-inventory.v1"
    assert inventory["packageCount"] == 171
    assert len(inventory["packages"]) == 171
    assert inventory["rVersion"] == "4.5.2"
    assert inventory["bioconductorVersion"] == "3.22"
    assert manifest["frozenRuntime"]["fullPackageInventory"] == "multiomics-engine/environment/reference-runtime-inventory.json"
    assert manifest["frozenRuntime"]["fullPackageCount"] == 171
    assert manifest["frozenRuntime"]["fullPackageInventorySha256"] == runtime["fullInstalledLibrary"]["sha256SortedPackageVersionLines"]
    assert manifest["frozenRuntime"]["cranSnapshot"] == "2026-10-16"

    required_files = [
        ".github/workflows/standalone-ci.yml",
        ".github/workflows/cox-reference.yml",
        ".github/workflows/longitudinal-reference.yml",
        ".github/workflows/public-airway-replay.yml",
        ".github/workflows/nutrimouse-mofa-validation.yml",
        ".github/workflows/nutrimouse-diablo-validation.yml",
        ".github/workflows/sacurine-zenodo-12raw-xcms.yml",
        ".github/workflows/external-time-survival-research-v1.yml",
        "docs/multiomics/EXTERNAL_TIME_SURVIVAL_RESEARCH_V1.md",
        "docs/multiomics/RESEARCH_V1_CONFIRMATORY_RELEASE_GATE.md",
        "docs/multiomics/ARTICLE_VALIDATION_EVIDENCE_MATRIX_20261016.md",
        "multiomics-engine/environment/reference-runtime-inventory.json",
        "static/multiomics/scientific-validation.md",
    ]
    missing = [p for p in required_files if not (ROOT / p).is_file()]
    assert not missing, f"Missing Research v1 evidence files: {missing}"

    issues = {row["issue"] for row in manifest["knownOpenScientificLimitations"]}
    assert issues == {20, 44, 47}
    assert manifest["scope"]["medicalDevice"] is False
    assert manifest["scope"]["clinicalDecisionSupportCertified"] is False
    assert manifest["scope"]["universalConfirmatoryValidity"] is False

    boundary = manifest["releaseBoundary"].lower()
    for required in ("does not certify", "statistical calibration", "biological", "clinical"):
        assert required in boundary

    combined = (freeze + "\n" + readme).lower()
    assert "research v1" in combined
    assert "medical-device" in freeze.lower() or "medical device" in freeze.lower()
    assert "fdr" in freeze.lower()
    assert "issue #44" in freeze.lower()
    assert "issue #47" in freeze.lower()
    assert "issue #20" in freeze.lower()
    assert "human" in freeze.lower()
    assert "pmx" in freeze.lower() and "separate" in freeze.lower()

    # The release note must explicitly refuse common overclaims.
    forbidden_release_claims = [
        "universally validated fdr",
        "clinically validated software",
        "publication ready automatically",
        "medical device certified",
    ]
    assert not any(x in freeze.lower() for x in forbidden_release_claims)

    if current_version != RELEASE_VERSION:
        assert "release/research-v1.0.0" in readme
        assert "post-v1" in readme.lower()
        assert current_version in readme
    print(
        "PASS: frozen Research v1 evidence remains internally consistent; "
        f"current source version is {current_version}."
    )

if __name__ == "__main__":
    main()
