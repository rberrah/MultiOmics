#!/usr/bin/env python3
"""Documentation-only scientific/graphical README release checks (stdlib)."""
from pathlib import Path
import re
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
IMAGE = ROOT / "docs/assets/multiomics-graphical-abstract.svg"


def main():
    root = ET.parse(IMAGE).getroot()
    ns = "{http://www.w3.org/2000/svg}"
    assert root.tag == ns + "svg"
    assert root.attrib.get("viewBox") == "0 0 1500 800"
    assert root.attrib.get("role") == "img"
    assert root.attrib.get("aria-labelledby") == "ga-title ga-desc"
    assert root.find(ns + "title") is not None
    assert root.find(ns + "desc") is not None
    assert "does not claim clinical PK validation" in (root.find(ns + "desc").text or "")
    assert len(root.findall(".//" + ns + "text")) > 12
    assert IMAGE.stat().st_size < 100_000
    for element in root.iter():
        assert element.tag.rsplit("}", 1)[-1].lower() not in (
            "script", "foreignobject", "image")
        assert all(not (name.endswith("href") and str(value).startswith(
            ("https://", "http://"))) for name, value in element.attrib.items())
    readme = README.read_text(encoding="utf-8")
    assert "](docs/assets/multiomics-graphical-abstract.svg)" in readme
    assert "## MultiOmics → Omics2PK: qualified boundary" in readme
    assert "not" in readme.lower() and "clinical" in readme.lower()
    assert "individual_signatures.csv" in readme
    for link in re.findall(r"\]\((docs/[^)]+)\)", readme):
        assert (ROOT / link.split("#")[0]).is_file(), link
    print("PASS: MultiOmics accessible SVG, research evidence boundary, all README local docs links")


if __name__ == "__main__":
    main()
