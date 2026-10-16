#!/usr/bin/env python3
"""No network/no R packages: acceptance and adversarial tests of replay comparator."""
from verify_multiomics_r_replay import verify

original = {
    "methods": {
        "transcriptomics_differential_deseq2": {
            "status": "ok",
            "top": [
                {"feature": "ENSG001", "log2FoldChange": 1.2345, "pvalue": 0.01, "padj": 0.03},
                {"feature": "ENSG002", "log2FoldChange": None, "pvalue": None, "padj": None}
            ],
            "referenceTableRows": 2,
            "referenceTableComplete": True,
        }
    },
    "preprocessing": {"transcriptomics": ["raw_counts", "independent_filtering"]},
    "analysisMode": "single_omic",
    "analysedLayers": ["transcriptomics"]
}
import copy
matching = copy.deepcopy(original)
assert verify(original, matching) == []
matching["methods"]["transcriptomics_differential_deseq2"]["top"][0]["log2FoldChange"] += 0.1
assert any("numerical mismatch" in x for x in verify(original, matching))
matching = copy.deepcopy(original)
matching["methods"]["transcriptomics_differential_deseq2"]["top"].pop()
assert any("row/list count" in x for x in verify(original, matching))
matching = copy.deepcopy(original)
matching["preprocessing"]["transcriptomics"] = ["raw_counts"]
assert any("preprocessing" in x for x in verify(original, matching))
matching = copy.deepcopy(original)
matching["methods"]["transcriptomics_differential_deseq2"]["top"][1]["padj"] = 0
assert verify(original, matching)
matching = copy.deepcopy(original)
matching["methods"]["transcriptomics_differential_deseq2"]["status"] = "error"
assert verify(original, matching)
assert verify({"methods": {}}, {"methods": {}})
print("PASS: R replay numerical comparator detects row, NA, preprocessing and method drift")
