#!/usr/bin/env python3
"""True two-reference Entrez Gene ID concordance on 45,101 real GSE22871 probes.

Read the exact SHA-pinned NCBI GPL1261.annot.gz, and independently generated
Bioconductor mouse4302.db/AnnotationDbi Entrez IDs for each original probe.
Never use genotype/PFOS labels to select mappings. No external scores.
"""
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

from qualify_gpl1261_annotation import PINNED_GPL1261_SHA256
from qualify_geo_gse22871 import PINNED_GEO_SOURCE_SHA256, require

MISSING = {"", "NA", "N/A", "---", "NULL", "NONE", "NAN", "?"}
DELIMITERS = re.compile(r"\s*///\s*|\s*;\s*|\s*\|\s*|,\s*")
DIGITS = re.compile(r"[0-9]+")
EXPECTED_PROBES = 45101


def tokens_of_entrez(value):
    """Return IDs, invalid original non-numeric tokens (never silently dropped)."""
    value = (value or "").strip()
    if value.upper() in MISSING:
        return (), ()
    clean, invalid = set(), set()
    for part in DELIMITERS.split(value):
        part = part.strip()
        if part.upper() in MISSING:
            # If mixed with bona fide IDs, retain that it is not fully mapped.
            if len(DELIMITERS.split(value)) > 1:
                invalid.add(part if part else "<empty>")
            continue
        if DIGITS.fullmatch(part):
            clean.add(str(int(part)))
        else:
            invalid.add(part)
    return tuple(sorted(clean, key=lambda x: (len(x), x))), tuple(sorted(invalid))


def classify(geo, bioc, geo_invalid=(), bioc_invalid=()):
    if geo_invalid or bioc_invalid:
        return "unparseable_identifier"
    a, b = set(geo), set(bioc)
    if not a and not b:
        return "both_unmapped"
    if not a:
        return "bioc_only"
    if not b:
        return "geo_only"
    if a == b:
        return "exact_entrez_set"
    if a & b:
        return "partial_entrez_overlap"
    return "disjoint_entrez_sets"


def headkey(name):
    return re.sub(r"[^a-z0-9]", "", name.lower().lstrip("\ufeff"))


def original_geo_entrez(archive):
    raw = Path(archive).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == PINNED_GPL1261_SHA256,
            "Historical original NCBI GPL1261 annotation checksum changed")
    results = {}
    symbols = {}
    header = None
    id_col = gene_col = sym_col = -1
    in_table = False
    finished = False
    with gzip.open(archive, "rt", encoding="utf-8", newline="") as stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if line.lower() in ("!platform_table_begin", "!annotation_table_begin"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.lower() in ("!platform_table_end", "!annotation_table_end"):
                finished = True
                break
            if not line or line.startswith("!"):
                continue
            row = next(csv.reader([line], delimiter="\t", quotechar='"'))
            if header is None:
                header = [headkey(k) for k in row]
                options = [v for v in ("entrezgeneid", "entrezgeneids", "geneid")
                           if v in header]
                require("id" in header and bool(options),
                        "Official GPL1261 lacks a named stable NCBI Entrez Gene ID column: "
                        + repr(header))
                id_col, gene_col = header.index("id"), header.index(options[0])
                sym_col = header.index("genesymbol") if "genesymbol" in header else -1
                continue
            require(len(row) == len(header), "Truncated official platform annotation row")
            probe = row[id_col].strip()
            require(probe and probe not in results,
                    "Invalid or duplicated genuine platform probe ID")
            results[probe] = row[gene_col]
            symbols[probe] = row[sym_col] if sym_col >= 0 else ""
    require(finished and len(results) == EXPECTED_PROBES,
            "GEO official Entrez column does not cover 45,101 original GPL1261 probes")
    return results, symbols


def bioconductor_reference(path):
    results = {}
    with open(path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        require(reader.fieldnames is not None and
                {"probe_id", "geo_symbols", "bioc_symbols", "bioc_entrez_ids",
                 "concordance_status"} <= set(reader.fieldnames),
                "Bioconductor real 45,101-probe reference table missing columns")
        for row in reader:
            probe = row["probe_id"]
            require(probe and probe not in results,
                    "Bioconductor assay key duplicated or missing")
            results[probe] = row
    require(len(results) == EXPECTED_PROBES,
            "Bioconductor reference missing original 45,101 assay probes")
    return results


def validate_evidence(source_dir, bioc_dir):
    source_dir = Path(source_dir)
    bioc_dir = Path(bioc_dir)
    source = json.loads((source_dir / "GPL1261_annotation_mapping_audit.json").read_text())
    cohort = json.loads((source_dir / "external_geo_intake_audit.json").read_text())
    bioc = json.loads((bioc_dir / "GPL1261_mouse4302_db_concordance_audit.json").read_text())
    require(source["official_annotation_gzip_sha256"] == PINNED_GPL1261_SHA256,
            "NCBI source audit is not the independently verified frozen GPL1261")
    require(source["external_cohort_source_sha256"] == PINNED_GEO_SOURCE_SHA256,
            "GEO source study does not match the immutable original GSE22871")
    require(cohort["n_biological_subjects"] == 30,
            "Unexpected independent-cohort animal count")
    require(bioc["probes_compared"] == EXPECTED_PROBES and
            bioc["score_rows_generated"] == 0 and
            bioc["source_geo_annotation_sha256"] == PINNED_GPL1261_SHA256 and
            bioc["source_external_cohort_sha256"] == PINNED_GEO_SOURCE_SHA256,
            "Independent Bioconductor audit is not linked to the exact real sources")
    require(bioc["actual_package_versions"]["mouse4302.db"] == "3.13.0",
            "Bioconductor mouse4302.db annotation version is not locked")


def run(source_dir, bioc_dir):
    source_dir, bioc_dir = Path(source_dir), Path(bioc_dir)
    validate_evidence(source_dir, bioc_dir)
    geo_ids, geo_symbols = original_geo_entrez(source_dir / "GPL1261.annot.gz")
    bio = bioconductor_reference(
        bioc_dir / "GPL1261_GEO_vs_mouse4302_db_all_probe_concordance.csv"
    )
    require(set(geo_ids) == set(bio),
            "NCBI and Bioconductor original 45,101 platform probe identities differ")
    counters = Counter()
    salvage = Counter()
    ambiguous = 0
    rows = []
    for probe in sorted(bio):
        geo, invalid_geo = tokens_of_entrez(geo_ids[probe])
        b_ids, invalid_b = tokens_of_entrez(bio[probe]["bioc_entrez_ids"])
        status = classify(geo, b_ids, invalid_geo, invalid_b)
        counters[status] += 1
        symbol_status = bio[probe]["concordance_status"]
        if status == "exact_entrez_set" and symbol_status != "exact_symbol_set_agreement":
            salvage[symbol_status] += 1
        if len(geo) > 1 or len(b_ids) > 1:
            ambiguous += 1
        rows.append({
            "probe_id": probe,
            "geo_entrez_raw": geo_ids[probe],
            "geo_entrez_ids": "|".join(geo),
            "bioc_entrez_ids": "|".join(b_ids),
            "geo_invalid_tokens": "|".join(invalid_geo),
            "bioc_invalid_tokens": "|".join(invalid_b),
            "geo_symbol_raw": geo_symbols[probe],
            "geo_symbols_existing_audit": bio[probe]["geo_symbols"],
            "bioc_symbols_existing_audit": bio[probe]["bioc_symbols"],
            "symbol_concordance": symbol_status,
            "entrez_concordance": status,
            "same_stable_id_despite_different_symbol": (
                status == "exact_entrez_set"
                and symbol_status != "exact_symbol_set_agreement"
            ),
            "multi_entrez_ambiguous": len(geo) > 1 or len(b_ids) > 1,
        })
    require(len(rows) == EXPECTED_PROBES and sum(counters.values()) == EXPECTED_PROBES,
            "Comparisons did not cover every genuine external assay probe")
    require(counters["exact_entrez_set"] > 0,
            "No genuine two-reference stable-ID matches, likely invalid parsing")
    require(counters["unparseable_identifier"] < 5000,
            "Too many nonnumeric GEO Entrez IDs: parsing may be incorrect")
    out_dir = bioc_dir / "original_NCBI_vs_Bioc_Entrez_ID"
    out_dir.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0])
    with (out_dir / "GPL1261_all_probe_entrez_concordance.csv").open(
            "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, columns)
        writer.writeheader()
        writer.writerows(rows)
    discord = [r for r in rows
               if r["entrez_concordance"] != "exact_entrez_set"]
    with (out_dir / "GPL1261_entrez_disagreements_all.csv").open(
            "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, columns)
        writer.writeheader()
        writer.writerows(discord)
    anchor = [r for r in rows if "cyp3a11" in geo_symbols[r["probe_id"]].lower()
              or re.search(r"(^|\|)cyp3a11(\||$)",r["bioc_symbols_existing_audit"],
                           re.IGNORECASE)]
    with (out_dir / "Cyp3a11_entrez_candidate_reference.csv").open(
            "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, columns)
        writer.writeheader()
        writer.writerows(anchor)
    summary = {
        "source_series": "GSE22871",
        "platform": "GPL1261",
        "source_probe_measurements": EXPECTED_PROBES,
        "source_original_gpl1261_sha256": PINNED_GPL1261_SHA256,
        "source_original_geo_series_sha256": PINNED_GEO_SOURCE_SHA256,
        "references": [
            "Original NCBI GEO GPL1261.annot.gz ENTREZ_GENE_ID (historical)",
            "Independent Bioconductor mouse4302.db 3.13.0 ENTREZID"
        ],
        "entrez_concordance_counts": dict(sorted(counters.items())),
        "unparseable_geo_or_bioc_id_count": counters["unparseable_identifier"],
        "ambiguous_multiple_ncbi_ids_probe_count": ambiguous,
        "same_entrez_ID_set_but_gene_symbol_labels_differ": sum(salvage.values()),
        "same_entrez_ID_by_symbol_difference_type": dict(sorted(salvage.items())),
        "cyp3a11_candidate_probes": len(anchor),
        "observed_biological_subjects": 30,
        "all_id_compared_without_outcome_labels": True,
        "gene_annotation_independent_of_phenotype": True,
        "no_gene_scale_transport_fitted": True,
        "no_technical_probe_collapse_selected": True,
        "external_scores_produced": 0,
        "external_pk_outcome_available": False,
        "status": "STABLE_ENTREZ_ID_COMPARISON_ONLY_EXTERNAL_PROJECTION_BLOCKED",
        "limitations": [
            "NCBI Entrez IDs are not invariant through time; retired/replaced IDs can disagree.",
            "Same Entrez ID does not establish probe-specificity or compatible expression scales.",
            "Historical official GPL1261 gene IDs and current Bioconductor annotation are distinct dated sources.",
            "Gene symbols and Entrez IDs can be ambiguous for probes matching multiple loci.",
            "No external phenotype labels used for selection, no missing gene imputation, no PK claim."
        ]
    }
    (out_dir / "GPL1261_entrez_id_concordance_audit.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("REAL NCBI GPL1261 vs Bioconductor Entrez concordance PASS: "
          f"all {EXPECTED_PROBES} original probes; exact={counters['exact_entrez_set']} "
          f"partial={counters['partial_entrez_overlap']} "
          f"disjoint={counters['disjoint_entrez_sets']} "
          f"ncbiOnly={counters['geo_only']} "
          f"biocOnly={counters['bioc_only']} "
          f"bothUnmapped={counters['both_unmapped']} "
          f"unparseable={counters['unparseable_identifier']} "
          f"stableID_agrees_symbol_differs={sum(salvage.values())}; "
          "external scores=0 (correct).")


def selftest():
    assert tokens_of_entrez("123 /// 456 /// 123") == (("123", "456"), ())
    assert tokens_of_entrez("---") == ((), ())
    assert tokens_of_entrez("123 /// foo") == (("123",), ("foo",))
    assert tokens_of_entrez("789;123") == (("123", "789"), ())
    assert classify(("23",), ("23",)) == "exact_entrez_set"
    assert classify(("23", "42"), ("23", "50")) == "partial_entrez_overlap"
    assert classify(("23",), ("42",)) == "disjoint_entrez_sets"
    assert classify((), ("5",)) == "bioc_only"
    assert classify(("5",), ()) == "geo_only"
    assert classify((), ()) == "both_unmapped"
    assert classify(("23",), ("23",), ("garbled",), ()) == "unparseable_identifier"
    print("PASS: fail-closed numeric Entrez parsing, all seven concordance groups, multi-target and invalid-ID retention")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        selftest()
    elif len(sys.argv) == 3:
        run(sys.argv[1], sys.argv[2])
    else:
        sys.exit("Usage: qualify_gpl1261_entrez_concordance.py --self-test | GEO_SOURCE_DIR BIOC_AUDIT_DIR")
