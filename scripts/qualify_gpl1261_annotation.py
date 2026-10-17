#!/usr/bin/env python3
"""Audit NCBI GEO GPL1261 probe-to-gene annotation without scoring phenotypes.

Input GSE22871 original processed probe matrix has already been verified by
qualify_geo_gse22871.py. Keep ambiguous and missing source mappings visible.
A mapped gene symbol is NOT cross-platform normalization/PK validation.
"""
import csv
import gzip
import hashlib
import io
import json
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from qualify_geo_gse22871 import PINNED_GEO_SOURCE_SHA256, require

ANNOT_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
             "GPL1nnn/GPL1261/annot/GPL1261.annot.gz")
ANNOT_LIMIT = 35_000_000
# Frozen byte-for-byte after real NCBI download and complete 45,101-row audit.
PINNED_GPL1261_SHA256 = "394c103cf4ca99500d06369e4fb2b9b9bffc2c10fd3cbeb0b0d64c51c3fc95f8"
INVALID = {"", "---", "NA", "N/A", "NULL", "NONE", "NAN", "?"}
MAX_PROBES = 55000


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def verified_matrix_probes(matrix_file):
    compressed = Path(matrix_file).read_bytes()
    require(sha(compressed) == PINNED_GEO_SOURCE_SHA256,
            "Prior qualified GSE22871 source bytes changed before annotation audit")
    header_seen = False
    done = False
    ids = set()
    with gzip.open(io.BytesIO(compressed), "rt", encoding="utf-8") as stream:
        for raw in stream:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                header_seen = True
                continue
            if not header_seen:
                continue
            if line == "!series_matrix_table_end":
                done = True
                break
            if line.startswith('"ID_REF"') or line.startswith("ID_REF\t"):
                continue
            probe = next(csv.reader([line], delimiter="\t", quotechar='"'))[0]
            require(bool(probe) and probe not in ids, "Repeated/empty original GSE probe")
            ids.add(probe)
    require(done and len(ids) == 45101, "Exact 45,101 qualified probe IDs missing")
    return ids


def download_annotation():
    request = urllib.request.Request(
        ANNOT_URL, headers={"User-Agent": "MultiOmics-source-audit/1.0",
                           "Accept": "application/gzip"}
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        chunks = []
        n = 0
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            n += len(chunk)
            require(n <= ANNOT_LIMIT, "NCBI GEO GPL1261 annotation compressed size exceeded")
            chunks.append(chunk)
    raw = b"".join(chunks)
    require(raw.startswith(b"\x1f\x8b"),
            "NCBI response was not the authentic GPL1261 annotation gzip")
    return raw


def canonical_column(name):
    return re.sub(r"[^a-z0-9]", "", name.lower().lstrip("\ufeff"))


def parse_symbols(cell):
    cell = (cell or "").strip()
    if cell.upper() in INVALID:
        return []
    # Delimiters explicitly used for multiple entries in official probe mappings.
    fields = re.split(r"\s*///\s*|\s*;\s*|\s+\|\s+", cell)
    fields = [v.strip() for v in fields if v.strip()]
    # Do not invent a single gene by selecting the first of ambiguous symbols.
    if any(v.upper() in INVALID or not re.fullmatch(r"[A-Za-z0-9_.-]+", v)
           for v in fields):
        return []
    return sorted(set(fields), key=str.lower)


def parse_gpl_annotation(lines, matrix_probes, *, min_source_rows=30000):
    in_table = False
    ended = False
    header = None
    rows = 0
    matched = set()
    annotated_ids = set()
    symbol_to_probes = defaultdict(set)
    probe_to_symbols = {}
    for raw in lines:
        line = raw.rstrip("\r\n")
        if not in_table:
            if line.lower() in ("!platform_table_begin", "!annotation_table_begin"):
                in_table = True
            continue
        if line.lower() in ("!platform_table_end", "!annotation_table_end"):
            ended = True
            break
        if not line or line.startswith("!"):
            continue
        fields = next(csv.reader([line], delimiter="\t", quotechar='"'))
        if header is None:
            header = [canonical_column(v) for v in fields]
            require("id" in header, "Official GPL annotation has no probe ID column")
            symbol_headers = ("genesymbol", "symbol", "genesymbols")
            options = [key for key in symbol_headers if key in header]
            require(bool(options),
                    "Official GPL annotation must include a named gene-symbol column: "
                    + repr(header))
            symbol_index = header.index(options[0])
            id_index = header.index("id")
            continue
        require(len(fields) == len(header), "Ragged GEO GPL platform annotation row")
        probe = fields[id_index].strip()
        require(probe and probe not in annotated_ids,
                "Duplicate or empty official probe annotation row")
        annotated_ids.add(probe)
        rows += 1
        require(rows <= MAX_PROBES, "Unexpected official GPL probe count")
        if probe not in matrix_probes:
            continue
        matched.add(probe)
        symbols = parse_symbols(fields[symbol_index])
        probe_to_symbols[probe] = symbols
        for sym in symbols:
            symbol_to_probes[sym].add(probe)
    require(header is not None and ended and rows >= min_source_rows,
            "Official GPL1261 table unavailable, truncated or too small")
    require(len(matched) >= min(35000, min_source_rows),
            "Too few original probes have source-backed official annotations")
    empty = sorted(matrix_probes - matched)
    ambiguous = sorted(probe for probe, symbols in probe_to_symbols.items()
                       if len(symbols) > 1)
    symbol_less = sorted(probe for probe, symbols in probe_to_symbols.items()
                         if not symbols)
    multi_probe_genes = sum(1 for probes in symbol_to_probes.values() if len(probes) > 1)
    return {
        "header": header,
        "official_annotation_rows": rows,
        "original_matrix_probe_rows": len(matrix_probes),
        "matrix_probe_ids_with_official_row": len(matched),
        "matrix_probes_without_annotation": empty,
        "matrix_probes_with_no_unique_parseable_gene_symbols": symbol_less,
        "multi_gene_ambiguous_probe_ids": ambiguous,
        "distinct_annotated_gene_symbols": len(symbol_to_probes),
        "gene_symbols_with_multiple_measurement_probes": multi_probe_genes,
        "probe_to_symbols": probe_to_symbols,
        "gene_to_probes": symbol_to_probes,
    }


def write_csv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as stream:
        wr = csv.writer(stream)
        wr.writerow(columns)
        wr.writerows(rows)


def qualify(outdir):
    out = Path(outdir)
    require((out / "GSE22871_series_matrix.txt.gz").exists(),
            "Original externally qualified GSE22871 source artifact is absent")
    probes = verified_matrix_probes(out / "GSE22871_series_matrix.txt.gz")
    content = download_annotation()
    require(sha(content) == PINNED_GPL1261_SHA256,
            "Official NCBI GPL1261 annotation bytes changed from qualified snapshot")
    archive = out / "GPL1261.annot.gz"
    archive.write_bytes(content)
    try:
        text = gzip.decompress(content).decode("utf-8")
    except (UnicodeDecodeError, OSError) as exc:
        raise ValueError("EXTERNAL-GEO-BLOCKED: Invalid official gzip/platform text") from exc
    report = parse_gpl_annotation(io.StringIO(text), probes)
    symbol_map = report.pop("probe_to_symbols")
    reverse = report.pop("gene_to_probes")
    write_csv(out / "GPL1261_probe_to_gene_all.csv",
              ["platform_probe_id", "gene_symbol", "probe_mapping_status"],
              (
                  [probe, sym, ("multiple_genes" if len(symbols) > 1 else "single_gene")]
                  for probe, symbols in sorted(symbol_map.items())
                  for sym in (symbols if symbols else [""])
              ))
    write_csv(out / "GPL1261_gene_to_probe_multiplicity.csv",
              ["gene_symbol", "n_platform_probes", "probe_ids"],
              ([gene, len(ids), "|".join(sorted(ids))]
               for gene, ids in sorted(reverse.items())))
    unmapped = report.pop("matrix_probes_without_annotation")
    missing_symbols = report.pop("matrix_probes_with_no_unique_parseable_gene_symbols")
    ambiguous = report.pop("multi_gene_ambiguous_probe_ids")
    write_csv(out / "GPL1261_unmapped_probe_ids.csv", ["platform_probe_id", "reason"],
              ([p, "missing_annotation_row"] for p in unmapped))
    write_csv(out / "GPL1261_unmapped_symbol_ids.csv", ["platform_probe_id", "reason"],
              ([p, "missing_or_unparseable_gene_symbol"] for p in missing_symbols))
    write_csv(out / "GPL1261_ambiguous_probe_ids.csv", ["platform_probe_id", "all_symbols"],
              ([p, "|".join(symbol_map[p])] for p in ambiguous))
    # Cyp3a11 is a PRESPECIFIED Nutrimouse anchor; do not claim equivalence from a name.
    anchor_probes = sorted({p for gene, ids in reverse.items()
                            if gene.lower() == "cyp3a11" for p in ids})
    report.update({
        "official_annotation_source": ANNOT_URL,
        "official_annotation_gzip_sha256": sha(content),
        "official_annotation_gzip_bytes": len(content),
        "external_cohort_source_sha256": PINNED_GEO_SOURCE_SHA256,
        "probes_without_annotation_count": len(unmapped),
        "probes_without_valid_symbol_count": len(missing_symbols),
        "multi_gene_ambiguous_probe_count": len(ambiguous),
        "cyp3a11_anchor_probes": anchor_probes,
        "cyp3a11_status": (
            "probe_name_candidate_only_requires_mapping_and_scale_bridge"
            if anchor_probes else "gene_symbol_not_found_in_official_platform_annotation"
        ),
        "cross_platform_score_projection": "BLOCKED_SCALE_AND_PROBE_COLLAPSE_NOT_VALIDATED",
        "per_animal_scores_generated": 0,
        "clinical_or_PK_validation": False,
        "limitations": [
            "This is a versioned GEO platform annotation source, not a validated human CYP phenotype.",
            "Multiple probes per gene and one probe mapping multiple genes are NOT silently collapsed.",
            "GPL1261 processed RMA and Nutrimouse's original gene measurements lack a qualified frozen scale bridge.",
            "Biological PPARalpha genotype/dose labels are NOT used for selecting probes or mappings.",
            "No external biological score is emitted, and no PK outcome is measured.",
        ],
    })
    (out / "GPL1261_annotation_mapping_audit.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8"
    )
    require(report["per_animal_scores_generated"] == 0,
            "No scores may be manufactured from mere gene symbol mappings")
    print(
        "REAL GEO GPL1261 annotation PASS: "
        f"{report['official_annotation_rows']} source rows; "
        f"{report['matrix_probe_ids_with_official_row']}/45101 original probes matched; "
        f"{report['distinct_annotated_gene_symbols']} distinct gene symbol labels; "
        f"{report['multi_gene_ambiguous_probe_count']} ambiguous probes; "
        f"Cyp3a11 candidate probe count={len(anchor_probes)}; "
        f"annotation_sha256={report['official_annotation_gzip_sha256']}; "
        "30-animal biological external projection BLOCKED."
    )


def selftest():
    fixture = [
        "^PLATFORM = GPL1261\n",
        "!platform_table_begin\n",
        "ID\tGene symbol\tENTREZ_GENE_ID\n",
        "probe_A\tCyp3a11\t1\n",
        "probe_B\tOther1 /// Other2\t2\n",
        "probe_C\t---\t\n",
        "!platform_table_end\n",
    ]
    result = parse_gpl_annotation(fixture, {"probe_A", "probe_B", "probe_C"},
                                  min_source_rows=3)
    assert result["official_annotation_rows"] == 3
    assert result["probe_to_symbols"]["probe_A"] == ["Cyp3a11"]
    assert result["probe_to_symbols"]["probe_B"] == ["Other1", "Other2"]
    assert result["multi_gene_ambiguous_probe_ids"] == ["probe_B"]
    assert result["matrix_probes_with_no_unique_parseable_gene_symbols"] == ["probe_C"]
    assert parse_symbols("Cyp3a11 /// Cyp3a11") == ["Cyp3a11"]
    assert parse_symbols("---") == []
    for bad in (
        [r.replace("Gene symbol", "NO_GENE_IDENTIFIER") for r in fixture],
        [r.replace("!platform_table_end", "!missing") for r in fixture],
        [r.replace("probe_C\t---", "probe_A\t---") for r in fixture],
    ):
        try:
            parse_gpl_annotation(bad, {"probe_A", "probe_B", "probe_C"}, min_source_rows=3)
        except ValueError:
            pass
        else:
            raise AssertionError("GEO annotation parser accepted ambiguous or incomplete input")
    print("PASS: official GPL probe→gene parsing, all ambiguous mappings preserved, missing genes flagged, broken official table rejected")


if __name__ == "__main__":
    try:
        if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
            selftest()
        elif len(sys.argv) == 2:
            qualify(sys.argv[1])
        else:
            raise ValueError("Usage: qualify_gpl1261_annotation.py --self-test | OUTPUT_DIR")
    except Exception as exc:
        sys.exit(str(exc))
