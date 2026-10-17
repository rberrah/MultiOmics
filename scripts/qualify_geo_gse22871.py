#!/usr/bin/env python3
"""Fail-closed external-cohort intake for authentic GEO GSE22871.

This script validates the original, measured processed Affymetrix matrix and
its 30 biological sample records. It deliberately DOES NOT project the
Nutrimouse gene-symbol model onto Affymetrix probe IDs.
"""
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sys
import tempfile
import urllib.request
from collections import Counter
from pathlib import Path

ACCESSION = "GSE22871"
PLATFORM = "GPL1261"
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/"
       "GSE22871/matrix/GSE22871_series_matrix.txt.gz")
GENOTYPES = ("PPARalpha-null", "wild-type")
DOSES = (0, 3, 10)
SAMPLE_IDS = tuple("GSM" + str(i) for i in range(565106, 565136))
EXPECTED_PROBES_MIN = 40000
EXPECTED_PROBES_MAX = 50000
MAX_GZIP_BYTES = 50_000_000


def require(ok, message):
    if not ok:
        raise ValueError("EXTERNAL-GEO-BLOCKED: " + message)


def parse_title(title):
    """Read only documented source metadata, never infer genotype from GSM order."""
    genotype = next((g for g in GENOTYPES if title.lower().startswith(g.lower() + " mouse")), None)
    match = re.search(r"\+\s*(0|3|10)\s*mg/kg\s*PFOS,\s*biological\s+rep([1-5])\s*$",
                      title, flags=re.IGNORECASE)
    require(genotype is not None and match is not None,
            "Unknown or ambiguous actual GSM title: " + repr(title))
    return genotype, int(match.group(1)), int(match.group(2))


def split_geo_line(line):
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def parse_matrix(stream, *, min_probes=EXPECTED_PROBES_MIN,
                 max_probes=EXPECTED_PROBES_MAX):
    meta = {}
    sample_ids = []
    observed_probes = set()
    probes = 0
    missing = 0
    sample_count = None
    sums = []
    minima = []
    maxima = []
    in_table = False
    ended = False
    for raw in stream:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if not in_table:
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line.startswith("!Sample_") or line.startswith("!Series_platform_id"):
                fields = split_geo_line(line)
                meta[fields[0]] = fields[1:]
            continue
        if not sample_ids:
            row = split_geo_line(line)
            require(len(row) == 31 and row[0] in ("ID_REF", "ID"),
                    "Matrix header must contain an ID_REF and exactly 30 original GSM assays")
            sample_ids = row[1:]
            require(len(set(sample_ids)) == 30 and set(sample_ids) == set(SAMPLE_IDS),
                    "GEO assay IDs do not match the 30 original GSE22871 GSM IDs")
            sample_count = len(sample_ids)
            sums = [0.0] * sample_count
            minima = [float("inf")] * sample_count
            maxima = [float("-inf")] * sample_count
            continue
        if line == "!series_matrix_table_end":
            ended = True
            break
        row = split_geo_line(line)
        require(len(row) == 31, "Ragged measured probe row in external dataset")
        probe = row[0]
        require(bool(probe) and not probe in observed_probes,
                "Missing or duplicated real Affymetrix probe ID")
        observed_probes.add(probe)
        probes += 1
        require(probes <= max_probes, "Too many probes; source file may have changed")
        for i, raw_value in enumerate(row[1:]):
            if raw_value in ("", "NA", "null", "NaN"):
                missing += 1
                continue
            try:
                value = float(raw_value)
            except ValueError as exc:
                raise ValueError("EXTERNAL-GEO-BLOCKED: non-numeric processed probe value") from exc
            require(math.isfinite(value), "Non-finite measured probe value")
            sums[i] += value
            minima[i] = min(minima[i], value)
            maxima[i] = max(maxima[i], value)
    require(ended and min_probes <= probes <= max_probes,
            "Incomplete series matrix or nonconforming number of original probe IDs")
    require(sample_count == 30, "External processed matrix has wrong biological sample count")
    require(missing == 0, "Preprocessed Affymetrix source has missing expression values")
    platform = meta.get("!Series_platform_id", [])
    require(PLATFORM in platform, "Source GEO platform is not GPL1261")
    orig_gsm = meta.get("!Sample_geo_accession", [])
    orig_title = meta.get("!Sample_title", [])
    require(len(orig_gsm) == 30 and len(orig_title) == 30,
            "Original GSM identities/titles are absent from the study metadata")
    require(orig_gsm == sample_ids,
            "Metadata GSM order differs from source matrix columns")
    groups = Counter()
    repeats = {}
    samples = []
    for i, (gsm, title) in enumerate(zip(orig_gsm, orig_title)):
        genotype, dose, rep = parse_title(title)
        groups[(genotype, dose)] += 1
        repeats.setdefault((genotype, dose), set()).add(rep)
        samples.append({
            "assay_id": gsm, "subject_id": gsm, "sample_id": gsm,
            "source_title": title, "genotype": genotype,
            "dose_mg_per_kg_day": dose, "replicate_label": rep,
            "platform_probe_namespace": PLATFORM,
            "expression_mean": sums[i] / probes,
            "expression_min": minima[i], "expression_max": maxima[i],
        })
    for genotype in GENOTYPES:
        for dose in DOSES:
            require(groups[(genotype, dose)] == 5 and
                    repeats[(genotype, dose)] == {1, 2, 3, 4, 5},
                    f"Unbalanced or duplicated actual animals: {genotype} PFOS={dose}")
    return {
        "dataset": ACCESSION,
        "n_biological_subjects": 30,
        "n_processed_affymetrix_probe_ids": probes,
        "platform_id": PLATFORM,
        "expression_scale": "processed_RMA_log2_microarray_intensity_source_declared",
        "tissue": "liver",
        "source_is_original_GEO_series_matrix": True,
        "source_GSMs_checked_against_GEO": True,
        "genotype_dose_groups": [
            {"genotype": g, "dose_mg_kg_day": d, "n": groups[(g, d)]}
            for g in GENOTYPES for d in DOSES
        ],
        "samples": samples,
        "analysis_status": "intake_qualified_external_projection_blocked",
        "projection": {
            "status": "BLOCKED_UNTIL_ANNOTATION_AND_SCALE_BRIDGE",
            "source_model": "Nutrimouse same-species mouse liver macroarray (120 gene names)",
            "external_data": "GSE22871 Affymetrix GPL1261 probe IDs",
            "reason": (
                "No independent, versioned many-to-many GPL1261 probe-to-gene mapping"
                " or cross-platform frozen measurement scale transfer has yet been validated."
            ),
            "score_rows_exported": 0,
            "independent_PK_outcomes_available": False,
        },
        "limitations": [
            "PFOS exposure and genetic status are experimental design factors, not PK outcomes.",
            "Only a separate public animal cohort intake has been validated.",
            "Original processed values are not original CEL microarray intensities.",
            "Do not combine animals across cohorts for training or select probes by heldout labels.",
            "The 4-gene Nutrimouse panel is a technical training-variance composite, not a biological signature.",
        ],
    }


def download_bytes():
    req = urllib.request.Request(URL, headers={
        "User-Agent": "MultiOmics-research-validation/1.0 (GEO public reproducibility)",
        "Accept": "application/gzip",
    })
    with urllib.request.urlopen(req, timeout=75) as stream:
        chunks = []
        size = 0
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            size += len(chunk)
            require(size <= MAX_GZIP_BYTES, "GEO source exceeds declared compressed size limit")
            chunks.append(chunk)
    data = b"".join(chunks)
    require(data[:2] == b"\x1f\x8b", "Fetched source is not original gzip GEO matrix")
    return data


def run(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    content = download_bytes()
    archive = out / "GSE22871_series_matrix.txt.gz"
    archive.write_bytes(content)
    with gzip.open(archive, "rt", encoding="utf-8", errors="strict",
                   newline="") as stream:
        report = parse_matrix(stream)
    report["source_url"] = URL
    report["source_archive_sha256"] = hashlib.sha256(content).hexdigest()
    report["source_archive_bytes"] = len(content)
    (out / "external_geo_intake_audit.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    fields = ["assay_id", "subject_id", "sample_id", "source_title", "genotype",
              "dose_mg_per_kg_day", "replicate_label", "platform_probe_namespace",
              "expression_mean", "expression_min", "expression_max"]
    with (out / "samples.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(report["samples"])
    require(report["projection"]["score_rows_exported"] == 0,
            "External model projection MUST remain blocked without qualified gene mappings")
    print("REAL GSE22871: 30/30 original animals, six verified groups, "
          f"{report['n_processed_affymetrix_probe_ids']} genuine probe rows "
          f"source_SHA256={report['source_archive_sha256']} "
          "and source data integrity PASS; Nutrimouse external scoring BLOCKED (correct).")


def selftest():
    titles = [
        (g, d, r, f"{g} mouse +{d}mg/kg PFOS, biological rep{r}")
        for g in GENOTYPES for d in DOSES for r in range(1, 6)
    ]
    gsm = list(SAMPLE_IDS)
    def synthetic_lines(items):
        begin = [
            '!Series_platform_id\t"GPL1261"',
            '!Sample_geo_accession\t' + "\t".join(f'"{x}"' for x in gsm),
            '!Sample_title\t' + "\t".join(f'"{x[3]}"' for x in items),
            '!series_matrix_table_begin',
            '"ID_REF"\t' + "\t".join(f'"{x}"' for x in gsm)
        ]
        rows = [f'"Affy_probe_{k}"\t' + "\t".join(str(8 + (k+i)/1000) for i in range(30))
                for k in range(15)]
        return [line + "\n" for line in begin + rows + ["!series_matrix_table_end"]]
    record = parse_matrix(synthetic_lines(titles), min_probes=10, max_probes=20)
    assert record["n_biological_subjects"] == 30
    assert record["projection"]["score_rows_exported"] == 0
    assert len(record["genotype_dose_groups"]) == 6
    wrong = list(titles)
    wrong[0] = titles[1]
    for mutate, label in [
        (synthetic_lines(wrong), "duplicate biological replicate"),
        ([x.replace("GPL1261", "GPL_OTHER") for x in synthetic_lines(titles)], "platform"),
        ([x.replace('"GSM565106"', '"GSM_WRONG"') for x in synthetic_lines(titles)], "assay ID"),
        ([x.replace("8.0", "NaN") for x in synthetic_lines(titles)], "nonfinite"),
    ]:
        try:
            parse_matrix(mutate, min_probes=10, max_probes=20)
        except ValueError:
            continue
        raise AssertionError("Expected fail-closed GEO negative control: " + label)
    print("PASS: GEO authentic-series parser 30 subjects, six groups, no pseudo-replicates,"
          " wrong-platform, changed-GSM and invalid-value refusal.")


if __name__ == "__main__":
    try:
        if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
            selftest()
        elif len(sys.argv) == 2:
            run(sys.argv[1])
        else:
            raise ValueError("Usage: script.py --self-test | OUTPUT_DIR")
    except Exception as exc:
        sys.exit(str(exc))
