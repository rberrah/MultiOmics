#!/usr/bin/env python3
"""Independently audit mapping of REAL frozen Nutrimouse genes to REAL GEO probes.

No target phenotype, expression-scale fitting, model transport or PK inference.
Passing identifiers is NECESSARY but NOT SUFFICIENT for cross-platform scoring.
Only original immutable Nutrimouse R fits and genuine NCBI/Bioc references qualify.
"""
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

from qualify_geo_gse22871 import PINNED_GEO_SOURCE_SHA256, require
from qualify_gpl1261_annotation import PINNED_GPL1261_SHA256

PINNED_NUTRIMOUSE_COMMIT = "ef3e760526623d9e91e945e4b50be48d764efc40"
TOTAL_GEO_PROBES = 45101


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stable_ids(value):
    if not value:
        return set()
    xs = value.split("|")
    require(all(re.fullmatch(r"[1-9][0-9]*", x) for x in xs),
            "Only literal nonzero numeric Entrez IDs are admitted")
    return set(xs)


def symbols(value):
    return {token.casefold() for token in re.split(r"[,;|\s/]+", value or "") if token}


def probe_candidates(gene, entrez, crosswalk):
    """Strict two-reference matching; do not pick best-correlated target probe."""
    matched = []
    disagreements = Counter()
    for row in crosswalk:
        if (gene.casefold() not in symbols(row["geo_symbols_existing_audit"])
                and gene.casefold() not in symbols(row["bioc_symbols_existing_audit"])):
            continue
        if row["entrez_concordance"] != "exact_entrez_set":
            disagreements["external_entrez_disagreement"] += 1
            continue
        if row["multi_entrez_ambiguous"].casefold() != "false":
            disagreements["multi_gene_probe"] += 1
            continue
        if (stable_ids(row["geo_entrez_ids"]) != {entrez}
                or stable_ids(row["bioc_entrez_ids"]) != {entrez}):
            disagreements["different_entrez_than_source"] += 1
            continue
        if (gene.casefold() not in symbols(row["geo_symbols_existing_audit"])
                or gene.casefold() not in symbols(row["bioc_symbols_existing_audit"])):
            disagreements["reference_gene_symbol_discordance"] += 1
            continue
        matched.append(row["probe_id"])
    return sorted(set(matched)), dict(sorted(disagreements.items()))



def near_source_candidates(gene, entrez_ids, crosswalk, accepted):
    """Diagnostics ONLY: rejected associations are never used to pick a probe."""
    out = []
    for row in crosswalk:
        gs = gene.casefold() in symbols(row["geo_symbols_existing_audit"])
        bs = gene.casefold() in symbols(row["bioc_symbols_existing_audit"])
        # Diagnostics span ALL 45k probes, even unrelated invalid/stale IDs.
        # An unrelated malformed source token is reported in the full
        # concordance archive rather than interrupting this nonselective scan.
        # Only the strictly accepted probe path uses fail-closed stable_ids().
        def diagnostic_ids(value):
            return {v for v in (value or "").split("|")
                    if re.fullmatch(r"[1-9][0-9]*", v)}
        gi = bool(entrez_ids & diagnostic_ids(row["geo_entrez_ids"]))
        bi = bool(entrez_ids & diagnostic_ids(row["bioc_entrez_ids"]))
        if not (gs or bs or gi or bi):
            continue
        out.append({
            "source_gene": gene,
            "source_entrez": "|".join(sorted(entrez_ids)),
            "probe_id": row["probe_id"],
            "geo_symbol": row["geo_symbols_existing_audit"],
            "bioc_symbol": row["bioc_symbols_existing_audit"],
            "geo_entrez_ids": row["geo_entrez_ids"],
            "bioc_entrez_ids": row["bioc_entrez_ids"],
            "entrez_concordance": row["entrez_concordance"],
            "reference_symbol_match": f"NCBI={gs};Bioc={bs}",
            "reference_id_match": f"NCBI={gi};Bioc={bi}",
            "multi_gene": row["multi_entrez_ambiguous"],
            "accepted_as_unique_strict_probe": "YES" if row["probe_id"] in accepted else "NO",
            "diagnostic_only_no_probe_rescue": "TRUE",
        })
    return sorted(out, key=lambda x:x["probe_id"])

def qualify(model_dir, source_ref_dir, geo_ref_dir, out_dir):
    model_dir, source_ref_dir, geo_ref_dir, out_dir = map(
        Path, (model_dir, source_ref_dir, geo_ref_dir, out_dir))
    model_file = model_dir / "frozen_panel_model.json"
    def_file = model_dir / "panel_feature_definition.json"
    identity_file = source_ref_dir / "frozen_nutrimouse_genes_source_entrez.csv"
    ident_report = json.loads((source_ref_dir / "frozen_source_gene_identity_audit.json").read_text())
    model = json.loads(model_file.read_text())
    feature_def = json.loads(def_file.read_text())
    out_dir.mkdir(parents=True, exist_ok=True)

    require(model["pinned_source_commit"] == PINNED_NUTRIMOUSE_COMMIT and
            ident_report["pinned_source_commit"] == PINNED_NUTRIMOUSE_COMMIT,
            "Nutrimouse reference commit drift")
    require(ident_report["model_sha256"] == digest(model_file) and
            ident_report["panel_definition_sha256"] == digest(def_file),
            "Frozen R model SHA mismatch against independent original source annotation")
    require(model["source_rda_sha256"] == ident_report["source_rda_sha256"],
            "Frozen panel does not come from independent real source bytes")
    require(len(model["training_subject_ids"]) == 30 and
            len(model["heldout_subject_ids"]) == 10 and
            not (set(model["training_subject_ids"]) & set(model["heldout_subject_ids"])) and
            model["phenotype_labels_used_for_selection_or_fitting"] is False and
            model["heldout_values_used_for_selection_or_fitting"] is False,
            "Frozen model violates subject-level train/test fit")
    genes = model["feature_ids"]
    require(isinstance(genes, list) and len(genes) == 4 and
            genes[0].casefold() == "cyp3a11" and
            len({x.casefold() for x in genes}) == 4 and
            genes == feature_def["selected_gene_ids"] ==
            ident_report["selected_feature_names"],
            "Frozen original training-selected panel identities changed")
    require(ident_report["features_selected_training_only"] and
            ident_report["source_mapping_is_pheno_blind"],
            "Annotation cannot depend on external or training phenotype")

    with identity_file.open(encoding="utf-8", newline="") as f:
        source_ids = list(csv.DictReader(f))
    require([x["gene_symbol"] for x in source_ids] == genes and
            len(source_ids) == 4, "Missing or reordered source panel annotation")
    require(all(x["source_annotation_version"] ==
                ident_report["source_annotation_version"] for x in source_ids),
            "Genuine org.Mm.eg.db source annotation changed")
    # Verify external publication's ORIGINAL accession, not the phenotype:
    # Martin et al. Hepatology 45 (2007) Table 4 p772: FAS -> NM_007988
    # NCBI mouse RefSeq NM_007988 uniquely resolves to Fasn/Entrez 14104.
    fas = [x for x in source_ids if x["gene_symbol"].casefold()=="fas"]
    require(len(fas)==1 and fas[0]["source_canonical_symbol"]=="Fasn"
            and fas[0]["source_entrez_ids"]=="14104"
            and fas[0]["original_primary_refseq"]=="NM_007988"
            and "DOI_10.1002/hep.21510" in fas[0]["original_primary_refseq_citation"]
            and ident_report["primary_FAS_refseq_verified_against_orgMmEgDb"] is True
            and ident_report["unresolved_historical_FAS_acronym"] is False,
            "Primary original FAS->Fasn accession not verified; reject false Fas/CD95 14102")
    require(len({x["source_canonical_symbol"].casefold() for x in source_ids})==4,
            "Alias resolution collapsed multiple frozen model features")

    geo_audit = json.loads((geo_ref_dir / "GPL1261_entrez_id_concordance_audit.json").read_text())
    require(geo_audit["source_original_gpl1261_sha256"] == PINNED_GPL1261_SHA256 and
            geo_audit["source_original_geo_series_sha256"] == PINNED_GEO_SOURCE_SHA256 and
            geo_audit["source_probe_measurements"] == TOTAL_GEO_PROBES and
            geo_audit["external_scores_produced"] == 0,
            "External GEO crosswalk was not independently qualified")
    path = geo_ref_dir / "GPL1261_all_probe_entrez_concordance.csv"
    with path.open(encoding="utf-8", newline="") as f:
        target_rows = list(csv.DictReader(f))
    require(len(target_rows) == TOTAL_GEO_PROBES and
            len({row["probe_id"] for row in target_rows}) == TOTAL_GEO_PROBES,
            "Incomplete/duplicated actual 45,101 target probes")

    candidates = []
    near_rows = []
    for x in source_ids:
        gene = x["gene_symbol"]  # Frozen model's ORIGINAL macroarray column.
        canonical = x["source_canonical_symbol"]  # Resolved ID for matching only.
        ids = stable_ids(x["source_entrez_ids"])
        if x["source_annotation_status"] == "unique" and len(ids) == 1:
            matched, discord = probe_candidates(canonical, next(iter(ids)), target_rows)
        else:
            matched, discord = [], {"source_gene_id_missing_or_multi": 1}
        possible = near_source_candidates(canonical, ids, target_rows, set(matched))
        near_rows.extend(possible)
        if len(matched) == 1:
            status = "gene_to_target_probe_strict_unique"
        elif len(matched) > 1:
            status = "multiple_valid_target_probes_no_aggregation_rule"
        else:
            status = "no_unambiguous_two_reference_target_probe"
        candidates.append({
            "frozen_feature_index": x["model_feature_index"],
            "source_gene_symbol": gene,
            "resolved_target_gene_symbol": canonical,
            "primary_original_refseq": x["original_primary_refseq"],
            "source_entrez_ids": x["source_entrez_ids"],
            "source_annotation_status": x["source_annotation_status"],
            "matching_geo_probe_ids": "|".join(matched),
            "strict_probe_count": len(matched),
            "technical_crosswalk_status": status,
            "excluded_mapping_disagreements": json.dumps(discord, sort_keys=True),
            "trained_zscore_direct_transport_permitted": "NO",
            "clinically_validated": "NO",
        })
    diag_file = out_dir / "all_possible_source_gene_GEO_probe_matches_no_selection.csv"
    require(bool(near_rows), "Even source CYP3A11 has no original GEO/Entrez evidence")
    with diag_file.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(near_rows[0]))
        w.writeheader()
        w.writerows(near_rows)
    # Gate on all four without optimizing any probe based on external outcomes.
    all_unique = all(x["strict_probe_count"] == 1 for x in candidates)
    rows_file = out_dir / "locked_panel_to_external_probe_crosswalk.csv"
    with rows_file.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(candidates[0]))
        writer.writeheader()
        writer.writerows(candidates)
    report = {
        "source_study": "Nutrimouse",
        "target_study": "GSE22871",
        "frozen_original_nutrimouse_commit": PINNED_NUTRIMOUSE_COMMIT,
        "frozen_original_nutrimouse_data_sha256": model["source_rda_sha256"],
        "source_train_only_model_sha256": digest(model_file),
        "source_model_definition_sha256": digest(def_file),
        "source_mouse_gene_annotation_version": ident_report["source_annotation_version"],
        "primary_FAS_gene_provenance": "Martin_2007_Table4_FAS_NM_007988_Fasn_Entrez14104",
        "primary_FAS_refseq_verified": True,
        "external_gpl1261_source_sha256": PINNED_GPL1261_SHA256,
        "external_original_gse22871_sha256": PINNED_GEO_SOURCE_SHA256,
        "n_source_train_mice": 30,
        "n_source_holdout_mice": 10,
        "n_independent_target_mice": 30,
        "panel_n_frozen_genes": 4,
        "strictly_mapped_feature_count": sum(x["strict_probe_count"] == 1 for x in candidates),
        "all_four_panel_genes_strictly_mapped": all_unique,
        "source_model_scale": "training_zscore_of_preprocessed_macroarray_expression",
        "target_data_scale": "log2_RMA_Affymetrix_GPL1261",
        "direct_source_zscore_to_target_is_justified": False,
        "cross_platform_reference_bridge_fitted": False,
        "overlapping_biological_replicate_reference_samples": 0,
        "external_phenotype_used_for_probe_matching": False,
        "n_external_signature_scores_exported": 0,
        "n_PK_concentration_outcomes": 0,
        "transport_status": ("IDENTIFIERS_MATCH_SCALE_UNCALIBRATED"
                             if all_unique else
                             "BLOCKED_UNMAPPED_OR_AMBIGUOUS_GENES_AND_SCALE"),
        "interpretation": (
            "Real independent identity feasibility only; original Nutrimouse "
            "z-scores cannot be evaluated on Affymetrix RMA log2 with these "
            "training parameters absent a separately validated transfer scale. "
            "No sign or significance of the source or external intervention is "
            "used to select or discard genes/probes."
        ),
        "candidate_table_sha256": digest(rows_file),
        "all_proximity_diagnostics_sha256": digest(diag_file),
        "n_diagnostic_probe_rows_including_rejections": len(near_rows)
    }
    (out_dir / "frozen_panel_external_transport_readiness.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("REAL frozen Nutrimouse 4-gene panel versus GEO original 45101 probes: ",
          report["strictly_mapped_feature_count"],"/4 strictly unique; status=",
          report["transport_status"],"; 0 claimed transferred signature scores; model=",
          report["source_train_only_model_sha256"][:12],sep="")
    for entry in candidates:
        gene=entry["source_gene_symbol"]
        close=[x for x in near_rows if x["source_gene"]==entry["resolved_target_gene_symbol"]]
        print(gene,entry["technical_crosswalk_status"],
              entry["matching_geo_probe_ids"],
              "source_rejected_mapping_reasons=",entry["excluded_mapping_disagreements"],
              "all_potential_probe_rows=",len(close))
        if not entry["strict_probe_count"]:
            for match in close[:12]:
                print("REJECTED_OR_AMBIGUOUS_REFERENCE",match["probe_id"],
                      "NCBIgene=",match["geo_symbol"],
                      "BiocGene=",match["bioc_symbol"],
                      "NCBIentrez=",match["geo_entrez_ids"],
                      "BiocEntrez=",match["bioc_entrez_ids"],
                      "ids=",match["entrez_concordance"],
                      "multigene=",match["multi_gene"])
    return report


def selftest():
    base = {"probe_id": "p1", "geo_entrez_ids": "123",
            "bioc_entrez_ids": "123", "entrez_concordance": "exact_entrez_set",
            "geo_symbols_existing_audit": "Cyp3a11",
            "bioc_symbols_existing_audit": "Cyp3a11",
            "multi_entrez_ambiguous": "False"}
    assert probe_candidates("Cyp3a11", "123", [base])[0] == ["p1"]
    other = dict(base,probe_id="p2")
    assert len(probe_candidates("Cyp3a11", "123", [base,other])[0]) == 2
    assert probe_candidates("Cyp3a11", "124", [base])[0] == []
    assert probe_candidates("Cyp3a11", "123", [dict(base,bioc_entrez_ids="124")])[0] == []
    assert probe_candidates("Cyp3a11", "123", [dict(base,bioc_symbols_existing_audit="Cyp3a110")])[0] == []
    assert probe_candidates("Cyp3a11", "123", [dict(base,multi_entrez_ambiguous="True")])[0] == []
    try:
        stable_ids("123|malicious")
        raise AssertionError("Invalid gene ID accepted")
    except ValueError:
        pass
    assert "cyp3a11" not in symbols("Cyp3a110")
    # Negative control: the source FAS acronym must NOT point at death receptor.
    assert probe_candidates("Fasn","14104",[dict(base,probe_id="x",
        geo_symbols_existing_audit="Fas",
        bioc_symbols_existing_audit="Fas",
        geo_entrez_ids="14102",bioc_entrez_ids="14102")])[0] == []
    print("PASS: independent source gene matching; wrong ID, target disagreement, duplicated/multi probes, malformed gene ID blocked")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        selftest()
    elif len(sys.argv) == 5:
        qualify(*sys.argv[1:])
    else:
        sys.exit("Usage: qualify_nutrimouse_geo_gene_bridge.py --self-test | FROZEN_R_DIR SOURCE_ID_DIR GEO_ENTREZ_DIR OUT_DIR")
