// Offline and dependency-free verifier packaged alongside an analysis dossier.
// Runs on the investigator's own computer; it reads but NEVER extracts or
// transmits potentially sensitive sample data.
export const RESEARCH_DOSSIER_VERIFY_PY = String.raw`#!/usr/bin/env python3
"""Check a MultiOmics research ZIP without extraction or network access.

Usage: python3 verifier_dossier.py multiomics_dossier_recherche_v1.zip

Checks: archive member integrity, SHA-256, exact original input hashes,
declared layer table completeness, absence of archived "fresh" results.
This checks FILE INTEGRITY ONLY; it does not validate statistics.
"""
import csv
import hashlib
import io
import json
import pathlib
import sys
import zipfile


def failure(message):
    raise ValueError("DOSSIER INVALIDE / INVALID DOSSIER: " + message)


def verify(archive_path):
    with zipfile.ZipFile(archive_path, "r") as z:
        infos = z.infolist()
        names = [i.filename for i in infos]
        if len(names) != len(set(names)):
            failure("duplicate archive members")
        if any(n.startswith("/") or ".." in pathlib.PurePosixPath(n).parts for n in names):
            failure("unsafe archive member path")
        if "provenance/manifest.json" not in names:
            failure("missing provenance/manifest.json")
        if z.testzip() is not None:
            failure("ZIP CRC-32 mismatch")
        m = json.loads(z.read("provenance/manifest.json"))
        if m.get("schema") != "pmx-multiomics-research-dossier/v1":
            failure("unknown research schema")
        if m.get("confirmatory_certification") is not False:
            failure("unexpected scientific certification")
        if m.get("clinical_validation") is not False:
            failure("unexpected clinical certification")
        indexed = m.get("files")
        if not isinstance(indexed, dict):
            failure("missing file digests")
        expected = set(indexed) | {"provenance/manifest.json"}
        if set(names) != expected:
            failure("missing or unaccounted archive members")
        for name, record in indexed.items():
            raw = z.read(name)
            if hashlib.sha256(raw).hexdigest() != record.get("sha256"):
                failure("SHA-256 differs: " + name)
            if len(raw) != record.get("bytes"):
                failure("file size differs: " + name)
        project = json.loads(z.read("source/projet_reproductible.json"))
        if project.get("format") != "pmx-multiomics-portable-project":
            failure("project format differs")
        if project.get("archivedResult") is not None:
            failure("archived result incorrectly treated as newly executed")
        for key, source in project.get("inputs", {}).items():
            digest = hashlib.sha256(source["content"].encode("utf-8")).hexdigest()
            if digest != source["sha256"]:
                failure("original input SHA-256 differs: " + key)
            if digest != m.get("input_file_sha256", {}).get(key):
                failure("manifest input SHA-256 differs: " + key)
        details = json.loads(z.read("provenance/analysis_complete.json"))
        expected_layers = m.get("browser_layer_rows")
        if not isinstance(expected_layers, dict) or not expected_layers:
            failure("missing measured layer counts")
        if set(expected_layers) != set(details.get("layers", {})):
            failure("layer types differ from full analysis")
        for layer, expected_count in expected_layers.items():
            if expected_count != len(details["layers"][layer].get("rows", [])):
                failure("declared result count differs: " + layer)
            csv_name = "tableaux/" + layer + ".csv"
            rows = list(csv.DictReader(io.StringIO(
                z.read(csv_name).decode("utf-8-sig"), newline="")))
            if len(rows) != expected_count:
                failure("CSV table truncated: " + layer)
            if not all(row.get("feature") for row in rows):
                failure("missing feature identifier: " + layer)
        ref = m.get("r_reference_full_tables", [])
        if not isinstance(ref, list):
            failure("invalid R reference manifest")
        if bool(ref) != (m.get("r_reference_full_feature_tables_exported") is True):
            failure("R reference completeness flag differs")
        seen_ref = set()
        for fit in ref:
            filename = fit.get("name")
            if not isinstance(filename, str) or not filename.startswith("reference_r/"):
                failure("invalid R reference member name")
            if filename in seen_ref or filename not in names:
                failure("duplicate or missing R reference CSV")
            seen_ref.add(filename)
            values = list(csv.DictReader(io.StringIO(
                z.read(filename).decode("utf-8-sig"), newline="")))
            if len(values) != fit.get("rows"):
                failure("R fitted feature rows truncated: " + filename)
            identifiers = [row.get("feature") for row in values]
            if any(not x for x in identifiers) or len(set(identifiers)) != len(identifiers):
                failure("nonunique or missing R reference feature identifiers")
        print("VERIFICATION OK / INTEGRITY PASS: %d layers; %d members"
              % (len(expected_layers), len(names)))
        print("NOTE: this proves file integrity, NOT statistical or clinical validity")


if __name__ == "__main__":
    try:
        if len(sys.argv) != 2:
            failure("usage: python3 verifier_dossier.py analysis.zip")
        verify(sys.argv[1])
    except (OSError, KeyError, ValueError, UnicodeDecodeError,
            zipfile.BadZipFile) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
`;

export const RESEARCH_DOSSIER_README = `MultiOmics / PMx Explain - dossier de recherche / research dossier

1. Dans le ZIP, les fichiers tableaux/*.csv contiennent TOUS les résultats
   exploratoires du navigateur. Ils ne sont pas automatiquement des résultats
   d'un moteur R de référence.
2. methods/methods.md explique les méthodes, les colonnes et les limites.
3. source/projet_reproductible.json contient les fichiers ORIGINAUX (données
   potentiellement personnelles ou sensibles). Ne pas le diffuser sans accord.
4. provenance/analysis_complete.json garde les sorties intégrales et les
   calculs R éventuels. provenance/manifest.json liste les SHA-256.
5. Vérification autonome (Python 3, aucune dépendance ni connexion réseau):
      python3 verifier_dossier.py multiomics_dossier_recherche_v1.zip
   Extraire verifier_dossier.py puis lancer la commande sur le ZIP intact.
6. Le vérificateur contrôle l'intégrité des fichiers et l'exhaustivité des
   tableaux, pas la validité biologique, statistique ou clinique.
7. Pour reproduire une analyse: importez le projet JSON dans MultiOmics,
   revérifiez la configuration et RELANCEZ l'analyse. Un résultat archivé
   n'est jamais présenté comme un nouveau calcul.
8. Conservez séparément les versions exactes du backend R, des bases/API,
   et les fichiers MS binaires non supportés par l'import texte.

English:
Read methods/methods.md, use CSVs for full browser results, preserve original
JSON data and all SHA-256 provenance. Python verification requires no packages
and makes no network request. R results, when present, are stored separately
in analysis_complete.json. This archive does NOT certify validity.
`;
