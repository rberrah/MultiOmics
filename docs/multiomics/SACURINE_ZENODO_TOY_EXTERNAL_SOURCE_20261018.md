# Nouveau jeu externe LC-MS : Sacurine toy v3 (Zenodo) — intégrité de la source

**Périmètre de cette qualification :** audit des métadonnées et des références MD5 **originales**, pas encore traitement xcms ni validation biologique.

## Pourquoi un nouveau point d'entrée ?

L'archive complète `sacurine-neg-subset_mzXML.zip` de Workflow4Metabolomics a été bloquée par des timeouts externes dans les tentatives documentées du 16 octobre. Il serait incorrect de présenter une prévalidation hors ligne comme une réussite instrumentale.

Un deuxième jeu **officiellement publié** est maintenant identifié : [Sacurine toy dataset, Zenodo v3, DOI 10.5281/zenodo.10130758](https://zenodo.org/records/10130758). Les auteurs sont Thévenot, Roux, Xu et Junot. Contrairement aux tables de pics Drosophila et aux quatre fichiers faahKO, ce jeu comprend **les trois types d'échantillons** et les ordres d'injection initiaux fournis dans la source elle-même.

| Source originale publiée | Nombre | Statut |
|---|---:|---|
| `HU_neg_*.mzML` | 6 | Urines biologiques |
| `QC1_*.mzML` | 3 | QC pools |
| `Blanc*.mzML` | 3 | Blancs analytiques |
| `sampleMetadata_12samp_completed.tsv` | 12 lignes | `injectionOrder`, `batch`, `class`, `sampleType` **réellement documentés** |

L'archive d'origine totalise environ 257,4 Mo en fichiers `.mzML` (et un fichier de métadonnées de 645 octets), plus petit que le sous-ensemble Sacurine W4M de 36 injections.

**Chronologie importante :** les ordres de source publiés vont de 27 à 185 : ce sont des injections espacées dans une acquisition plus grande, **pas** douze passages consécutifs. Les trois QC sont aux ordres 27, 105 et 185 ; les trois blancs aux ordres 29, 94 et 173. Une correction de dérive instrumentale robuste ou une évaluation quantitative du carryover ne peut pas être déclarée réussie simplement parce que ces ordres existent.

## Travail d'intégrité effectué en phase A

Les scripts `scripts/qualify_sacurine_toy_zenodo.py` et `scripts/test_sacurine_zenodo_preflight.py` :

- fixent la **version v3** du dépôt Zenodo (ne pas suivre le « latest » sans nouvel audit) ;
- exigent le fichier `sampleMetadata_12samp_completed.tsv` original et son empreinte MD5 publiée `f88044818627dee797c75d5a310f5c8e` ;
- vérifient les 12 identités de fichiers `.mzML`, leurs 12 MD5 publiés et les rôles précisément **6 biologiques / 3 QC / 3 blancs** ;
- refusent les erreurs `class`/`sampleType`, ordre d'injection absent ou dupliqué, mauvais batch, doublons, champs corrompus, inventaire incomplet et fichiers chromatographiques aux MD5 modifiées ;
- enregistrent explicitement `raw_file_md5_verified=false` si les chromatogrammes ne sont **pas** téléchargés, `raw_peak_processing=NOT_EXECUTED`, `clinical_certification=false` ;
- ne sauvegardent dans le dépôt **aucun fichier d'origine** ou information médicale individuelle.

Le workflow `.github/workflows/sacurine-zenodo-source-preflight.yml` récupère en ligne **le vrai TSV**, et refuse de remplacer ce téléchargement par la fixture synthétique des tests unitaires. Une erreur réseau est un échec de la qualification de source, pas une erreur xcms.

## Phase B exécutée après validation de la phase A

1. Télécharger les **12 fichiers mzML originaux** depuis la même version Zenodo, vérifier chacune des empreintes MD5 de l'éditeur avant tout calcul.
2. Lire les spectres réellement acquis sous `xcms`, figer les paramètres CentWave, l'alignement chromatographique et le regroupement avant interprétation, exporter les aires de pics intégrées (ce ne sont pas des AUC PK).
3. Construire la matrice features×injections avec les **vrais identifiants**, les trois rôles, les vrais ordres et des manquants correctement distingués.
4. Quantifier le niveau de fond des blancs, la proportion de pics présentant un signal dans les pools et le RSD QC, sans revendiquer une correction de dérive validée sur seulement trois pools très espacés.
5. Confronter le résultat à un pipeline scientifique de référence indépendant, avec contrôle de l'instrument, de l'échantillonnage et des métadonnées. Conserver les résultats négatifs.
6. N'extrapoler à aucune analyse PK/PD sans **mesures de concentrations PK appariées individuellement**, absentes de cette preuve analytique.

**État après phase B :** les métadonnées, rôles, identifiants, ordres publiés et les 12 chromatogrammes ont été authentifiés, puis la chaîne xcms a été exécutée et auditée sur cette source réelle. Cela permet de revendiquer une **qualification d'exécution raw-MS sur ces 12 injections**, avec les résultats QC négatifs conservés. Il reste impossible de conclure à une correction validée de dérive/carryover, à une identité métabolique confirmée, à des biomarqueurs cliniquement validés ou à un passage Omics→PK démontré.

Cette qualification **n'affecte pas** les scripts de production du moteur MultiOmics ou l'interface débutant.
