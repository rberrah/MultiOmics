# Origine tissulaire et cellulaire des omiques — contrat MultiOmics → Omics2PK
*Statut : qualification analytique, sans fonction enzymatique inférée ni validation clinique. Applicable aux omiques **bulk** comme aux autres matrices.*

## Pourquoi ?
**La provenance est une partie de la mesure.** Pour CYP3A4/5, un ARN extrait de PBMC ne constitue pas une mesure de CYP3A4/5 hépatique ou intestinal. Un plasma acellulaire ne contient pas les mêmes cellules qu'un échantillon de sang total ; une protéine détectée ne prouve pas une enzyme fonctionnelle. L'expression `bulk` varie avec la composition cellulaire et l'état de chaque sous-population. Sans données sources, l'organe de prélèvement et l'organe contribuant à la PK restent distincts.

## Extension strictement additive de l'export V1
`exchange/samples.csv` conserve désormais les champs sources s'ils sont réellement présents (d'autres champs non connus restent inchangés) :
- **Source :** `tissue` (champ existant), `specimen_type` (e.g. PBMC, serum, plasma, liver_biopsy), `sample_matrix`, `tissue_ontology_id` (UBERON, CL, EFO **réels**).
- **Cellularité :** `cell_type_or_bulk` (`bulk`, `single_cell`, `single_nucleus`, `sorted_cells`, `cell_line`, `cell_free` ou `unknown`), `cell_composition_status` (`measured`, `estimated`, `not_measured`, `unknown`), `cell_composition_method`.
- **Protocole :** `assay_platform`, `source_preprocessing`, `source_protocol_id`.
- **Contexte PK distinct :** `pharmacology_target_organ` (organe contribuant à l'effet étudié, seulement si renseigné), `relative_to_dose` (`pre_dose`, `post_dose`, `at_dose`, `not_applicable`, `unknown`). Un `T0` non documenté n'est **jamais** converti en « prédose ».

Les champs omis restent **vides** (inconnus). Les champs explicites sont recopiés sous un format stable, **sans inventer** une origine hépatique, une ontologie, une cellule ou une heure. Le test de cohérence interdit des valeurs tissulaires incompatibles entre deux analyses issues du **même** `sample_id`. Le schéma reste `1.0.0` ; les champs ajoutés ne modifient pas les valeurs numériques ni les statistiques. Le parseur original d'Omics2PK accepte déjà les colonnes optionnelles.

`exchange/manifest.json` ajoute un **résumé agrégé** `validation.source_location_provenance` : nombres de prélèvements dont la provenance est renseignée, proportions de sources inconnues, nombre d'ARN/plasma circulants dont l'organe cible déclaré est le foie. Ces discordances sont classées `n_cross_organ_proxy_not_validated`, **pas** interdites pour toutes les analyses : elles requièrent une **validation de proxy** sur des donneurs correctement appariés.

L'analyse, les tests et l'export **n'attribuent jamais** de fonction enzyme/transporteur, d'identité cellulaire déduite d'un gène, ou de transfert prouvé sang→foie. Les essais in vitro et in vivo doivent conserver l'origine biologique et les conditions d'administration.

## Tests
`node scripts/test_multiomics_tissue_provenance.mjs` couvre :
- ARN PBMC + organe cible foie : **proxy non validé**, pas « enzyme CYP3A hépatique ».
- Métabolomique de plasma acellulaire : origine différente de PBMC et sang total.
- ARN hépatique de biopsie : tissu source explicitement attesté.
- Valeurs et identifiants de prélèvement source exportés sans altération.
- Origine inconnue reste inconnue ; impossible d'inférer « prédose » d'un seul `T0`.
- Même `sample_id` mais deux origines contradictoires : refus avant export.

## Prochaine couche scientifique
La cartographie **Reactome/KEGG** et l'intégration MultiOmics doivent annoter chaque réaction et variable avec l'organe réellement échantillonné, les identifiants chimiques/géniques et leur assurance analytique ; sans lien vérifié, afficher `PROXY_NOT_VALIDATED` ou `TARGET_TISSUE_MISMATCH`. La comparaison clinique effective appartient à **Omics2PK** sur un endpoint PK réel, avec un modèle de référence clinique + PGx et validation externe.

Sources : [GTEx open public expression](https://gtexportal.org/home/downloads/adult-gtex/), [GTEx tissue-specific effects](https://doi.org/10.1126/science.aaz1776) ; [Metformin/Rotroff 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4906013/), source de PK observée.

**Note confidentialité :** les exports peuvent contenir des données humaines. Ne pas les publier sans autorisation.
