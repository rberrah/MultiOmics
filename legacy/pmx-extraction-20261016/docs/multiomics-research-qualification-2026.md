# Qualification scientifique MultiOmics — 12 octobre 2026

Ce dossier sépare **ce qui est effectivement implémenté et testé**, ce qui
relève de la **validation biologique externe**, et ce qui reste **non qualifié**.
Réussir un workflow GitHub Actions ne certifie pas une nouvelle étude.

## 1. Réplication publique indépendante et provenance

Workflow : `.github/workflows/multiomics-public-benchmark.yml`

La suite charge plusieurs jeux de données publics épinglés à des commits Git
précis et exécute les hypothèses biologiques de référence déjà documentées.
Elle couvre des **sources distinctes** :

- **Nutrimouse** (mixOmics) : ARN et lipides, signal lié à PPARalpha ;
- **TCGA sein** (mixOmics) : ARN et protéines, direction HER2 ;
- **AgingHFCD** (xOmicsShiny) : ARN, protéines et métabolites ; concordance de
  signes contre les sorties de référence publiées ;
- **STATegra** (PaintOmics) : mesures biologiques répétées du métabolome ;
  plan de 36 échantillons et cohérence des directions ;
- **NCI-60** (IntLIM) : relation différentielle gène–métabolite.

Le nouveau `scripts/audit_public_omics_evidence.mjs` **bloque** en cas de jeu
attendu absent, provenance Git divergente, fichier exporté non haché, absence
des contrôles négatifs ou récupération insuffisante des repères préspécifiés.
La sortie `qualification-scope.json` énonce les limites de cette preuve.

**Portée réelle :** reproduction de signaux connus sur jeux publics traités.
Ce n'est pas la validation d'un modèle prédictif figé entre cohortes, ni une
reconstruction des résultats RNA-seq depuis des FASTQ, ni une certification
clinique. Le nombre d'échantillons et les techniques disponibles diffèrent
selon les jeux. La validation externe d'un biomarqueur nécessite une cohorte
indépendante, une hypothèse et un modèle gelés et des mesures de calibration.

## 2. Méthodes et périmètre du plan expérimental

| Plan | Ce qui est disponible | Conclusion confirmatoire automatique |
| --- | --- | --- |
| Deux groupes indépendants | Contrôles d'identifiants, groupes, lots, moteurs R par omique | **Non** : éligibilité technique seulement, puis revue humaine |
| Apparié | Test descriptif/écran apparié, liens sujets | **Non** : revue du modèle sujet/contraste |
| Répété / longitudinal | Parcours R longitudinal limité et diagnostics spécifiques | **Non** : dépendance intra-sujet et hypothèses à vérifier |
| Crossover | Bloqué en l'absence d'un modèle période/séquence/carryover | **Non** |
| Cohortes partiellement appariées | Intersection des vrais sujets et diagnostics de couverture | **Non** pour une intégration qui suppose de faux appariements |

Un panneau FR/EN dans le guide de plan d'étude expose ces frontières, sans
ajouter d'option aux questions de démarrage. Le contrôle R strict des
comptages RNA-seq reste indépendant de cette aide : un panneau réussi n'est
pas une autorisation d'inférence.

## 3. Rapport de méthodes reproductible pour un chercheur

La page propose désormais, dans **Autres exports**, un fichier JSON
`multiomics_methods_scientific_audit.json`, également inclus dans le rapport
HTML complet. Il rassemble :

- objectif, design et recouvrement des échantillons ;
- effectif des sujets, méthode indiquée pour chaque omique et statut des
  moteurs R de référence ;
- nombres de variables mesurées, retenues et exclues lorsque vérifiables ;
- niveau de manquants, avertissements sur l'imputation MNAR ;
- statut d'exécution et limites de validité.

**Confidentialité :** le JSON de méthodes ne contient pas la matrice de
comptages, la table individuelle des patients, ni les prédictions. Pour
réanalyser à l'identique, utiliser séparément le **projet portable** qui
archive les entrées textuelles avec SHA-256, le protocole et les sorties.
Cet archivage peut contenir des informations sensibles : il ne faut pas le
partager publiquement.

**Limite reproductibilité :** ni le JSON de méthodes, ni le projet portable ne
contiennent un verrouillage complet `renv`/conteneur, les RAW instrumentaux,
ou des copies figées de KEGG/Reactome. Le statut scientifique n'est jamais
« publication approuvée ».

## 4. Audit de fiabilité pour les autres omiques

Un script professionnel facultatif analyse les **aires de pics MS** sans
imputer ni sélectionner automatiquement de variables :

```bash
Rscript multiomics-engine/run_ms_quantitative_assurance.R \
  peak_areas.csv sample_manifest.csv output_directory
```

Le fichier `peak_areas.csv` commence par `feature_id`, puis contient une
colonne par `assay_id`. Le manifeste CSV comporte au minimum :
`assay_id,sample_type,condition,injection_order` et, si connu, `batch`.
Les valeurs `sample_type` sont exclusivement `biological`, `pooled_qc`
ou `blank`. Les identifiants doivent correspondre exactement aux colonnes
de la matrice, sans disparition silencieuse d'un échantillon.

Le rapport produit :

- taux de données manquantes par variable biologique ;
- RSD des pooled-QC, quand trois pooled-QC sont disponibles ;
- signal des blancs / médiane biologique, avec au moins deux blancs ;
- déséquilibre de manquants entre deux groupes déclarés ;
- dérive descriptive sur l'ordre d'injection, avec au moins cinq pooled-QC ;
- alerte si le lot technique est totalement confondu avec le groupe.

### Extension de contrôle : preuve mesurable par variable

Le contrôle MS fournit également les **nombres de mesures réellement disponibles
pour chaque variable** dans les échantillons biologiques, pooled-QC et blancs.
La présence de trois injections pooled-QC ne suffit pas si la variable
concernée est absente de deux d'entre elles : un RSD non calculable est
étiqueté `flag_qc_not_evaluable`, jamais « conforme ». De même, si deux
blancs sont déclarés mais ne mesurent pas cette variable, son bruit de fond
est `flag_blank_not_evaluable`.

Cas critique corrigé : si la médiane biologique d'une variable est zéro,
mais le signal des blancs est strictement positif, le rapport émet
`positive_blank_with_zero_biological_median` et `flag_blank_signal`.
Auparavant, un dénominateur nul pouvait seulement produire un ratio manquant,
sans alerte de contamination. Si les deux médianes sont nulles, le rapport
reste **non évaluable** pour ce ratio : absence de signal ne prouve pas absence
de contamination.

Ces statuts sont exportés dans `ms_quantitative_qc_features.csv`, avec
la colonne `pooled_qc_precision_status`, `blank_background_status`
et `pooled_qc_drift_status`. Ils ne calculent pas automatiquement
une correction des données.

Le rapport de méthodes téléchargé depuis le navigateur corrige également
un biais de synthèse : si la médiane des fractions manquantes vaut zéro mais
qu'un échantillon possède des observations manquantes, l'alerte reste
présente. La sortie publique rapporte uniquement le nombre de profils
concernés et leur taux maximal, jamais leurs identifiants.

Ces seuils (RSD 30 %, signal des blancs 20 %, manquants 20 %, déséquilibre
30 points) sont des **écrans de revue par défaut**, à justifier selon le
protocole instrumental ; ils ne prouvent ni mécanisme MCAR/MAR/MNAR ni
validité de l'identification moléculaire. L'absence de blancs ou de pooled-QC
donne une **preuve manquante**, jamais un passage silencieux. Les zéros et
les valeurs manquantes sont préservés. Aucun résultat statistiquement
« intéressant » ne déclenche un changement de prétraitement.

Un audit R complémentaire est disponible dans
`multiomics-engine/proteomics_quantitative_assurance.R`. Il distingue
explicitement `lfq_ms`, `tmt` et `affinity`, oblige à déclarer l'échelle
(`raw` ou `log2`), et vérifie que les sujets biologiques ne sont pas
comptés plusieurs fois. Il conserve, par protéine, le pourcentage de
valeurs manquantes et son déséquilibre entre groupes ; un lot confondu
est signalé. Les limitations propres à chaque technologie sont distinctes :
inférence peptide–protéine et censure pour LFQ-MS, canal de référence et
compression des ratios pour TMT, spécificité de liaison et version du
panel pour affinité. Cet audit ne procède **ni à l'imputation**, **ni à
l'inférence différentielle**, **ni à l'identification des protéines** et
ne peut qualifier une étude sans contrôle des métadonnées.

**Protéomique :** les données LFQ-MS, iTRAQ/TMT et les plateformes par
affinité ont des contrôles de mesure différents ; on ne peut appliquer les
seuils du spectre LC-MS à toutes les technologies. Les associations de
protéines restent conditionnées à la normalisation, aux lots, à la qualité des
identifiants et aux comparaisons indépendantes.

## Test réel du site publié : ne pas confondre HTML et outil fonctionnel

Le contrôle de disponibilité GitHub Pages vérifie les routes publiques, la
**fin effective de l'hydratation Svelte**, le lancement d'une démonstration
avec une seule omique et la présence des résultats, sur ordinateur puis sur
mobile. L'affichage du bouton par le HTML pré-rendu ne prouve pas que son
gestionnaire JavaScript est prêt. Les ressources HTTP et les erreurs navigateur
sont consignées lorsque la démonstration ne démarre pas.

Ce test doit être exécuté **après la réussite du workflow de déploiement** :
tester la version publique avant de déployer une PR aurait comparé le nouveau
code à l'ancienne version visible. Une courte période de propagation du CDN
est tolérée par trois essais bornés ; trois échecs restent bloquants et ne
doivent jamais être présentés comme une réussite scientifique. Les branches
sont testées séparément par leurs propres suites Svelte/Playwright locales.

## Prochain critère non négociable avant usage confirmatoire large

Une validation **indépendante, préenregistrée, sur des données brutes ou de
référence clinique adaptées à chaque technologie**, avec étude de taille
suffisante, scénarios négatifs, covariables et critères gelés. Une interface
simple facilite la réalisation ; elle ne remplace pas l'examen du plan
expérimental par un statisticien ou un biologiste.
