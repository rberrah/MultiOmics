# Guide MultiOmics : deux prélèvements d'un même sujet (RNA-seq apparié)

**À qui s'adresse ce guide ?** Aux chercheurs qui disposent de comptages RNA-seq bruts mesurés dans deux conditions chez les **mêmes sujets biologiques** (patients, animaux ou lignées cellulaires) et souhaitent récupérer des résultats complets, traçables et interprétables.

## 1. Comprendre le plan en 20 secondes

Si les mêmes quatre sujets sont mesurés avant et après traitement, le nombre d'unités biologiques **indépendantes est 4, pas 8**. Le modèle DESeq2 correct pour deux conditions appariées est `~ subject_id + condition` (ou `~ cell + condition` pour des lignées), avec le **groupe de référence indiqué explicitement**. Le sujet est utilisé comme bloc fixe ; l'effet du traitement est estimé *à l'intérieur* des sujets.

Ce parcours **n'est pas adapté** à un crossover multi-périodes, à un suivi longitudinal à plus de deux temps ou à des traitements dont les effets de période ou de séquence doivent être estimés. Ces plans nécessitent d'autres modèles.

## 2. Vos deux fichiers

1. **Tableau des échantillons (.csv).** Une ligne par mesure `assay_id`. Colonnes : `subject_id`, `sample_id`, `assay_id`, `omic`, `condition`, `timepoint`, `batch`, `technical_replicate`, `sample_type`. Pour chaque même sujet, garder **le même `subject_id`**, mais donner des identifiants différents aux deux prélèvements biologiques et aux mesures. `sample_type=biological` pour tous les prélèvements du modèle. Ne pas inscrire les contrôles QC comme s'il s'agissait de sujets.
2. **Matrice RNA-seq (.csv).** Une ligne par gène, première colonne `feature_id` ; une colonne par `assay_id` correspondant exactement au tableau. Importer des **comptages bruts, entiers, non négatifs, non normalisés**. Ne jamais mettre des TPM/FPKM sous l'étiquette « raw_counts ».

Dans MultiOmics, sélectionner **« Deux mesures liées chez les mêmes sujets »** puis télécharger le **modèle fictif de quatre paires**. Remplacer les identifiants fictifs et ne jamais interpréter les valeurs de démonstration comme des données biologiques.

## 3. Ce que le logiciel doit vérifier avant DESeq2

- **Deux conditions exactement**, au moins quatre sujets appariés, et chaque sujet mesuré dans les deux conditions, sans échantillon biologique dupliqué pour une condition.
- Identifiants exacts entre matrice et métadonnées, aucune mesure orpheline ni suppression silencieuse d'échantillons ; les réplicats techniques restent distincts des unités biologiques.
- Comptages valides et matrice de conception **de rang complet** ; une nuisance (batch, covariable) non estimable ou confondue avec la condition ou le sujet empêche l'ajustement.
- Un contraste **groupe comparé versus groupe de référence** sélectionné volontairement en analyse confirmatoire.
- **R requis**, car le calcul exploratoire du navigateur n'est pas un moteur d'inférence RNA-seq appariée.

Le seuil de 4 sujets appariés est une règle minimale de sécurité logicielle, **pas** une garantie de puissance statistique ni de fiabilité biologique.

## 4. Obtenir une analyse reproductible

Le backend R local produit les véritables estimations DESeq2/limma-voom. Après calcul, télécharger le dossier de recherche ZIP ; distinguer :

- `tableaux/` : calculs exploratoires dans le navigateur.
- `reference_r/` : résultats de référence R **exhaustifs uniquement pour les méthodes dont toutes les lignes sont présentes**.
- `methods/` : méthodes et conditions de calcul, incluant le nombre de sujets indépendants et le contraste.
- `provenance/` : versions, résultats complets, manifestes SHA-256.
- `source/` : fichiers originaux et paramètres pour recalculer la même étude.
- `verification/verifier_dossier.py` : contrôle d'intégrité hors ligne.

Commande : `python3 verifier_dossier.py multiomics_dossier_recherche_v1.zip`.

Si le résultat R excède la limite du serveur JSON, le logiciel doit le signaler comme **incomplet dans le ZIP**. La conservation sur le disque du serveur R local est alors possible sur option explicite `PMX_MULTIOMICS_RESULTS_DIR`. Ne pas partager sans vérifier les données sensibles.

**Pour utiliser le backend R :** R 4.5.2 et Bioconductor 3.22 (DESeq2, edgeR, limma, statmod, jsonlite, plumber) ; lancer le serveur local avec `Rscript multiomics-engine/run_backend.R`. Dans les réglages MultiOmics, choisir **« R requis »** ; aucun jeu de données n'est transmis au serveur R sans cette sélection.

## 5. Interpréter avant de citer

- `log2FoldChange` dépend de l'ordre du contraste. Positif = expression accrue dans le groupe comparé.
- `pvalue` n'est pas la probabilité que l'hypothèse soit vraie.
- `padj` est une p-value ajustée pour tests multiples sur la famille de gènes définie par l'analyse. Une q-value faible ne certifie pas une identité génique ou une causalité.
- Les gènes filtrés **avant** le modèle, les valeurs `NA` du modèle, les lots, les exclusions et la version des packages doivent être rapportés.
- Ne pas confondre la répétition technique et une répétition biologique ; déclarer explicitement le **nombre de sujets indépendants**.

## 6. Exemple public vérifié : Airway / GEO GSE52778

Dans `airway` (Himes et al. 2014), huit prélèvements proviennent de **quatre lignées cellulaires** mesurées sous témoin `untrt` et traitement `trt`. Un test public indépendant de la démonstration utilisateur compare `run_backend_analysis` au DESeq2 natif, avec `~subject_id+condition`, et contrôle l'identité des p-values/q-values, l'exhaustivité du CSV et la direction de CRISPLD2 (ENSG00000103196). Voir [PR #86](https://github.com/rberrah/rberrah.github.io/pull/86) et [PR #84](https://github.com/rberrah/rberrah.github.io/pull/84).

La validation du même package DESeq2 avec deux modes d'appel mesure la fidélité de l'intégration, **pas la précision universelle de DESeq2**, ni une validation clinique indépendante.

---

### English quick guide

Use one `subject_id` for the two conditions belonging to the same patient or biological cell line. Give each sample unique `sample_id` and `assay_id`, and make sure assay columns match the sample sheet exactly. Select the paired design; load **raw integer RNA-seq counts**, choose the comparison versus reference group for confirmatory inference and explicitly opt into the **local R engine**. The reference model includes a **subject block** (`~subject_id+condition`). Four matched subjects are four independent biological units — not eight. After analysis, export the complete R method CSVs, source files, versions and SHA-256 manifest, and run the independent Python verifier. Invalid pairs, unestimable batches or missing contrasts must block fitting. This is scientific software quality assurance, not a clinical or universal FDR certification.
