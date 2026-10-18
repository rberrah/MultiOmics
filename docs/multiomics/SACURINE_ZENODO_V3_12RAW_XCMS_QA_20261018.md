# Qualification externe instrumentale — 12 acquisitions Sacurine v3 réellement mesurées

**Statut de cette phase : protocole et code créés ; le résultat instrumental final
dépend de l'exécution du workflow sur les 12 chromatogrammes originaux.**
Une CI de code PASS et une origine MD5 vérifiée ne suffisent pas à déclarer
une qualification instrumentale réussie.

## Données originales et indépendance

Source unique, version fixe : [Sacurine toy dataset Zenodo v3 — DOI
10.5281/zenodo.10130758](https://zenodo.org/records/10130758).
Le jeu contient six spectres mzML biologiques d'urines humaines,
trois véritables pools QC et trois véritables blancs analytiques,
plus `sampleMetadata_12samp_completed.tsv`. La MD5 publiée du TSV est
`f88044818627dee797c75d5a310f5c8e`. Les douze autres MD5 publiées
sont figées dans `scripts/qualify_sacurine_toy_zenodo.py` et revérifiées
par `scripts/verify_sacurine_zenodo_12raw.py`.

La phase antérieure [source seule](SACURINE_ZENODO_TOY_EXTERNAL_SOURCE_20261018.md)
a effectivement téléchargé et authentifié le TSV. **Cette phase B**
doit télécharger les vrais `.mzML` (~257 Mo), avant toute extraction,
puis traiter **les douze injections de la même acquisition**, sans fusionner
des sources venant des autres jeux Sacurine Zenodo.

La source donne des ordres d'injection absolus **27 à 185** : ne pas les
remplacer par les entiers 1–12. Il s'agit d'injections dispersées dans
une série beaucoup plus vaste, donc aucune correction fiable de la
dérive instrumentale ou du carryover ne peut être revendiquée avec
seulement trois QC et trois blancs. Seul le lot source `ne1` est admis.

## Méthode scientifique figée avant les résultats

- Acquisition, 12/12 : accès HTTPS Zenodo version figée ou originaux
  locaux institutionnels, empreintes MD5 et SHA-256 par fichier,
  taille plafonnée, réutilisation de cache non authentifié interdite.
  **Aucun chromatogramme originel** n'est ajouté au dépôt ni à l'artefact.
- CentWave natif `xcms`, `ppm=25`, `peakwidth=10–60 s`,
  `snthresh=8`, `prefilter=c(3,100)`, correction du temps de
  rétention Obiwarp (`binSize=0.6`), correspondance PeakDensity
  (`bw=8`, `minFraction=0.5`, `minSamples=1`,
  `binSize=0.025`, `ppm=15`), remplissage
  `ChromPeakAreaParam` et intégration des pics `into`.
  Ce sont **des paramètres de test logiciel pré-définis**, pas
  une optimisation d'instrument.
- Les **blancs ne définissent pas les groupes de features** lors
  de la correspondance, mais leurs signaux originaux restent présents
  dans la matrice quantitative exportée.
- Référence R/Bioconductor 4.5.2 / 3.22, fingerprint strict de
  l'environnement natif de 206 packages. Vérification de
  l'intégrité numérique entre matrice native et CSV exporté,
  définitions de features complètes conservées en RDS.
- Contrôle QC indépendant : pour chaque feature comportant trois
  valeurs QC finies et une moyenne positive,
  `RSD = sd(QC)/mean(QC)`.
- Contrôle des blancs indépendant : médiane des intensités
  de trois blancs / médiane des intensités de six urines,
  seulement si au moins deux valeurs finies existent dans
  chaque catégorie et si la médiane biologique est strictement
  positive. Les autres cas restent explicitement non évaluables.
- Les indicateurs `RSD > 30 %` et ratio blanc/urine `≥20 %`
  sont **uniquement des alertes descriptives**.
  Ils ne sont ni des tests de contamination identitaire,
  ni des p-values, ni un ajustement FDR.
- Aucune analyse par groupe (sexe, âge ou BMI de la source)
  n'est effectuée et aucune affectation de groupe n'est inventée.
  La matrice n'est pas jumelée à des données de concentration PK.

## Critères de réussite et échecs prévus

1. Source 12/12 authentifiée, type des échantillons 6/3/3,
   lot unique `ne1`, chronologie source non modifiée.
2. xcms exécuté sur **les vrais spectres** ; features non vides,
   identifiants et définitions uniques et correspondance numérique
   intégrale de l'export au résultat R.
3. Trois blancs exclus de la définition initiale des groupes
   chromatographiques mais conservés dans les valeurs quantitatives.
4. Recalcul manuel indépendant des RSD et ratios
   pour **chaque feature**, concordant aux méthodes QA du backend,
   et présence explicitement mesurée des colonnes `NA`
   lorsque les ratios ne sont pas évaluables.
5. Les écrans de dérive et comparaison de batch restent bloqués
   pour seulement trois QC / un seul lot, bien que les ordres
   d'injection soient sourcés.
6. Les paramètres, versions, SHA de code, MD5 et SHA-256 d'origine
   sont conservés dans des fichiers légers ; **aucun fichier mzML
   ou attribut clinique individuel** n'est publié comme artefact.
7. Si une source est indisponible, corrompue, ou un résultat R
   devient incompatible avec le backend, le job échoue et ne
   publie aucun faux succès instrumental.

## Exécution et limites

CI autonome :
`.github/workflows/sacurine-zenodo-12raw-xcms.yml`.
Les sorties quantitatives compactes sont archivées pendant 90 jours,
avec des avertissements permanents et sans raw spectra originales.

Même après une exécution positive, le niveau atteint restera
**qualification logicielle et instrumentale descriptive indépendante**,
pas validation de l'identification de molécules, de la précision
interlaboratoire, d'une association biologique, du FDR, de
biomarqueurs utilisables en clinique ni des liens sujets Omics→PK.

Pour un article utilisant MultiOmics, cet audit permet de documenter
un cas LC-MS externe avec blancs et QC mesurés mais ne remplace pas
le contrôle de qualité du projet de recherche individuel.
