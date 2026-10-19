# Sacurine v3 — comparaison indépendante xcms ↔ OpenMS

**Protocole préspécifié le 19 octobre 2026 avant l'exécution OpenMS.**

Cette qualification répond à la partie « second preprocessing implementation » de
l'issue #47. Elle réutilise exactement les **12 mzML Sacurine v3 authentifiés**
(DOI 10.5281/zenodo.10130758), leurs MD5 éditeur et leurs vrais rôles
6 urines / 3 pooled-QC / 3 blancs. Aucun résultat biologique, âge, sexe ou BMI
n'est utilisé pour choisir les paramètres.

## Deux implémentations réellement distinctes

- référence déjà qualifiée : R/Bioconductor `xcms`, paramètres Research v1
  gelés ;
- comparateur : **OpenMS 3.6.0**, exécutables officiels dans le conteneur
  `ghcr.io/openms/openms-executables:3.6.0`. Le digest réellement tiré est
  archivé avec l'évidence et sera ensuite épinglé pour les rejeux.

OpenMS `FeatureFinderMetabo` attend des données LC-MS centroidées. Le mode de
pic des 12 fichiers est donc lu directement depuis les CV mzML **avant**
feature finding. Si les MS1 sont profile, `PeakPickerHiRes` est appliqué ;
s'ils sont centroidés, aucune recentroïdisation n'est effectuée ; tout mode
mixte/ambigu bloque le workflow.

## Pipeline OpenMS fixé

1. `FeatureFinderMetabo` par injection avec les **paramètres stables par
   défaut** de la release 3.6.0. Aucun paramètre n'est ajusté à partir du
   nombre de features xcms, du RSD, des blancs ou du chevauchement.
2. Alignement des 12 feature maps par `MapAlignerPoseClustering`, avec comme
   référence le **premier pooled-QC dans l'ordre d'injection publié**. Ce choix
   est source-only et ne dépend d'aucun résultat.
3. Groupement inter-runs par `FeatureLinkerUnlabeledQT` avec les paramètres
   stables par défaut OpenMS 3.6.0.
4. Export complet par `TextExporter`. Les features observées uniquement dans
   les blancs ne définissent pas l'univers comparé, mais les intensités des
   blancs restent conservées pour les features définies par biologique/QC.
5. RSD pooled-QC : `sd/mean` uniquement lorsque les 3 QC sont finis et la
   moyenne >0. Ratio blanc/urine : médiane blanc / médiane urine lorsque >=2
   valeurs finies existent dans chaque catégorie et médiane urine >0.
   Les autres cellules restent `NA`.

## Comparaison préspécifiée

Les identifiants de features ne sont **jamais** comparés entre logiciels.
Une correspondance exige un voisin **unique et réciproque** à <=20 ppm et
<=30 secondes, règle déjà utilisée dans la sensibilité xcms Research v1.

Sont rapportés sans seuil de succès a posteriori :

- nombre total de features xcms et OpenMS ;
- matches réciproques uniques ;
- aucune candidate / plusieurs candidates / non-réciproque ;
- fraction matchée sur **les deux dénominateurs** ;
- différences signées m/z (ppm) et RT ;
- RSD et ratio blanc/urine OpenMS ;
- concordance des écrans RSD >30 % et blanc/urine >=20 % uniquement pour les
  paires où les deux métriques sont calculables.

## Ce qu'un PASS signifie

Le PASS du workflow signifie uniquement : **source authentique, deux
implémentations exécutées, provenance complète et comparaison numérique
préspécifiée produite**. Un faible chevauchement ou un mauvais QC ne fait pas
échouer artificiellement la qualification : il est conservé comme résultat
scientifique négatif.

Même un chevauchement élevé ne valide ni identité chimique, adducts/isomères,
carryover, correction de dérive, précision interlaboratoire, biomarqueur,
FDR, association biologique ni relation Omics→PK.
