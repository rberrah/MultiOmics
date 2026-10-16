# MultiOmics — étude des facteurs de taille et du risque de faux positifs RNA-seq

Date du protocole : 14 octobre 2026, avant exécution. Ce chantier ne modifie ni les résultats des patients, ni le moteur d'inférence, ni les seuils du site.

## Motivation

Les campagnes indépendantes #52, #75 et #77 rapportent une fréquence de faux appels sous nul global au-dessus de 5 % dans certains scénarios. Les contrôles précédents portent sur le filtrage, la dispersion et la référence normale/t du test de Wald, mais ne testent pas directement le rôle des facteurs de taille estimés par rapport aux vrais facteurs de librairie du générateur.

## Hypothèse et conditions fixées avant résultats

Générateur binomial négatif, 400 gènes tous nuls, expression entre 45 et 300, dispersion Uniforme(0,08 ; 0,32), facteurs de librairie log-normaux et effet de lot log(1,35). Deux plans distincts : 40 sujets (20/20) et 80 sujets (40/40), avec lots parfaitement équilibrés entre conditions et identifiés dans le modèle.

40 simulations indépendantes par plan, chacune analysée de deux façons en utilisant les mêmes matrices et étiquettes :

- DESeq2 standard : estimation interne des facteurs de taille et des dispersions, Wald, filtrage indépendant et correction BH.
- DESeq2 avec véritables facteurs multiplicatifs de librairie issus du générateur (normalisés par moyenne géométrique). Les dispersions restent estimées par DESeq2, puis Wald, filtrage et BH suivent la même procédure.

Même préfiltre, défini sans groupe : au moins 10 comptages dans 20 % des échantillons, minimum deux. Les résultats sont appariés par jeu simulé et jamais optimisés pour produire le taux le plus favorable.

## Critères et limites

Sous un nul global, une découverte quelconque est nécessairement fausse, donc la fréquence de jeux avec une découverte est la FWER empirique et l'estimateur de la FDR. Donner, pour chaque bras, les IC binomiaux exacts 95 %, les p-valeurs brutes < 0,05, les q-valeurs absentes, le nombre de gènes testés, la dispersion estimée/vraie et l'erreur de normalisation. La différence appariée est présentée avec bootstrap de 2 000 rééchantillonnages, purement descriptive.

La campagne de 40 répétitions par plan est diagnostique et manque de précision pour certifier un niveau 5 %. Les vrais facteurs de librairie ne sont pas disponibles chez les chercheurs, et les dispersions recalculées avec facteurs imposés empêchent d'isoler un effet causal pur de la normalisation.

Aucune substitution automatique, validation clinique ou certification de FDR. Si des anomalies persistent : accroître N de simulations, varier le générateur indépendamment, tester la présence de vrais effets, documenter le mécanisme de dispersion, étudier d'autres tailles d'échantillons et confronter les méthodes R natives à leurs hypothèses.

## Reproduction

Exécuter le workflow GitHub Actions multiomics-rnaseq-oracle-normalization.yml. Le script est multiomics-engine/validation/rnaseq_normalization_oracle_audit.R. Tous les fichiers per_dataset.csv, summary.csv, paired_differences.csv, manifest.json et errors.csv (s'il existe) doivent être conservés. Une CI verte ne prouve pas le contrôle du FDR.
