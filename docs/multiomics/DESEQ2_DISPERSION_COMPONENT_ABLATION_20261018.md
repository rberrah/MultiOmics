# Étude préspécifiée : décomposition des dispersions DESeq2 et test de Wald

**18 octobre 2026, protocole écrit avant lecture des résultats.**
Suivi de l'[issue #44](https://github.com/rberrah/MultiOmics/issues/44).

## Question scientifique

Sous trois générateurs synthétiques NB à coefficient de condition
strictement nul, MultiOmics a observé une inflation du FDR DESeq2.
L'introduction de la vraie dispersion `α` génératrice améliore nettement
la fréquence des faux appels dans une étude exploratoire antérieure
(32 études/scénario), sans démontrer le mécanisme exact. Le fabricant
de la méthode décrit trois composants distincts : estimation
**par gène**, tendance dispersion-moyenne et estimateur **MAP final**
régularisé. [DESeq2 : documentation des composantes
natives](https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

Le présent protocole examine lequel de ces trois composants diffère
le plus de la vérité *synthétique*, sans remplacer ni corriger DESeq2.

## Quatre bras sur les mêmes matrices de comptes

Tous les bras utilisent l'API native `DESeq2::nbinomWaldTest`,
avec les mêmes 40 sujets simulés, 400 gènes d'origine, table de lot,
contraste traité/témoin nul et préfiltrage `counts≥10` chez au moins
20 % des sujets. La normalisation DESeq2 par facteur de taille est
**estimée une seule fois par cohorte** et réutilisée à l'identique.

| Bras | α utilisé dans le même test Wald |
|---|---|
| `map_final` | Dispersions finales DESeq2 après estimation et régularisation MAP |
| `gene_wise_CR` | Estimations spécifiques aux gènes, ajustées Cox–Reid, avant l'étape de modération |
| `fitted_trend` | Valeurs `dispFit` de la tendance moyenne–dispersion ajustée |
| `true_nb_alpha` | Valeur de génération NB exacte de chaque gène (disponible uniquement en simulation) |

Les estimations et la tendance proviennent du **même objet**
`DESeq2::estimateDispersions`, et les quatre appels
`nbinomWaldTest` utilisent les mêmes paramètres et le même
`betaPrior=FALSE`. Cela est **une ablation statistique de nuisance**,
pas quatre logiciels ou validations indépendantes.

Afin de ne pas attribuer aux dispersions un effet du traitement
automatique des valeurs aberrantes, le filtre `cooksCutoff=FALSE`
est appliqué **aux quatre bras** lors de `results`. Le filtrage
indépendant des gènes reste activé et BH est recalculé
indépendamment depuis chaque famille de p-values natives
ayant une q-value finie. **Ce contexte diffère du comportement
DESeq2 habituel avec filtrage Cook**, et les résultats ne doivent
pas être confondus avec les taux de FDR des campagnes précédentes.

## Horizon, estimands et règles de décision fixés

- Trois plans source : équilibré ; forte variation de profondeur/
  dispersion ; groupes 12/28 avec lot partiellement confondu.
- Nouvelles répétitions **4001–4032**, exactement
  32 études indépendantes par plan ; non-réutilisation des graines
  1–128, 1001–1128, 2001–2032, 3001–3032.
- Quatre bras × 32 études × trois scénarios = **384 vrais ajustements
  DESeq2 Wald natifs**. Chaque étude produit quatre lignes appariées.
- Estimand de sécurité : fréquence des **études** avec au moins
  une fausse découverte BH à `q≤0,05`, égale à la FWER **uniquement
  au nul global**. IC Clopper–Pearson à 95 % individuels,
  larges avec 32 études, **non interprétables comme des garanties
  simultanées**.
- Diagnostics mécanistiques : erreur quadratique moyenne logarithmique
  `α_est/α_true` par étude ; statistique moyenne `Wald Z²` ;
  fractions de p-values brutes ≤0,05 et ≤0,001 ; nombres BH et
  niveau de filtrage indépendant. **Les gènes ne sont pas
  des réplicats biologiques indépendants**.
- Les douze cellules et toutes les études sont exigées ; aucun
  arrêt adaptatif et aucun choix de paramètres en fonction des
  résultats. Toutes les divergences, features écartées et
  estimations non finies provoquent un échec contrôlé.

## Portée et limites des conclusions

Un changement des faux appels entre `gene_wise_CR` et `map_final`
serait **compatible avec un rôle de la régularisation** dans ces
générateurs. Une meilleure correspondance au simulateur ne signifie
pas qu'un estimateur soit optimal sur une cohorte biologique réelle.
Le bras `true_nb_alpha` est un diagnostic inaccessible en pratique.

Cette étude n'isole pas une causalité universelle : la loi NB du
générateur, le design, les valeurs aberrantes/Cook, l'approximation Wald
et le filtrage indépendant restent des dimensions à examiner.
Une validation ultérieure **par un générateur écrit indépendamment**
demeure nécessaire.

Les scripts, graines, versions, résultats par étude et manifestations
d'erreur restent archivés. R 4.5.2, Bioconductor 3.22 et la fermeture
de **171 packages** sont fixés. Aucune modification des paramètres de
production, des p-values, des seuils de significativité ou de
l'interface MultiOmics n'est introduite ; statut final
**« revue statistique indépendante requise »**.
