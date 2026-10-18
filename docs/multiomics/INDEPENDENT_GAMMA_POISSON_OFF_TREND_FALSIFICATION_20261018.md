# Falsification préspécifiée : source Gamma–Poisson indépendante et dispersion hors tendance

**Protocole figé le 18 octobre 2026, avant les calculs de cette branche.**

## Hypothèse à mettre à l'épreuve

La [campagne DESeq2 par composantes de dispersion](DESEQ2_DISPERSION_COMPONENT_ABLATION_20261018.md)
a montré une erreur logarithmique particulièrement faible de la tendance
dispersion–moyenne, sur un simulateur dont la vérité était elle-même
presque une fonction déterministe de la moyenne d'expression.
Cette construction peut artificiellement favoriser la tendance.

Le présent test ne reprend **ni les comptes ni la fonction génératrice
`rnbinom()` de l'étude précédente**. Une autre implémentation,
`independent_gp_null_study()`, construit une loi Gamma–Poisson :
`lambda(g,s) ~ Gamma(shape=1/alpha_g, scale=alpha_g*mu(g,s))`,
puis `Y(g,s) ~ Poisson(lambda(g,s))`.
Ainsi `E[Y]=mu` et `Var[Y]=mu+alpha_g*mu²` ; une expérience
unitaire distincte vérifie numériquement ses moments.

Pour chaque gène, `mu_g` est tiré uniformément en log entre 30 et 500,
mais **la dispersion α vraie présente en plus un résidu lognormal
indépendant** : `alpha_g = (0.04 + 9/mu_g) * exp(epsilon_g)`,
avec `epsilon_g ~ N(0,0.65²)`. Dans le scénario forte dispersion
la dispersion est multipliée par 2,5. Par construction, les α réelles
ne se trouvent plus toutes sur la tendance des moyennes. Les lots
possèdent un effet aléatoire propre aux gènes ; l'effet de condition
est **exactement nul**. C'est un scénario de falsification de la
performance relative de `dispFit`, **pas** une représentation
exhaustive de données humaines.

## Expérience appariée

Trois plans : groupes équilibrés et profondeurs modérées,
profondeurs hétérogènes et grandes dispersions, groupes 12/28 avec
lot partiellement confondu sans colinéarité.
Chaque étude comprend **400 gènes × 40 échantillons simulés**.
Trente-deux nouvelles études par scénario, répétitions **5001 à 5032**,
graines distinctes de toutes les campagnes antérieures ; au total
**96 études indépendantes** et **384 véritables ajustements R DESeq2**.

Les quatre modèles partagent exactement la **même** matrice de comptes,
les mêmes gènes après filtre, covariables, facteurs de taille
DESeq2, contraste `~batch+condition`, test Wald et BH :

| Bras de la même méthode native | Dispersion α utilisée |
|---|---|
| `map_final` | Estimation DESeq2 finale après MAP |
| `gene_wise_CR` | Estimateur gène-par-gène `dispGeneEst` |
| `fitted_trend` | Tendance ajustée `dispFit` |
| `true_nb_alpha` | Vérité du **nouveau générateur** (inconnue en données réelles) |

`cooksCutoff=FALSE` dans **tous** les bras pour isoler le rôle
des dispersions, à distinguer du paramétrage standard de production ;
filtrage indépendant activé de manière identique.
Les q-values natives sont recalculées par l'auditeur BH indépendant.

## Résultats mesurés et limites d'interprétation

Les sorties enregistrent pour chaque étude/bras la fréquence de faux
appels BH, la fraction de p-values brutes petites, les moments
de la statistique Wald, la RMSE logarithmique des dispersions
par rapport à la vérité et l'amplitude empirique du résidu hors tendance.
Les résumés sont calculés sur **les 32 études indépendantes** de
chaque cellule avec intervalles binomiaux exacts individuels à 95 %,
sans pseudoréplication des 400 gènes. Les effectifs sont trop réduits
pour certifier un FDR nominal à 5 % ou un effet différentiel de méthode.

**Résultat potentiellement falsifiant :** la tendance `dispFit`
n'est plus systématiquement plus proche de la dispersion réelle
ou ne réduit plus autant les faux appels par rapport à `map_final`.
Les résultats positifs comme négatifs seront conservés. Si les
relations persistent, une explication simpliste par le seul
alignement structurel au simulateur devient moins convaincante,
sans prouver une supériorité de la tendance sur données biologiques.

Cette étude n'est **pas un jeu externe de patients**, ni une
validation interlaboratoire, ni un test de biomarqueurs, ni un
benchmark universel ; deux générateurs artificiels partagent encore
des hypothèses NB et les données restent sans effet de condition.

## Reproduction et frontières des preuves

- Source `multiomics-engine/rnaseq_independent_gamma_poisson_null.R`.
- Tests de moments, de graines, de dispersion hors tendance et
  de rôle des échantillons dans
  `tests/test_independent_gamma_poisson_null.R`.
- Programme native R `validation/rnaseq_independent_gp_dispersion_components.R`.
- Réconciliation 12 cellules et 384 lignes, refus d'études
  recyclées ou manquantes, dans `rnaseq_gp_dispersion_summary.R`
  et `tests/test_independent_gp_dispersion_summary.R`.
- Workflow `.github/workflows/deseq2-independent-gamma-poisson.yml` :
  R 4.5.2 / Bioconductor 3.22 et 171 packages natifs épinglés,
  artefacts détaillés à conserver au-delà de la durée de GitHub Actions.
- Statut scientifique final attendu : **revue statistique indépendante
  requise**. Ne changer aucun estimateur de dispersion, test Wald,
  p-value, seuil BH, code de production ni étiquette de certification
  d'après les résultats d'une seule campagne.

L'[issue FDR #44](https://github.com/rberrah/MultiOmics/issues/44)
reste ouverte pour une validation indépendante et l'exploration de
Cook/filtrage, des plans biologiques et des autres générateurs.
