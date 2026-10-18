# Diagnostic indépendant : dispersion NB et statistique de Wald DESeq2

**Protocole fixé le 18 octobre 2026, avant lecture des nouveaux résultats R.**
Suivi de l'[issue scientifique #44](https://github.com/rberrah/MultiOmics/issues/44).
Ce travail teste un mécanisme possible de l'inflation du FDR **dans les
générateurs artificiels**, sans changer les analyses des utilisateurs.

## Situation motivant le test

La qualification RNA-seq (128 études nulles indépendantes par scénario) avait
montré des taux de faux appels DESeq2 **22/128, 23/128 et 18/128**, avec
des bornes simultanées à 95 % toutes supérieures au seuil BH nominal de 5 %.
Un test complémentaire de 96 études appariées n'a pas supprimé l'excès
observé lorsque les facteurs de profondeur véritables étaient fournis.

Le contrôle numérique BH des p-values natives a réussi, mais ne
prouve pas une distribution correcte des p-values avant ajustement.

## Trois bras sur EXACTEMENT le même jeu simulé

La fonction génératrice d'origine utilise 400 gènes, 40 échantillons,
des paramètres de dispersion négative binomiale **connus**, des effets
gène×lot réels et un coefficient de condition planté exactement nul.
Elle est relancée avec une nouvelle série de graines (3001–3032)
**distincte** des séries initiales 1–128, 1001–1128 et 2001–2032.

| Bras natif R/Bioconductor | Facteurs de taille | Dispersion NB par gène | Test |
|---|---|---|---|
| `native_default` | Estimés par DESeq2 | Estimée et modérée par DESeq2 | Wald + BH |
| `oracle_true_dispersion` | Estimés par DESeq2 | **Vraie α génératrice**, connue seulement en simulation | Wald + BH |
| `oracle_true_depth_and_dispersion` | **Vraie profondeur génératrice** | **Vraie α génératrice** | Wald + BH |

Les trois bras utilisent **les mêmes 40 sujets, la même matrice de
comptes entiers, la même table des lots, le même contraste
`condition treated vs control`, le même préfiltrage** de production
et le même type de test de Wald. L'API native de DESeq2 est utilisée :
`DESeq(dds)` pour le premier bras, puis `sizeFactors(dds)`,
`dispersions(dds)` et `nbinomWaldTest(dds)` pour les oracles.
Aucune p-value n'est calculée par une formule substitutive.
L'ajustement BH est néanmoins **rejoué indépendamment** à partir
des p-values natives de chaque bras.

Le véritable générateur fournit des paramètres par gène. Leur
alignement est contrôlé par **identifiant de gène**, et non par
position arbitraire dans la matrice.

## Taille et mesures, fixées avant lecture des résultats

- Trois mécanismes nuls déjà étudiés : équilibré, profondeur/dispersion
  élevées, groupes 12/28 et confondement de lot partiel.
- **32 études indépendantes × 3 scénarios × 3 bras = 288 vrais
  ajustements DESeq2**, répartis en neuf cellules, sans arrêt adaptatif.
- **Unité de Monte-Carlo : étude complète**, et non gène.
- Primaire : fréquence par cellule des études avec au moins un faux
  appel `padj≤0,05` ; au nul global, FDR = FWER.
- Incertitude : intervalle binomial exact individuel à 95 % par cellule,
  présenté comme descriptif ; aucun intervalle individuel n'est
  assimilé à une preuve simultanée de contrôle du taux nominal.
- Mécanisme : RMSE des `log(α_est/α_true)`, RMSE des vrais facteurs de
  taille, moyenne par étude des statistiques de Wald au carré
  (référence théorique descriptive `E[Z²]=1` pour `Z∼N(0,1)`),
  fréquence des p-values brutes ≤0,001 / 0,01 / 0,05 et nombres BH.
  Les gènes sont corrélés : les comptes de p-values entre gènes ne
  sont **pas** des réplications biologiques indépendantes.

**Interprétation conditionnelle :** si l'excès d'erreurs DESeq2 diminue
sensiblement lorsque la vraie dispersion NB est injectée, l'estimation
des dispersions devient une explication *possible* dans **ce simulateur**,
pas un défaut universel de DESeq2. Si l'excès persiste même avec
dispersion et profondeur vraies, il faut étudier les limites du
test de Wald ou du générateur, des valeurs aberrantes/Cook's,
du design et du filtre indépendant.

Avec seulement 32 études par cellule, les IC seront larges.
L'expérience ne pourra **pas** fournir une certification FDR à 5 %.

## Reproduction et préservation des résultats négatifs

- Source figée R 4.5.2, Bioconductor 3.22, DESeq2 natif ;
  contrôle de l'empreinte de l'environnement R de 171 packages.
- Scripts `rnaseq_deseq2_dispersion_wald_metrics.R`,
  `tests/test_rnaseq_deseq2_dispersion_wald_metrics.R` et
  `validation/rnaseq_deseq2_dispersion_wald_oracle.R`.
- Workflow `.github/workflows/deseq2-dispersion-wald-oracle.yml` :
  vrais 288 ajustements, archive par étude et CSV des neuf cellules,
  manifestes et versions, y compris en présence de résultats
  d'inflation défavorables.
- Tout jeu manquant, toute graine réutilisée, p/q native altérée ou
  association de dispersions par gène invalide l'exécution.
- Aucun changement de paramètres de production, aucune sélection
  automatique de méthode, aucune assertion de certification biologique
  ou clinique.

Cette expérience qualifie le **diagnostic logiciel et statistique
d'un mécanisme simulé**, pas une validation universelle des
associations RNA-seq, ni des futurs modèles Omics→PK.
