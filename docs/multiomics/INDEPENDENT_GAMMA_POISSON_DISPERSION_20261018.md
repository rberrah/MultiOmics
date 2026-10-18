# Préspécification : générateur Gamma–Poisson indépendant et dispersions hors tendance

**18 octobre 2026, figé avant toute exécution DESeq2 du nouveau générateur.**
Suivi de l'[issue #44](https://github.com/rberrah/MultiOmics/issues/44).
Niveau : **contre-épreuve de simulation**, non validation biologique externe.

## Pourquoi une nouvelle source indépendante ?

Les trois générateurs négatifs binomiaux précédents utilisaient
`rnbinom` avec des dispersions presque déterministes de la
moyenne d'expression par l'index du gène. Une tendance
moyenne–dispersion peut donc sembler artificiellement optimale.
[Résultats précédents à conserver](DESEQ2_DISPERSION_COMPONENT_NATIVE_FINDINGS_20261018.md).

Le nouveau générateur **ne lit pas et n'appelle pas la fonction
d'origine**. Il échantillonne une intensité latente
`lambda ~ Gamma(shape=1/alpha, scale=mu*alpha)`, puis
`Y ~ Poisson(lambda)`. Ainsi `E(Y)=mu` et
`Var(Y)=mu+alpha*mu²`, ce qui constitue une
**construction alternative de la loi NB**. Des tests statistiques
sur les moments de 40 000 tirages vérifient séparément cette identité.

### Source hors tendance et plans

- 400 gènes et 40 échantillons ; **aucun effet de condition**.
  Moyennes de gènes hétérogènes tirées en log-uniforme entre 30 et 500.
- Dispersion véritable :
  `alpha_g=(0.04+9/mu_g) × exp(epsilon_g)`,
  avec `epsilon_g ~ N(0,0.65²)`. Les résidus par gène sont
  non réductibles à une tendance déterministe.
- Effets de batch indépendants entre gènes :
  `beta_batch,g ~ N(0,0.32²)`. Le lot conditionne la moyenne,
  mais l'effet de condition demeure nul.
- Trois plans fixés :
  `gp_balanced_heterogeneous` (20/20, profondeur .85–1.15),
  `gp_depth_heterogeneous` (20/20, profondeur log-uniforme .25–3.5
  et alphas ×2,5), `gp_unbalanced_batch` (12/28, lot partiellement
  confondu mais modèle identifiable).
- Reprises **5001–5032 par scénario**, graines spécifiques :
  32 études complètes par plan, soit 96 jeux indépendants.
  Chaque jeu est analysé par quatre bras, **384 vrais
  ajustements R/DESeq2**.

## Analyses appariées et métriques

Le moteur R utilise le même code de `estimateDispersions` et
`nbinomWaldTest`, le même filtre de comptes et les mêmes facteurs
de taille pour chaque étude, comme dans la PR #55 :
dispersions finales MAP, Cox–Reid `dispGeneEst`, tendance
`dispFit`, et vraie alpha de ce **nouveau** générateur.
Le filtre de Cook reste désactivé **dans les quatre bras**
pour ne pas mélanger les effets ; le filtrage indépendant et
BH restent natifs puis audités numériquement.

Estimand primaire : taux des **études entières** ayant au moins
une découverte fausse BH à q≤0,05 (FDR = FWER seulement au
nul global), avec IC binomiaux individuels exacts à 95 %.
Métriques secondaires : erreur log dispersion vs vérité,
moment moyen `Wald Z²`, fréquences des p-values brutes,
comptes de gènes filtrés et audit indépendant des q-values.
Les gènes ne constituent jamais des observations
Monte-Carlo indépendantes.

Le protocole impose douze cellules complètes, appariement
exact des 32 graines et mêmes identifiants de gènes par bras.
Toute erreur de source, p-value, API native, tendance
non estimable ou variance non finie invalide le run.
Aucun arrêt anticipé adaptatif, aucun ajustement
a posteriori des paramètres en fonction du FDR observé.

## Portée des conclusions

Si la tendance ne reste plus la meilleure estimation de
`alpha` quand des résidus biologiques artificiels sont
ajoutés, cela corroborera la **dépendance de la conclusion
antérieure au générateur**, sans démontrer que le MAP soit
universellement meilleur. Si le FDR simulé reste inflationniste,
il faudra étudier les statistiques de Wald, filtres Cook,
régularisation et une validation RNA-seq réelle.

Un générateur Gamma–Poisson indépendant est une **nouvelle
implémentation et famille de paramètres simulés**, pas
une source empirique externe, et les ajustements reposent
toujours sur le même package DESeq2. Avec 32 études/cellule,
les IC restent trop larges pour certifier FDR=5 %.

**Aucun changement de production, de seuil p/q,
d'interface, de licence ou de politique confirmatoire
ne découlera automatiquement de cette expérience.**
