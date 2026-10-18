# Résultats natifs DESeq2 : dispersions par gène, tendance, MAP et vérité

**18 octobre 2026.** Analyse de sensibilité de DESeq2 avec de vrais
ajustements R sur des données **simulées**, sans modification de
l'algorithme de production de MultiOmics.

## Méthode et traçabilité

[PR #55](https://github.com/rberrah/MultiOmics/pull/55), fusionnée
dans `main` au SHA `341cce4c9d3453760ccf4e12e941684ec076ebdc`.
[Workflow R natif](https://github.com/rberrah/MultiOmics/actions/runs/39161667198),
[artefact 384 études×bras](https://github.com/rberrah/MultiOmics/actions/runs/39161667198/artifacts/12148968286).

96 cohortes nulles NB : 32 par scénario avec les répétitions 4001–4032,
400 gènes et 40 sujets simulés chacune. Quatre vrais appels
`DESeq2::nbinomWaldTest` par étude, tous utilisant **les mêmes comptes,
la même normalisation, le même lot, le même préfiltrage et le même
contraste**, avec trois vecteurs de dispersion internes (Cox–Reid
par gène `dispGeneEst`, tendance `dispFit`, MAP `dispersions`)
ou l'alpha de génération exact. `cooksCutoff=FALSE` est appliqué
**à tous** les bras pour isoler l'effet des dispersions ;
ce contexte **n'est pas identique** à l'analyse DESeq2 par défaut.
Le filtrage indépendant reste activé et **384/384 recalculs BH**
indépendants sont concordants.

## Nombre d'études avec ≥1 faux appel BH, dénominateur de 32 par cellule

| Nul global | MAP final | Cox–Reid par gène | Tendance ajustée | Vraie dispersion |
| --- | ---: | ---: | ---: | ---: |
| Équilibré | 5 (15,63 %) | 5 (15,63 %) | 3 (9,38 %) | 3 (9,38 %) |
| Forte variation profondeur/dispersion | 2 (6,25 %) | 3 (9,38 %) | 2 (6,25 %) | 0 (0 %) |
| Groupes 12/28 et lot partiellement confondu | 5 (15,63 %) | 6 (18,75 %) | 2 (6,25 %) | 2 (6,25 %) |

La FDR empirique au nul global est la **fraction des études avec au
moins un faux appel** (FWER ici). Les intervalles exacts individuels
95 % sont larges : par exemple 0/32 admet une borne supérieure
d'environ 10,89 % ; 2/32 un intervalle d'environ 0,77–20,81 % ;
5/32 environ 5,28–32,79 %. Ils ne démontrent pas un FDR nominal
contrôlé à 5 %, encore moins universellement.

## Mesures des estimations de dispersion et du Wald

| Scénario | RMSE log α MAP | RMSE log α par gène | RMSE log α tendance | Moyenne Z² MAP | Moyenne Z² tendance | Moyenne Z² vérité |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Équilibré | 0,209 | 0,260 | **0,046** | 1,039 | 1,027 | 0,998 |
| Forte dispersion | 0,195 | 0,231 | **0,049** | 1,038 | 1,030 | 1,004 |
| Groupes déséquilibrés | 0,212 | 0,264 | **0,052** | 1,053 | 1,043 | 1,007 |

`RMSE log α` compare les vraies dispersions par gène aux valeurs
internes natives DESeq2 (sauf l'oracle où RMSE = 0 par construction).
Les résultats montrent descriptivement une tendance estimée plus
proche de la vérité de **ces générateurs** et des moments du Wald
légèrement plus proches de leur valeur asymptotique nominale.

**Risque de biais majeur du banc d'essai :** le générateur source
définit une moyenne croissante selon l'index du gène
(`exp(seq(log(45),log(350),...))`) et une dispersion lisse croissante
avec ce même index (`seq(.08,.20,...)` ou `seq(.30,.70,...)`).
Le modèle d'une **tendance moyenne–dispersion** est donc
structurellement favorisé : très faible erreur de la tendance
ne signifie pas qu'elle soit préférée au MAP pour des données
réelles ayant des variations individuelles de dispersion hors tendance.

## Conclusion scientifique

Le signal de surdétection des faux appels peut dépendre de
l'estimation de la dispersion et de son ajustement final, mais
**aucune méthode n'est sélectionnée a posteriori** dans l'outil.
Ces données ne démontrent ni défaut universel de DESeq2, ni
supériorité générale de la tendance, ni validité biologique et
encore moins clinique de nouvelles p-values.

**Contre-épreuve indispensable :** un générateur Gamma–Poisson
écrit indépendamment, avec dispersions aléatoirement hétérogènes
autour de la tendance et même évaluation native des quatre bras.
Voir [protocole de falsification indépendant](INDEPENDENT_GAMMA_POISSON_DISPERSION_20261018.md)
et [issue #44](https://github.com/rberrah/MultiOmics/issues/44).
