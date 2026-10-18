# Résultat quantifié : origine possible de l'inflation du FDR RNA-seq DESeq2

**18 octobre 2026 — preuve expérimentale informatique négative/positive conservée.**
Documentation du diagnostic **préspécifié avant le calcul** :
[DESEQ2_TRUE_DISPERSION_WALD_DIAGNOSTIC_20261018.md](DESEQ2_TRUE_DISPERSION_WALD_DIAGNOSTIC_20261018.md).

## Traçabilité

- [PR scientifique #53 fusionnée](https://github.com/rberrah/MultiOmics/pull/53)
  dans `main` au commit `8ce26764d4797e294995ac1787ec66cb43248509`,
  tête effectivement qualifiée `d96852f1c95bafcf381085dcf020c8272663b646`.
- [Job GitHub Actions de 288 ajustements R](https://github.com/rberrah/MultiOmics/actions/runs/39154798726),
  entièrement passé.
- [Artefact d'étude native](https://github.com/rberrah/MultiOmics/actions/runs/39154798726/artifacts/12146483687) :
  CSV intégral de **288 lignes**, CSV de neuf cellules,
  manifeste, session R et empreinte de la bibliothèque de 171 packages
  du conteneur Bioconductor **3.22, R 4.5.2**.
  L'artefact GitHub est temporaire ; résultats essentiels recopiés
  ci-dessous dans Git pour conservation durable.
- Original : données RNA-seq **synthétiques** négatives binomiales,
  400 gènes × 40 sujets simulés, trois plans exacts et 32 **nouveaux
  jeux indépendants par plan** (répétitions 3001–3032, exclusives
  des campagnes précédentes).

## Résultat : vrais ajustements DESeq2, trois modèles appariés

Les trois bras traitent **chaque matrice source identique**, avec même
lot, contraste nul de condition, filtrage de comptes, test de Wald et
seuil BH q≤0,05. Les paramètres `α` de dispersion de chaque gène et
les facteurs de profondeur de chaque sujet sont **connus du seul
générateur** et donc inutilisables directement comme méthode clinique
ou biologique par défaut.

| Nul global (32 études indépendantes) | DESeq2 normal : dispersion et profondeur estimées | Oracle : **dispersion α vraie**, profondeur estimée | Oracle : α et profondeur vraies |
| --- | ---: | ---: | ---: |
| Équilibré | **4/32 (12,50 %)** | **1/32 (3,13 %)** | **1/32 (3,13 %)** |
| Profondeur et dispersion hétérogènes | **5/32 (15,63 %)** | **1/32 (3,13 %)** | **2/32 (6,25 %)** |
| Groupes 12/28, lot partiellement confondu | **3/32 (9,38 %)** | **1/32 (3,13 %)** | **1/32 (3,13 %)** |

L'estimand est la fraction des **études complètes** avec au moins une
découverte fausse (`FDR = FWER` uniquement dans un nul global).
Les IC exacts binomiaux individuels à 95 % sont larges :
**1/32 : 0,079 %–16,22 %**, **2/32 : 0,77 %–20,81 %**,
**4/32 : 3,51 %–28,99 %**, **5/32 : 5,28 %–32,79 %** et
**3/32 : 1,98 %–25,02 %**. Ils ne constituent pas une démonstration
simultanée de calibration à 5 %.

Sur les **96 jeux appariés**, les résultats par paire montrent
**9 études avec faux appels seulement sous DESeq2 normal**,
**aucune étude avec faux appels seulement sous dispersion α vraie**
et **3 études avec faux appels dans les deux bras**.
La réduction descriptive est présente dans chacun des plans ;
elle est néanmoins mesurée sur un simulateur et reste sujette
à l'incertitude et au mécanisme exact de génération.

## Indices mécanistiques secondaires

| Scénario | Moyenne de `Wald Z²` normale | Moyenne de `Wald Z²` avec dispersion α vraie | Erreur RMSE des `log(α estimée/α vraie)` |
| --- | ---: | ---: | ---: |
| Équilibré | **1,0647** | **1,0115** | **0,2129** |
| Forte profondeur/dispersion | **1,0698** | **1,0311** | **0,1983** |
| Groupes déséquilibrés et lot | **1,0510** | **1,0095** | **0,2096** |

La statistique normale de référence sous un nul asymptotique aurait
`E[Z²] ≈ 1`. Ici, la moyenne des moments calculée par jeu puis
moyennée entre études se rapproche de 1 lorsque les dispersions
vraies sont utilisées. Cela **renforce l'hypothèse** que l'estimation
ou la régularisation des dispersions contribue à la déviation du test
dans ces mécanismes artificiels.

Ce constat ne distingue **pas encore** les estimations brutes,
leur shrinkage, les effets du modèle sur les SE, les valeurs aberrantes
de Cook's, le filtrage indépendant, la loi réelle des counts ou
les hypothèses du générateur. Une comparaison de grande taille
ou un générateur écrit indépendamment demeure nécessaire.

Le recalcul indépendant de **toutes les q-values BH natives** a réussi ;
la correction multiple n'est donc pas responsable d'une
erreur d'implémentation numérique dans ce jeu.

## Décision scientifique explicite

**Hypothèse soutenue** : les dispersions estimées peuvent contribuer
à l'inflation des faux positifs DESeq2 dans les simulateurs
négatifs binomiaux étudiés.

**Non établi** : une erreur universelle de DESeq2 sur de vraies données
RNA-seq, une supériorité générale edgeR/limma, une garantie de FDR 5 %
par l'oracle ou une correction directement utilisable en recherche
clinique. Il n'y a pas de calibration confirmatoire universelle
de MultiOmics par un workflow vert.

**Politique MultiOmics** : conserver les modèles, les p-values, les
q-values, les seuils et l'interface de production inchangés ;
maintenir une revue méthodologique des analyses confirmatoires ;
garder [l'issue #44](https://github.com/rberrah/MultiOmics/issues/44)
ouverte. Le prochain protocole devra dissocier la dispersion
gène-par-gène brute, sa modération et la statistique de Wald,
puis confronter ces résultats à un générateur indépendant.
