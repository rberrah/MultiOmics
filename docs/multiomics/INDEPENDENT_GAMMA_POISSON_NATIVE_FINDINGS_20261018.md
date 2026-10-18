# MultiOmics — résultats réels R : falsification de l'avantage de tendance DESeq2

**18 octobre 2026 — résultat de simulation informatique indépendante, avec
deux rejeux R natifs concordants.** Méthodes et horizon définis dans
[le protocole préspécifié](INDEPENDENT_GAMMA_POISSON_DISPERSION_20261018.md).

## Source, versions et preuves conservées

- Les vrais modèles R de [PR #56](https://github.com/rberrah/MultiOmics/pull/56)
  ont été exécutés à partir de la tête de source
  `92df0bf11cc5695ca1861265e66c0c994dba8d88`,
  antérieure au nettoyage ultérieur des workflows redondants.
- Première exécution :
  [GitHub Actions #39163552036](https://github.com/rberrah/MultiOmics/actions/runs/39163552036),
  [artefact A](https://github.com/rberrah/MultiOmics/actions/runs/39163552036/artifacts/12149703311).
- Deuxième exécution :
  [GitHub Actions #39163551827](https://github.com/rberrah/MultiOmics/actions/runs/39163551827),
  [artefact B](https://github.com/rberrah/MultiOmics/actions/runs/39163551827/artifacts/12149708557).
- Les deux archives contiennent chacune **384 lignes** et 12 cellules.
  Comparaison indépendante Python par
  `(scenario,replicate,seed,arm)` : **384 lignes appariées sans
  duplication**, **16 champs numériques/booléens communs comparés**,
  **0 désaccord**, et différences maximales des q-values BH vérifiées
  **égales à 0**. Cette concordance est un contrôle de rejeu ;
  les deux ajustements partagent DESeq2 et **ne sont pas deux validations
  biologiques ou interlaboratoires indépendantes**.
- Version : R **4.5.2**, Bioconductor **3.22**, 171 dépendances figées.
  Le générateur est **Gamma–Poisson latent puis Poisson**, distinct
  du code source `rnbinom()` original ; 400 gènes et 40 sujets
  **synthétiques**, effet réel de condition strictement nul.
- Les 96 nouvelles cohortes originales utilisent les répétitions
  **5001–5032** dans trois plans ; les quatre bras sont appariés
  sur **exactement** les mêmes matrices, facteurs de taille,
  modèles `~batch+condition` et familles BH. Dans **tous les bras**,
  `cooksCutoff=FALSE` : pas de comparaison directe en taux absolus
  avec le paramétrage utilisateur DESeq2 par défaut.

## Résultats primaires : fréquence des études avec ≥1 faux appel BH, q≤0,05

Au **nul global**, cette fréquence est égale à la FWER et à
l'estimateur du FDR. L'unité de réplication est **l'étude entière**,
pas les 400 gènes.

| Générateur Gamma–Poisson, 32 études | MAP final DESeq2 | Cox–Reid par gène | Tendance `dispFit` | Vraie dispersion α |
|---|---:|---:|---:|---:|
| Équilibré, dispersions hétérogènes | **2/32 (6,25 %)** | 2/32 (6,25 %) | **15/32 (46,88 %)** | 1/32 (3,13 %) |
| Variation de profondeur et forte dispersion | **3/32 (9,38 %)** | 4/32 (12,50 %) | **25/32 (78,13 %)** | 1/32 (3,13 %) |
| Groupes 12/28, lot partiellement confondu | **5/32 (15,63 %)** | 6/32 (18,75 %) | **25/32 (78,13 %)** | 3/32 (9,38 %) |

IC binomiaux exacts individuels à 95 % :
`2/32 [0,8;20,8] %`, `3/32 [2,0;25,0] %`,
`5/32 [5,3;32,8] %`, `15/32 [29,1;65,3] %`,
`25/32 [60,0;90,7] %`. **Ces IC ne sont pas simultanés** ;
ils démontrent une grande différence observée entre les
générateurs, mais ne certifient pas un FDR nominal à 5 % pour MAP
ou pour les autres bras.

Les différences au sein des mêmes cohortes confirment que ce
n'est pas un changement d'échantillons qui explique la comparaison
MAP–tendance : **57/96 études ont un faux appel exclusivement
avec la tendance**, **2/96 exclusivement avec MAP**,
et les autres discordances sont conservées dans le fichier par étude.
Il s'agit d'un résultat descriptif, pas d'un test hiérarchique
confirmatoire préenregistré.

## Diagnostics de dispersion et de test

| Source | RMSE moyenne `log(α estimée/α vraie)` : MAP | Cox–Reid par gène | Tendance `dispFit` |
|---|---:|---:|---:|
| Équilibré | **0,2375** | 0,2579 | **0,6944** |
| Forte dispersion | **0,2266** | 0,2436 | **0,6887** |
| Groupes déséquilibrés | **0,2370** | 0,2566 | **0,6822** |

L'estimation MAP a ici une erreur de dispersion considérablement
plus faible que la tendance seule : c'est **l'inverse** de l'étude
précédente où la vérité `α_g` suivait presque parfaitement
une fonction régulière de la moyenne d'expression.

Même avec de nombreux faux appels, la moyenne des statistiques
de Wald au carré reste parfois proche de 1
(tendance : **1,004** dans le scénario équilibré malgré
**15/32** études avec des faux appels). **Le deuxième moment
central ne suffit donc pas à caractériser les queues des p-values
et la probabilité des événements rares** ; examiner la queue
extrême des tests est indispensable.

## Interprétation et limites obligatoires

**Résultat robuste dans les simulateurs :** l'avantage apparent de
la tendance `dispFit` de DESeq2 est **dépendant de la relation
moyenne–dispersion imposée par le générateur**. En introduisant
un résidu de dispersion individuel `epsilon~N(0,0.65²)`,
l'estimateur MAP final peut être nettement meilleur que la
tendance et les taux de faux appels de cette dernière se dégradent
fortement. Ce résultat **falsifie toute extrapolation automatique**
du bon comportement de la tendance depuis la simulation
structurellement favorable.

**Non démontré** : supériorité universelle du MAP sur données RNA-seq
réelles, FDR contrôlé à 5 %, cause unique de l'inflation dans les
études originales, garantie pour chaque étude ou validation Omics→PK.
Les deux générateurs demeurent artificiels et NB, les effets de lot,
les tailles d'échantillon et les dispersions sont stylisés.

Les 384 ajustements ont passé le contrôle numérique BH ;
**aucune modification des p-values ou du moteur DESeq2 de
production n'est justifiée**. Une étude ultérieure doit examiner
les queues nulles, le traitement Cook, les dispersions sur données
biologiques publiques et des modèles de comptage non-NB
sans choisir de méthode d'après l'expérience la plus favorable.

L'[issue scientifique #44](https://github.com/rberrah/MultiOmics/issues/44)
reste **ouverte** ; la présente étude doit être citée dans le dossier
Research v1 comme preuve de **capacité à détecter des limites
inférentielles**, et non comme une certification universelle.
