# MultiOmics — preuves indépendantes de sensibilité xcms et diagnostic FDR DESeq2

**Études exécutées le 18 octobre 2026 ; résultats scientifiques négatifs conservés.**
Ce document contient les chiffres vérifiés, pas seulement la réussite des workflows.
**Il ne certifie ni les biomarqueurs, ni un FDR contrôlé à 5 % dans toute étude.**

## 1. Sensibilité xcms sur 12 spectres externes originaux

**Sources :** Zenodo [Sacurine toy v3, DOI 10.5281/zenodo.10130758](https://zenodo.org/records/10130758).
Six véritables acquisitions d'urine, trois injections pooled-QC et trois blancs ; les
empreintes MD5 originales des douze mzML, la feuille de métadonnées et les
positions originales d'injection 27–185 ont été vérifiées avant calcul.
Méthodes : `xcms 4.8.0` sous R 4.5.2 / Bioconductor 3.22, avec fermeture de
206 packages natives vérifiée. Aucune source d'intensités simulées n'a été utilisée
pour cette qualification instrumentale.

La [PR #50](https://github.com/rberrah/MultiOmics/pull/50),
fusionnée au commit `a1547508acd52243f45103a2ef0d147f94b06519`,
a comparé les **mêmes spectres** et les **mêmes** paramètres
d'alignement et de regroupement dans deux traitements xcms :

| Paramètre | Primaire | Variante préspécifiée |
| --- | ---: | ---: |
| `centWave ppm` | 25 | 20 |
| `centWave snthresh` | 8 | 10 |
| `peakwidth` | 10–60 s | 10–60 s |
| Alignement, regroupement et intégration | Identiques | Identiques |
| Features groupées | **3 135** | **2 982** |

Le calcul exige des appariements **uniques et réciproques** de features
à partir des médianes natives `mzmed` et `rtmed`, dans les fenêtres
fixées de **20 ppm et 30 s**. Les identifiants xcms d'un traitement
ne sont jamais supposés équivalents à ceux d'un autre.

| Contrôle indépendant | Valeur |
| --- | ---: |
| Features réciproquement uniques | **2 727** |
| Part des 3 135 features primaires appariées | **86,99 %** |
| Part des 2 982 features variantes appariées | **91,45 %** |
| Features primaires sans candidate | **190** |
| Features primaires avec plusieurs candidates ambiguës | **154** |
| Features primaires sans réciprocité univoque | **64** |
| Paires avec RSD pooled-QC évaluables dans les deux traitements | **2 497** |
| Concordance du dépistage RSD QC >30 % sur ces paires | **97,36 %** |
| Paires avec ratios blanc/urine évaluables dans les deux traitements | **368** |
| Concordance du dépistage ratio blanc/urine >20 % sur ces paires | **98,10 %** |

**Source vérifiable des données numériques :**
[GitHub Actions #39150564276](https://github.com/rberrah/MultiOmics/actions/runs/39150564276),
[table complète de 3 135 appariements, JSON et provenance](https://github.com/rberrah/MultiOmics/actions/runs/39150564276/artifacts/12144919816).
L'artefact GitHub est soumis à expiration ; les chiffres et méthodes
de cette page sont donc conservés durablement sous contrôle de version.

**Interprétation permise :** stabilité substantielle des associations m/z–RT,
au sens précis d'un appariement univoque, sous cette **perturbation combinée**
des paramètres de détection et concordance élevée des écrans analytiques
parmi les features **évaluables dans les deux conditions**.

**Interprétation interdite :** équivalence de molécules identifiées,
stabilité interlaboratoire, puissance biologique, sensibilité clinique ou
validité d'une sélection automatique des 2 727 paires.
Les seuls paramètres `ppm` et `snthresh` ont changé ensemble : on
ne peut pas attribuer l'effet observé à l'un ou l'autre.

Surtout, les résultats antérieurs sur le traitement primaire
([étude xcms 12-raw](SACURINE_ZENODO_V3_REAL_INSTRUMENT_FINDINGS_20261018.md))
restent défavorables sur le plan analytique : **1 843/2 850** RSD QC
évaluables >30 %, **389/535** ratios de fond mesurables ≥20 %, et
**2 600/3 135** ratios de blanc non évaluables. Une haute concordance de
drapeaux de qualité **ne fait pas disparaître ces drapeaux défavorables**.
Trois QC espacés dans une séquence incomplète ne suffisent pas à valider
une correction de dérive instrumentale.

Le traitement utilisé pour ce test de sensibilité est **encore xcms** :
ce n'est pas une comparaison indépendante avec OpenMS, MZmine ou une
référence annotée. La [qualification instrumentale #47](https://github.com/rberrah/MultiOmics/issues/47)
reste ouverte.

---

## 2. Diagnostic mécanistique du FDR DESeq2 par facteurs de profondeur connus

La campagne préalable avec 128 études nulles/scénario avait montré
une inflation DESeq2 dans les trois générateurs RNA-seq :
**22/128, 23/128 et 18/128** études avec au moins une découverte
BH à `q≤0,05`, les bornes inférieures des intervalles simultanés à
95 % dépassant 5 %. Ce résultat est documenté
[dans son analyse initiale](RNASEQ_HIGH_REP_FDR_FINDINGS_20261018.md)
et ne signifie pas que DESeq2 est universellement incorrect.

La [PR #51](https://github.com/rberrah/MultiOmics/pull/51),
fusionnée au commit `bcda814102b5cae8ae6e9c550295cd5a9b691cec`,
a étudié **une seule hypothèse de mécanisme** : la contribution des
facteurs de taille inférés sur les données. Les 96 nouvelles cohortes
nulles synthétiques utilisent des graines 2001–2032, indépendantes
des précédentes. Chacune est ajustée deux fois par **le véritable DESeq2** :
l'algorithme habituel, puis les **facteurs de profondeur exacts connus
du générateur**, inaccessibles sur des données biologiques réelles.

Le design `~batch+condition`, les **mêmes 40 sujets simulés** et les
400 gènes d'origine, le préfiltrage, les méthodes d'estimation de la
dispersion, le test de Wald, le filtrage indépendant et BH sont
**strictement identiques** dans les deux bras.

| Générateur totalement nul | DESeq2 habituel | Facteurs de profondeur oracle |
| --- | ---: | ---: |
| Équilibré | **5/32 (15,63 %)** | **4/32 (12,50 %)** |
| Forte variation de profondeur et dispersion | **4/32 (12,50 %)** | **5/32 (15,63 %)** |
| Groupes 12/28, lot partiellement confondu | **5/32 (15,63 %)** | **6/32 (18,75 %)** |

L'erreur quadratique moyenne des `log(size factors)` inférés face à
la vérité du générateur était, en moyenne : **0,03693**, **0,10601**
et **0,03730** selon le scénario ; elle était vérifiée **nulle**
dans les trois bras oracle. L'appariement étude-par-étude montre
respectivement 4, 4 et 4 cohortes contenant des faux appels dans
**les deux bras**, plus de rares cas propres à un seul bras.

**Preuves natives :**
[GitHub Actions #39151110619](https://github.com/rberrah/MultiOmics/actions/runs/39151110619),
[96 études appariées, table R réelle, IC binomiaux et environnement](https://github.com/rberrah/MultiOmics/actions/runs/39151110619/artifacts/12145178446).
L'intégrité des p-values/q-values a été vérifiée numériquement avec
une implémentation BH indépendante.

**Conclusion mesurée :** la connaissance exacte de la profondeur ne
supprime pas l'excès apparent de fausses découvertes dans les trois
simulations. L'estimation des facteurs de taille **n'explique donc pas,
à elle seule**, les résultats négatifs préalables.

**Limites :** seulement 32 études indépendantes par cellule, donc
incertitude Monte-Carlo élevée ; ce résultat n'identifie pas la cause
dominante, ni une anomalie causale particulière de DESeq2 sur des
données réelles. Les propriétés du générateur, les estimations
de dispersion, la statistique de Wald, le filtrage indépendant et les
paramètres d'ajustement doivent faire l'objet de tests mécanistiques
supplémentaires, idéalement avec un générateur indépendant.
L'[issue #44](https://github.com/rberrah/MultiOmics/issues/44)
reste ouverte.

---

## Décision pour Research v1

Ces études ajoutent **des preuves scientifiques contrôlables**, y
compris les données négatives : authentification instrumentale,
stabilité descriptive de coordonnées xcms, contrôles QC/blancs,
appariement des cohortes nulles et vérification de calculs R natifs.

**Elles ne constituent pas un certificat de FDR universel, d'identité
moléculaire, de précision multi-instrument ni de valeur prédictive
Omics→PK.** Les analyses confirmatoires doivent être adaptées au
plan de chaque étude et revues par un spécialiste ; les données
exploratoires calculées dans le navigateur ne doivent pas être
présentées comme des résultats R de référence.

Les moteurs de production et seuils des utilisateurs **n'ont pas été
modifiés** en fonction de ces benchmarks. Cette décision préserve la
traçabilité et évite une correction post hoc uniquement destinée à
faire passer un benchmark à 5 %.
