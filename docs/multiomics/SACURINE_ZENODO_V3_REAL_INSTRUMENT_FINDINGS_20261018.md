# MultiOmics — résultats instrumentaux externes réels Sacurine v3 (18 octobre 2026)

**Résultat observable de l'exécution scientifique, non une revendication de qualité biologique globale.**

## Provenance et méthode

- Source réelle indépendante : Zenodo Sacurine v3 [DOI 10.5281/zenodo.10130758](https://zenodo.org/records/10130758), **6 acquisitions d'urine, 3 injections pooled-QC, 3 blancs** du même sous-ensemble. Les 12 MD5 publiées, leurs SHA-256 calculées et les **ordres d'injection originaux** 27 à 185 sont consignés.
- Code validé de la [PR #48](https://github.com/rberrah/MultiOmics/pull/48), commit final `6c25058bcd637dd8e0d43d4bc3deaaaee63874c5`, squash `e34d82f541f5249a8d5023888c2b1dde048c31e6`. Exécution [GitHub Actions #39146123154](https://github.com/rberrah/MultiOmics/actions/runs/39146123154), [artefact original](https://github.com/rberrah/MultiOmics/actions/runs/39146123154/artifacts/12144095556) conservé temporairement par GitHub pendant 90 jours.
- Bioconductor figé (R 4.5.2, Bioconductor 3.22, xcms 4.8.0). Les paramètres `centWave / obiwarp / PeakDensity`, les définitions de feature `xcms`, l'intégration `into` et les critères descriptifs [sont définis a priori](SACURINE_ZENODO_V3_12RAW_XCMS_QA_20261018.md). **La reproductibilité numérique du CSV face à la matrice R d'origine a été vérifiée.**
- Un deuxième run GitHub avant le dernier commit [#39145982784](https://github.com/rberrah/MultiOmics/actions/runs/39145982784) est également vert et livre les mêmes valeurs ; des tests ultérieurs ne remplacent pas le numéro de run/artefact effectivement examiné.

## Résultats natifs observés

| Mesure | Observé |
|---|---:|
| Spectres authentifiés et traités | **12 / 12** |
| Pics chromatographiques natifs `xcms` | **41 588** |
| Features moléculaires *non identifiées* | **3 135** |
| Fractions de la matrice xcms sans aire finie | **26,86 %** |
| Features avec 3 valeurs QC et RSD calculable | **2 850 / 3 135** |
| Features avec RSD QC **> 30 %** | **1 843 / 2 850 = 64,67 %** des évaluables |
| Features avec ratio médian blanc/urine calculable | **535 / 3 135** |
| Features avec ratio blanc/urine **≥ 20 %** | **389 / 535 = 72,71 %** des évaluables |
| Features dont le ratio de fond des blancs n'est **pas évaluable** | **2 600 / 3 135 = 82,93 %** |
| Features sur lesquelles **les deux contrôles sont évaluables** | **426 / 3 135** |
| Features satisfaisant simultanément les *deux écrans descriptifs* parmi ces 426 | **21** |

**Le chiffre « 21 » n'est pas une sélection analytique validée.** Ces 21 features sont seulement celles présentant un RSD QC ≤30 % et un ratio blanc/urine <20 % **dans le sous-ensemble où les deux ratios peuvent être calculés**. Les features aux données de blanc non évaluables ne peuvent pas être qualifiées de « fond propre » : `NA` **n'est pas** zéro.

L'observation `PASS` signifie : téléchargement des vrais fichiers, intégrité, correspondance échantillons/metadata, exécution native et concordance numérique. Les valeurs de QC mesurées **signalent un problème potentiel de stabilité et de fond analytique**. Elles ne démontrent pas que les 1 843 features RSD >30 % sont toutes biologiquement fausses ou que les 389 ratio ≥20 % proviennent d'une contamination moléculaire confirmée. Le résultat dépend notamment des paramètres de détection choisis sans optimisation a posteriori.

## Requalification indépendante de la sortie

En plus de la relecture native sous R, `scripts/audit_sacurine_12raw_evidence.py` relit chaque feature avec la **bibliothèque standard Python uniquement** :

- vérification du nombre de sources, rôles et SHA source, nombre d'acquisitions et concordance du nombre de features entre la preuve R, la matrice exportée et les enregistrements CSV ;
- recalcul des effectifs `RSD>0.30`, `ratio>=0.20`, des manquants et de l'évaluabilité avec un contrôle strict du calcul feature-par-feature et des flags du moteur R ;
- conservation du nombre de métriques **non évaluables**, distinction explicite des mesures manquantes et signaux nuls ;
- contrôles négatifs où sont falsifiés une valeur QC, un ratio, la signature source, une classe de prélèvement ou un statut de certification.

Un second workflow **ne se contente pas de réimprimer le JSON R** : il refuse de qualifier comme favorables les features à blancs non mesurés, et génère un rapport `scientific_QC_review_required` sans supprimer automatiquement une feature. Les vraies données d'origine ne sont jamais copiées dans l'artefact public ; seuls rapports de source et indicateurs de features anonymes sont archivés.

## Décision pour le chercheur

**État du sous-ensemble Sacurine : qualification de l'intégrité et de l'exécution instrumentale démontrée ; revue analytique obligatoire avant l'inférence biologique.**

Les trois QC sont trop peu nombreux et trop éloignés dans la séquence originale pour valider une **correction** de dérive. Il n'y a qu'un seul lot ; aucune comparaison inter-lots n'est qualifiée. Le workflow ne démontre ni identité moléculaire, ni reproductibilité interlaboratoire, ni contrôle du FDR de tests biologiques, ni relation individuelle omics→PK. Les spectres d'origine contiennent des prélèvements humains ; seul le matériau de provenance et de QC prévu est conservé dans l'artefact.

Cette étude fournit une preuve externe utile pour **un article de logiciel** et pour vérifier le comportement de MultiOmics sur une matrice instrumentale mesurée. Elle ne constitue **pas** une validation clinique ou universelle du traitement LC-MS.
