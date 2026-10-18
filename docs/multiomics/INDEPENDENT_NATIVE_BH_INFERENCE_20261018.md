# MultiOmics — qualification indépendante de l'inférence RNA-seq

**Date : 18 octobre 2026.** Protocole fixé avant les exécutions de cette
branche. **Aucun choix de méthode RNA-seq, seuil patient, modèle de
production ou correction a posteriori n'est modifié.**

## Pourquoi cette qualification supplémentaire ?

Les campagnes sur des données négatives et à vérité mixte avaient déjà
montré que « la CI est verte » ne signifie pas « FDR biologiquement
contrôlé à 5 % ». Les essais RNA-seq Bioconductor indépendants étaient
exécutés, mais leur contrôle automatique portait surtout sur le nombre
global de découvertes, la puissance, le sens des effets et une alarme
grossière (FDR simulé > 0,35).

Cette campagne teste **une autre hypothèse** : les valeurs ajustées q
restituées par les vraies méthodes sont-elles numériquement conformes
à l'ajustement BH de leurs *p-values brutes sélectionnées pour BH* ?
C'est un contrôle de **l'intégrité de l'inférence et des tableaux**.
Un BH correctement calculé n'implique aucunement des p-values bien
calibrées sous une mauvaise hypothèse ou un mauvais plan.

## Trois validations distinctes — critères avant analyse

| Couche | Contrôle | Portée permise |
| --- | --- | --- |
| A — intégrité numérique R | Ajustement BH recalculé par **formule step-up indépendante de `stats::p.adjust`** depuis les p-values **de chaque résultat R d'origine** : DESeq2, edgeR-QL et limma-voom ; identifiants identiques à ceux du jeu simulé, absence de doublon, tolérance absolue 1e-8 | Intégrité du tableau q pour ce jeu, **pas** validation du modèle |
| B — oracle algorithmique | 2 048 jeux indépendants × 200 p-values nulles Uniform(0,1), alpha=0,05, graine 20261018, BH selon une implémentation indépendante, confrontée à `stats::p.adjust`, probabilité de ≥1 découverte et intervalle binomial exact 95 % | Contrôle positif du comportement BH sous des p-values uniformes **imposées**, pas générées par un modèle RNA-seq |
| C — témoin d'échec | 2 048 jeux indépendants × 200 p-values sous Beta(0,5 ; 1), toutes déclarées nulles, BH indépendant, même alpha | Doit montrer une inflation forte ; atteste que le contrôle n'est pas un PASS automatique |

Critères informatiques préfixés :

- A : écart maximal entre BH indépendant et q natif **≤1e-8**,
  identifiants/longueurs identiques aux sources simulées, aucune valeur
  ajustée finie sans p-value finie. Un défaut fait échouer la CI.
- B : fréquence empirique de découverte au moins une fois comprise dans
  **0,05 ± 0,03**. Cette tolérance grossière de non-régression n'est
  **pas** un test certifiant exactement 5 %.
- C : taux de découverte faux **>0,20**, afin qu'une rupture volontaire
  d'uniformité déclenche un résultat négatif visible. Tous les résultats
  et IC binomiaux sont archivés, y compris en cas d'échec.
- A/B/C : aucune sélection d'un « meilleur » logiciel sur ces résultats
  ni adaptation de paramètres à leurs performances observées.

## Précaution scientifique importante : le filtrage DESeq2

DESeq2 peut produire une p-value **finie** mais une `padj=NA` après
filtrage indépendant. Le contrôle B-H ne doit donc **pas** recalculer les
q-values sur la totalité des p-values brutes disponibles. Il reconstruit
BH par la formule mathématique sur les p-values dont la q-value native est finie, et archive
`excluded_raw_p_with_no_q`. Ce détail ne démontre pas que le
filtrage indépendant est sans effet sur le FDR : cette question reste
ouverte et nécessite ses propres simulations.

Les comparaisons par méthode emploient leurs propres sorties natives :
`pvalue/padj` (DESeq2), `PValue/FDR` (edgeR-QL) et
`P.Value/adj.P.Val` (limma-voom). Les mêmes matrices de simulation,
étiquettes et contrastes du benchmark existant sont conservés.

## Exécution réelle et artefacts attendus

La qualification figée R 4.5.2 / Bioconductor 3.22 exécute :

```bash
Rscript multiomics-engine/tests/test_reference_inference_audit.R
Rscript multiomics-engine/validation/rnaseq_bh_oracle.R tmp/reference-methods
Rscript multiomics-engine/tests/bioconductor_null_calibration.R
Rscript multiomics-engine/tests/bioconductor_mixed_truth_reference.R
```

Dans `tmp/reference-methods/` :

- `bh_oracle_independent_null_per_study.csv` et
  `bh_oracle_independent_null_summary.csv`,
  plus le manifeste comportant les véritables versions R ;
- `bioconductor_null_native_bh_integrity.csv` ;
- `bioconductor_mixed_truth_native_bh_integrity.csv` ;
- les tables originales de calibration négative et vérité mixte,
  jamais écrasées par les données d'oracle.

Les intervalles Monte-Carlo décrivent l'incertitude *entre simulations*.
L'oracle Uniform(0,1) ne teste pas l'estimation NB, la dispersion, la
normalisation, les lots, les petits effectifs, le MNAR ni le plan des
cohortes réelles. Les campagnes natives à seulement huit répétitions
par scénario restent des **détecteurs d'erreurs grossières**.

## Décision méthodologique

**Autorisé après succès :** documenter l'exécution réelle de trois
méthodes, la traçabilité numérique BH, les vrais q-values par scénario
simulé, les éventuels faux appels et une référence négative démontrant
que le contrôle statistique sait échouer.

**Interdit :** « FDR universel à 5 % », « analyse confirmatoire
certifiée », « biomarqueur validé », « résultats cliniques qualifiés »,
« supériorité de la méthode qui a le moins de faux appels ».

Étapes ultérieures : davantage de répétitions indépendantes avec
précision Monte-Carlo préspécifiée et coûts calculatoires documentés,
générateurs externes au modèle NB, splits sujets réels, contrôles
indépendants sur les designs appariés et les effets de lot, puis validation
biologique multi-cohorte. La production n'est pas modifiée.
