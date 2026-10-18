# MultiOmics — protocole de précision Monte-Carlo du FDR RNA-seq

**Préenregistrement logiciel : 18 octobre 2026, avant les calculs des nouveaux jeux.**

## Question posée

Les qualifications antérieures ont identifié des taux de faux positifs
potentiellement supérieurs au seuil nominal dans certains simulateurs RNA-seq
à loi binomiale négative, mais les répétitions des tests quotidiens (huit
par scénario) sont trop peu nombreuses pour une estimation précise. Le
contrôle numérique Benjamini–Hochberg de la PR #42 ne permet **pas** à lui
seul d'attribuer ou d'exclure ce problème.

Cette étude augmente les **répétitions de jeux de données indépendants**,
sans changer le simulateur, les méthodes R, les effets plantés, les plans,
les règles de sélection ou les seuils. **La CI ne sera pas interprétée comme
un certificat de FDR à 5 %.**

## Plan et horizon fixes (non adaptatifs)

| Famille | Trois scénarios figés | Jeux *par scénario* | Méthodes sur les mêmes jeux |
| --- | --- | ---: | --- |
| Nul global | équilibré standard ; variation profondeur + dispersion ; lots partiellement confondus avec groupe 12/28 | 128 | DESeq2, edgeR-QL, limma-voom |
| Vérité mixte | 1 %, 5 %, 20 % de gènes à effet planté +/−1,5 log2FC | 64 | DESeq2, edgeR-QL, limma-voom |

Soit **384+192=576 jeux biologiques synthétiques nouveaux**, et
**1 728 fits R** (trois méthodes par jeu). Le générateur et les scripts
de benchmark existants sont réutilisés, en ajoutant seulement les
*décalages de répétitions* ci-dessous :

- Nul global : `MULTIOMICS_NULL_REPLICATE_OFFSET=1000`,
  indices 1001–1128, graines `780000 + scenario_index × 10000 + index`.
- Vérité mixte : `MULTIOMICS_MIXED_REFERENCE_REPLICATE_OFFSET=1000`,
  indices 1001–1064, graines `620000 + scenario_index × 10000 + index`.
- Aucun indice ne recouvre 1–8 de la CI figée. Les graines sont
  distinctes entre générateurs et scénarios, identiques **entre méthodes**
  pour une même simulation (appariement correct).
- Aucune répétition additionnelle conditionnée par les résultats n'est
  autorisée. Tout fit manquant ou graine dupliquée invalide le rapport.

Version R, Bioconductor et packages exacts contrôlés comme dans la CI figée
(R 4.5.2, Bioconductor 3.22, image Docker épinglée par SHA-256).

## Estimands et incertitude

### Nul global

L'unité indépendante est un **jeu simulé complet**. Pour chaque méthode et
scénario, le FDR empirique sous le nul global est la fréquence des jeux
ayant au moins une découverte BH à q≤0,05 (FWER = FDR sous le nul global).

- Estimation : `k/128`, avec `k` = nombre de jeux comportant au moins
  une fausse découverte.
- Intervalle individuel exact Clopper–Pearson 95 %.
- Intervalle **simultané familial** 95 % sur les neuf cellules
  scénario × méthode, via Clopper–Pearson au niveau
  `1 − 0,05/9` pour chacune. Bonferroni reste valide malgré la
  dépendance des méthodes sur les mêmes jeux.
- Le statut `simultaneous_inflation_detected` signifie :
  borne inférieure simultanée strictement >0,05 dans **ce générateur**.
  Si l'intervalle contient 0,05 : `inconclusive`, non « validé ».
- Ce protocole n'est pas un test de supériorité interméthodes, et aucune
  méthode ne sera choisie a posteriori sur ces neuf cellules.

### Vérité mixte

- FDR simulé : `mean(V/max(R,1))` entre les **64 jeux indépendants**
  de chaque cellule, non pas un pool de gènes assimilés à des réplications.
- Puissance : fraction des effets biologiques réellement plantés détectés,
  dénominateur fixé par le générateur, même si un gène est filtré.
- IC descriptifs bootstrap à 95 % rééchantillonnant des études entières,
  graine fixée, ainsi que bornes simultanées Hoeffding–Bonferroni sur
  les neuf cellules. Cette borne est très prudente à N=64 et sa largeur
  ne doit pas être masquée.

Un passage « vert » signifie exécution numérique correcte et contrôles
d'intégrité respectés. **Des taux FDR élevés doivent être conservés
tels quels, ne pas être corrigés par un seuil inventé.** Le test logiciel
conserve seulement le garde-fou grossier hérité du benchmark ; il ne
certifie jamais le taux nominal à 5 %.

## Calcul et artefacts

Fichier autonome :
`.github/workflows/rnaseq-monte-carlo-precision.yml`.
Le workflow se déclenche en PR sur les scripts pertinents, lors de la
fusion vers `main` sur les mêmes chemins, ou manuellement. Il est distinct
des huit workflows usuels afin d'éviter de facturer 1 728 fits à chaque
modification d'interface.

Il lance dans cet ordre la validation des fonctions sans Bioconductor,
les versions R figées, les 384 jeux nuls, les 192 jeux à vérité mixte,
puis la validation *fail-closed* du nombre de jeux et des graines.

Artefacts d'intérêt :

- `bioconductor_null_replicates.csv` et
  `bioconductor_mixed_truth_per_study.csv` : tous les jeux, graines,
  scénarios et découvertes ;
- `bioconductor_null_native_bh_integrity.csv` et
  `bioconductor_mixed_truth_native_bh_integrity.csv` : contrôle indépendant
  des q-values retournées par les vrais packages ;
- `rnaseq_precision_null_simultaneous_95.csv`,
  `rnaseq_precision_mixed_95.csv` : incertitude Monte-Carlo ;
- `RNASEQ_MONTE_CARLO_PRECISION_REPORT.md` et
  `rnaseq_precision_run_manifest.txt` : synthèse, limites, SHA et méthodes ;
- `bioconductor_*_environment.txt` : versions R/Bioconductor et prémisses.

Un échec du workflow ne justifie jamais la suppression d'un résultat
scientifique défavorable. L'artefact `if: always()` conserve les données
partielles après une erreur.

## Non démontré par cette campagne

La précision Monte-Carlo croissante ne garantit pas que le simulateur NB
soit le mécanisme approprié aux données réelles. Restent nécessaires :
générateurs hors NB, étude du MNAR, designs biologiques appariés, effets
de lot expérimentaux, données externes LC-MS mesurées (blancs/QC/ordre
d'injection documentés), reproductibilité interlaboratoire et revue
méthodologique des analyses individuelles. Aucun argument de marquage CE
ou d'utilisation clinique n'est avancé.
