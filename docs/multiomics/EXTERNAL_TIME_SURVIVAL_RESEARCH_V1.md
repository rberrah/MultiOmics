# Research v1 — qualification externe réelle des voies longitudinales et survie

**Protocole figé le 19 octobre 2026 avant lecture des résultats feature-wise.**

Cette étape répond à l'issue #61. Elle ne cherche ni un biomarqueur favorable ni une
réplication biologique sélectionnée a posteriori. L'objectif est plus étroit :
vérifier sur **de vraies matrices omiques publiques et de vrais endpoints/métadonnées**
que les tables natives projetées par MultiOmics reproduisent exactement des ajustements
R écrits indépendamment, tout en conservant les diagnostics défavorables.

## Survie — NKI

Source : package Bioconductor `breastCancerNKI`, cohorte NKI publiée par van't Veer /
van de Vijver. Endpoint fixé avant calcul : **DMFS** via `t.dmfs` et `e.dmfs`.

- seuls les sujets avec temps DMFS positif et statut 0/1 complet sont inclus ;
- aucune covariable n'est choisie après inspection des résultats ;
- famille de test : les **64 premiers identifiants de features par ordre lexical**
  ayant des mesures finies et non constantes dans la population DMFS ;
- `run_cox_survival_reference(..., ties="efron")` est comparé feature par feature
  à un `survival::coxph` écrit séparément ;
- coefficient, erreur standard, p, BH, HR, IC sur beta et `cox.zph` doivent concorder ;
- une violation PH reste visible et ne fait jamais disparaître la feature.

Cette preuve établit une **parité d'exécution sur données réelles**, pas la validité
pronostique de la famille de 64 features ni un contrôle universel du FDR.

## Longitudinal — GSE191328

Source : GEO **GSE191328**, cohorte découverte infliximab de Mishra et al.
(Genome Medicine 2022, PMID 36153599). Le dépôt GEO documente le patient,
la cohorte, le traitement, le statut de rémission à S14 et le temps de prélèvement.

Design fixé avant calcul :

- traitement : discovery–infliximab uniquement ;
- groupes stables : **rémission S14 vs non-rémission S14** ;
- visites : **baseline et 2 semaines** uniquement ;
- seuls les patients possédant exactement ces deux prélèvements RNA-seq sont éligibles ;
- les comptages featureCounts publics sont téléchargés individuellement et hashés ;
- preprocessing avant modèle mixte : `edgeR::filterByExpr`, TMM, logCPM
  (`prior.count=2`) ; les comptages bruts ne sont donc jamais présentés à
  MultiOmics comme une matrice confirmatoire longitudinale ;
- famille : 64 premiers identifiants Ensembl par ordre lexical après filtrage
  d'expression, sans regarder les coefficients ou p-values ;
- deux visites impliquent par construction un **random intercept** ; le coefficient
  groupe × temps de `run_lmer_matrix` est comparé au même modèle
  `lmerTest::lmer(value ~ condition*time + (1|subject_id), REML=FALSE)`
  écrit séparément.

Le statut de rémission est observationnel. Cette qualification ne transforme pas
le contraste en effet causal du traitement, ne valide pas la linéarité de trajectoires
plus longues et ne certifie ni puissance ni calibration fréquentiste.

## Critère Research v1

Les deux jobs doivent :

1. partir de sources publiques nommées et refuser une modification majeure du contrat source ;
2. préserver la famille de features préspécifiée complète ;
3. concorder numériquement avec les ajustements R indépendants ;
4. archiver les diagnostics et résultats négatifs ;
5. déclarer explicitement `statisticalCalibrationCertified=false`,
   `biologicalValidation=false` et `clinicalCertification=false`.

Un PASS permet de dire : **« les voies Cox et longitudinales natives sont qualifiées
sur une exécution externe réelle pour ces contrats précis »**. Il ne permet pas
de dire : « MultiOmics est universellement confirmatoire ».
