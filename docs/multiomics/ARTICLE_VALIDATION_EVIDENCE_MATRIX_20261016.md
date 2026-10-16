# MultiOmics — état des preuves pour un éventuel article utilisant la plateforme

**Date : 16 octobre 2026.** Ce tableau est un inventaire de preuves disponibles, non une revendication de validation confirmatoire universelle ni la rédaction d'un article sur l'outil. Pour publier une **analyse réalisée avec MultiOmics**, la première cible est un rapport *Methods/Results* régénérable et attaché à la cohorte effectivement étudiée.

## Matrice de qualification (preuve testée, interprétation limitée)

| Domaine | Source indépendante et protocole | Ce qui est établi | Ce qui n'est PAS établi |
| --- | --- | --- | --- |
| RNA-seq / appariement | **Airway / GSE52778**, vrai DESeq2 R/Bioconductor, 17 199 lignes, ZIP/Rmd et 30 × 8 valeurs visualisées; [CI Airway](https://github.com/rberrah/MultiOmics/actions/workflows/public-airway-r-reference.yml) | Rejeu numérique du design et des tables spécifiés, provenance code/R et vérification de l'export scientifique | Calibration de tous les designs, effet clinique, transférabilité |
| Intégration multi-omique | **Nutrimouse** RNA/lipides, 40 animaux, MOFA2 et DIABLO réels, contrôle des paires, permutations négatives, ablations mono-omiques | Exécution mécanisme/statistique de méthodes biologiquement réelles et diagnostics de facteurs/erreurs | Supériorité générale du multi-omique, indépendance inter-cohortes, biomarqueurs validés |
| LC-MS blancs | **Sacurine Zenodo 1346742** : deux urines et deux blancs, vrais profils xcms (étude séparée) | Vérification de la soustraction/détection des signaux de fond pour les acquisitions disponibles | Plage validée d'erreur, QC temporel, sélectivité de marqueurs |
| LC-MS QC | **Sacurine Zenodo 3244991 / Galaxy GTN** : six urines + trois vrais QC poolés, données mzML d'origine; [PR #18](https://github.com/rberrah/MultiOmics/pull/18), [CI 38953745195](https://github.com/rberrah/MultiOmics/actions/runs/38953745195) | MD5 sources originaux, 9 assays xcms, **101 526 pics chromatographiques**, **3 347 features regroupées**, **2 957 features avec RSD QC évaluable**, alignement RDS/CSV | Aucun blanc **dans ce jeu**, ordre d'injection non vérifié, pas de conclusion sur la dérive, pas de biomarqueur ou FDR externe |
| Calibrage RNA-seq | Générateurs nuls / mixed truth DESeq2, edgeR et limma ; vraie exécution Bioconductor R figée | Diagnostics des faux appels, erreurs de direction, puissance dans les scénarios prédéfinis | **Contrôle FDR 5 % universel non établi** ; certains scénarios avaient des faux appels excessifs, précision Monte-Carlo limitée |
| Reproductibilité du ZIP | Vérification SHA-256, sources originales, CSV R exhaustifs quand disponibles, replay Rmd, comparaison navigateur/R | Capacité de vérifier certains calculs avec un moteur réel et une archive déclarée | La simple présence d'un script ne prouve pas un modèle R identique ou une interprétation biologique |
| Interopérabilité PK | ZIP d'origine MultiOmics → lecteur indépendant **Omics2PK figé**; jeux de démonstration *synthétiques* avec 8 sujets, 16 prélèvements, 49 analyses et T0/T12 conservés; [PR #21](https://github.com/rberrah/MultiOmics/pull/21) | Sérialisation/import, hachage, identité sample/assay/subject, préservation des temps catégoriels et refus des scores de groupe | **Aucun effet PK réel**; pas de score gelé validé sur une cohorte clinique indépendante |
| Accessibilité, usage quotidien | Parcours navigateur Playwright FR/EN desktop/mobile ; quatre démonstrations, CSV MS, export | Fonctionnement du parcours automatisé sur les versions testées | Pas de test utilisateur indépendant ni d'évaluation empirique d'utilisabilité |

Ces études ne forment **pas** une cohorte combinée. Surtout, les neuf acquisitions QC et les quatre acquisitions blancs de Sacurine viennent de **deux enregistrements distincts**; ne jamais les mélanger pour estimer contamination, dérive ou FDR.

## Critères minimaux avant une analyse scientifique réutilisable dans un article

1. **Design** : hypothèse, contraste et direction déclarés avant la modélisation ; unité biologique indépendante (sujet, culture, animal), modalités et répétitions techniques identifiables. Éviter pseudoréplication et ajustement de covariables non identifiable.
2. **Vraies données et autorisations** : accession publique/consentement, schéma de données, dates/échelles originales, statut des manquants, prétraitements et métadonnées expérimentalement documentés, sans inférence à partir des noms.
3. **Méthode exécutée** : préciser distinctement Svelte/JS exploratoire, R DESeq2/limma/xcms/fgsea/MOFA2/DIABLO effectivement exécuté, paramétrage, versions R/Bioc/paquets et random seeds ; confirmation analytique uniquement dans le domaine scientifiquement qualifié.
4. **Comparateur indépendant** : jeux publics externes non utilisés pour adapter les hyperparamètres, comparaison de résultats natifs ligne par ligne ; pour MS, QCs + blancs + ordre d'injection authentique **dans le même batch** si l'on veut parler de dérive/contamination.
5. **Statistics** : conserver toutes les features testées, les exclusions, tailles d'effet et leurs incertitudes, p-values et corrections multiplicité, stratégies NA ; simulations nuls + power avec intervalles Monte-Carlo préspécifiés ; aucune correction post hoc choisie après inspection des mêmes données.
6. **Multi-omique** : splits par sujet et, si possible, par cohorte ; prétraitements ajustés uniquement sur train ; baselines mono-omiques et permutation d'appariement ; aucune inférence de mécanisme causal depuis des facteurs latents.
7. **Sorties auditables** : archive input+méthodes+CSV exhaustifs+vraies traces R+SHA-256+script de rejeu ; vérification sur machine et environnement séparés, contrôle des erreurs et des runs bloqués.
8. **Validité extérieure et lecteur** : réplication sur une cohorte effectivement indépendante et protocole novice/biostatisticien pour vérifier la compréhension des statuts et l'export si publication du logiciel elle-même.

## Priorités scientifiques après la PR #18

**Priorité S1 : compléter l'étude LC-MS sur une série unique documentant vrais blancs, QC poolés, biologique et ordre d'injection**, puis comparaison pré-déclarée à un workflow externe tel que Workflow4Metabolomics avec métriques de correspondance des features (m/z, rt), distributions QC RSD (avec incertitudes), taux de valeurs manquantes et reproductibilité de l'export.

**Priorité S2 : valider prospectivement le contrôle de la calibration**, plutôt que d'augmenter simplement le nombre de scénarios. À chaque scénario figer l'hypothèse nulle, R, N, dispersions et déséquilibres, `N_simulations`, seuil de décision et incertitude avant calcul. Conserver et publier les résultats défavorables.

**Priorité S3 : montrer une étude d'usage reproductible de bout en bout** sur une vraie cohorte accessible : import du jeu public, analyse complète, tableaux R natifs (si moteur réellement exécuté), figure, méthodes et ZIP vérifié indépendamment. L'absence de données externes nécessaires doit apparaître comme **BLOCKED**, et non comme une étude « PASS ».

### Tableau de décision pour rédiger

- **Article utilisant MultiOmics** : raisonnable si l'étude biologique choisie respecte le contrat ci-dessus, avec données réelles, comparateurs et traces indépendantes.
- **Article logiciel / note technique** : architecture + provenance + dossiers reproduits + études positives/négatives déjà possibles à présenter de manière strictement bornée ; pas de prétention « plateforme confirmatoire universelle ».
- **Article de validation scientifique large** : **non prêt**. Manquent, en particulier, la généralisation multi-laboratoire, un contrôle expérimental MS complet, validation externe des signatures et des tests utilisateurs.

## Publication logicielle et site public

Les expériences se qualifient sur le **SHA exact**. PMx Explain monte le paquet statique d'un **autre** commit source figé ; son déploiement doit être validé séparément. Les six workflows scientifiques/logiciels d'une release candidates sont exécutés sur tout push `main` sans filtrage de chemins, afin qu'un changement d'interface ne contourne pas ou ne bloque pas la qualification de release. Cette mesure n'enlève ni l'obligation de signaler les échecs ni la vérification du site réellement publié.
