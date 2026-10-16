# RNA-seq R1 — dossier de qualification du parcours confirmatoire

**Statut (11 octobre 2026) : garde-fous du mode confirmatoire en cours de
validation logicielle ; pas de certification statistique, biologique ou
clinique.** Ce document précise ce qui est réellement vérifiable, ce qui
reste à prouver, et quand le pipeline doit refuser le mode confirmatoire.

## Domaine de validité V1 restreint

- Type d'échantillon explicitement documenté comme biologique.
  Les blancs, QC poolés et injections techniques hors sujet ne sont
  **jamais supprimés silencieusement** dans le parcours confirmatoire :
  la présence de ces éléments bloque l'analyse jusqu'à constitution
  d'un manifeste d'exclusions préspécifié et vérifiable.
- Une seule omique : **comptages RNA-seq bruts non négatifs et entiers**.
  La forme numérique ne prouve pas l'origine réelle des comptages :
  provenance FASTQ/BAM, alignement, version annotation et comptage
  nécessitent des preuves supplémentaires.
- Deux groupes **biologiquement indépendants**, au moins trois sujets
  distincts par groupe. Cette limite est un seuil de prévention des
  erreurs de logiciel, absolument **pas** un seuil de puissance.
- Le navigateur doit transmettre réellement `analysisIntent` au moteur R.
  En mode confirmatoire, le chercheur choisit explicitement le groupe
  comparé et le groupe de référence. L'ordre alphabétique des noms
  n'est **jamais** utilisé pour définir le sens scientifique du log2FC.
  L'absence du choix, son ambiguïté ou des niveaux inconnus bloquent
  l'ajustement R ; les choix sont conservés dans l'export reproductible.
- Hypothèse d'intérêt : différence de condition, un contraste explicitement
  nommé et préspécifié. Régression additive avec intercept, covariables
  observées, sans interaction, appariement, cluster ou série temporelle.
- Un sujet biologique unique par `sample_id`. Un sujet réparti entre
  plusieurs lignes considérées à tort comme indépendantes est bloqué.
  Les réplicats techniques du même prélèvement peuvent être agrégés
  uniquement lorsque leurs `subject_id`, `condition` et `batch`
  sont cohérents.
- **DESeq2 de référence effectivement exécuté** avec audit strict des
  identifiants d'analyses (assays) *avant* agrégation et audit des
  entrées et de l'identifiabilité *avant* ajustement statistique.
- Correction multiple documentée par famille de tests ;
  le résultat ne devient jamais automatiquement une hypothèse confirmée.

## Contrôles logiciels automatisés intégrés

| Étape | Rejet de l'analyse demandée si… | Trace |
| --- | --- | --- |
| CSV → métadonnées | Un `assay_id` manque, est dupliqué ou n'est pas présent dans les mesures | `assayIntegrityContract` |
| Réplicats techniques | Deux lignes d'un même prélèvement indiquent des sujets, conditions ou lots incompatibles | Refus avant agrégation |
| Sujets indépendants | Un `subject_id` est absent ou appartient à plusieurs échantillons considérés indépendants | `confirmatoryInputContract` |
| Mesures RNA | Valeurs manquantes, négatives, fractionnaires, bibliothèque vide, IDs de features dupliqués | Refus strict |
| Jointure | Les ensembles d'IDs RNA et métadonnées ne correspondent pas exactement | Refus, aucune intersection silencieuse |
| Contraste et plan | Moins de deux groupes, moins de trois sujets indépendants par groupe, formule non additive ou contraste ambigu | Refus strict |
| Covariables | Batch variant omis du modèle, covariable manquante, alias ou défaut de rang, trop peu de degrés de liberté | Refus strict |
| Intégration navigateur | Une exécution R réussie est annoncée sans les **deux** contrats précédents | Affichage `blocked` |
| Conclusion | Les contrôles passent | **Éligible pour revue indépendante uniquement** ; `certified=false` |

Les tests sont dans
`multiomics-engine/tests/test_rnaseq_confirmatory_preflight.R`,
`scripts/test_multiomics_diagnostics.mjs` et le workflow
`.github/workflows/multiomics-rnaseq-confirmatory.yml`.

## Matrice de validation statistique R1 — reste à réaliser

| Scénario préspecifié | Estimands et preuves attendues | État |
| --- | --- | --- |
| Hypothèse nulle, 2 groupes équilibrés | Proportion de jeux avec ≥1 découverte BH (FDR=FWER sous nul global), IC binomial exact à 95 % | 12 répliques exécutées sur Bioconductor : 0/12 découvertes pour les 3 moteurs ; intervalle trop large pour conclure à 5 % |
| Nul avec profondeur variable et forte dispersion | Groupes équilibrés ; profondeur ~10× ; dispersion NB élevée ; BH à 5 % | 8 jeux exécutés : DESeq2 1/8, edgeR 0/8, voom 0/8 ; aucune validation nominale revendiquée |
| Nul avec répartition des groupes déséquilibrée | 12 témoins / 28 traités ; lot technique associé au groupe mais non complètement confondu | 8 jeux exécutés : DESeq2 2/8, edgeR 1/8, voom 0/8 ; puissance et couverture restent non établies |
| Vérité mixte : 1/5/20 % de gènes modifiés | FDP par jeu, FDR moyenne, puissance, biais log2FC, nombre de gènes conservés | 3 scénarios de 12 jeux exécutés ; FDR estimés parfois > 5 %, revue méthodologique requise |
| Effets faibles, forts et de signes opposés | Sensibilité, estimation et couverture des IC selon effectifs | À étendre |
| Batch mesuré, partiellement confondu | Exactitude du contraste / vérification du rang, variations de FDR | À étendre |
| Batch parfaitement confondu | Refus explicite, **aucune estimation identifiable** | Préflight ajouté |
| Sujet dupliqué, étiquettes mélangées, rang dégénéré | Échec explicite, aucun calcul concluant silencieusement | Préflight ajouté |
| Cohorte RNA-seq publique externe, non inspectée auparavant | Reproduction d'un effet primaire préspécifié, prédictions/contrastes, QC indépendants | Non établi |
| Reproductibilité dans un environnement verrouillé | Docker/renv, versions R/DESeq2/annotation, checksums, réanalyse byte-identique ou tolérances numériques | Pas encore garanti |

### Comment quantifier honnêtement l'incertitude du FDR ?

**Problème scientifique identifié et corrigé :** si 8 ou 12 simulations indépendantes ne donnent aucune fausse découverte, le bootstrap percentile du FDR renvoie un intervalle **[0 %, 0 %]** parce qu'il ne peut rééchantillonner que les zéros observés. **Cet intervalle ne démontre absolument pas un FDR réel nul ou maîtrisé.** Il est conservé comme descriptif, mais ne sert plus à décider de la maîtrise nominale.

L'audit R `multiomics-engine/mixed_truth_reference_metrics.R` fournit désormais une **borne simultanée de Hoeffding, avec correction de Bonferroni**, en complément du bootstrap. Pour chaque scénario et méthode, elle encadre la moyenne du FDP entre **études indépendantes** à niveau de confiance familial 95 %, en tenant compte de toutes les cellules scénario–méthode du même rapport. Cette garantie mathématique utilise uniquement le fait que chaque FDP se trouve entre 0 et 1 et que les études simulées sont **indépendantes et identiquement distribuées dans chaque scénario fixé**. La corrélation entre méthodes appliquées aux mêmes études est autorisée via la correction simultanée.

- Avec un faible nombre de simulations, la borne supérieure reste large **même lorsque toutes les valeurs de FDP observées sont nulles**. C'est voulu : l'absence de fausses découvertes observées ne constitue pas une calibration confirmée.
- Un intervalle bootstrap de `[0 %, 0 %]` doit être affiché comme **descriptif insuffisant**, jamais comme « méthode validée à 5 % ».
- Les indicateurs d'alerte du stress test utilisent maintenant les **bornes conservatrices** et distinguent `fdr_calibration_unresolved_insufficient_independent_studies` d'une borne inférieure au-dessus du seuil.
- La borne porte exclusivement sur la **distribution artificielle simulée** ; elle ne garantit rien pour une nouvelle cohorte réelle, ni pour une simulation sélectionnée après inspection, ni pour tous les mécanismes d'expression biologique.

Ce garde-fou n'obtient pas magiquement une validation statistique : il **empêche une preuve trompeuse par insuffisance de simulations**. La précision suffisante devra être fixée **avant** une campagne Monte-Carlo indépendante plus vaste et mise en relation avec des jeux de données publics.

### Résultats R effectivement observés le 11 octobre 2026

Deux workflows GitHub Actions ont **exécuté** les moteurs Bioconductor de référence sur les données simulées. Les résultats ci-dessous sont des **mesures de simulation** et non des taux d'erreur attendus pour une vraie cohorte. Les logs, les sorties par simulation et les versions R ont été conservés dans les artefacts des workflows.

**Nul global (au moins une découverte BH q ≤ 0,05 par jeu)** — [workflow 38173002502](https://github.com/rberrah/rberrah.github.io/actions/runs/38173002502) :

| Scénario (jeux indépendants) | DESeq2 | edgeR QL | limma-voom |
| --- | ---: | ---: | ---: |
| 20/20, dispersion modérée (12) | 0/12 | 0/12 | 0/12 |
| 20/20, forte dispersion et profondeur variable (8) | 1/8 | 0/8 | 0/8 |
| 12/28, lot technique partiellement confondu (8) | 2/8 | 1/8 | 0/8 |

**Vérité mixte (FDR moyen = moyenne des proportions de fausses découvertes par jeu)** — [workflow 38174106412](https://github.com/rberrah/rberrah.github.io/actions/runs/38174106412) :

| Vrais gènes modifiés parmi 400 (12 jeux par cas) | DESeq2 | edgeR QL | limma-voom |
| --- | ---: | ---: | ---: |
| 4 (1 %) | 11,11 % | 8,33 % | 6,67 % |
| 20 (5 %) | 6,76 % | 3,36 % | 2,64 % |
| 80 (20 %) | 5,45 % | 4,33 % | 3,47 % |

Tous les gènes modifiés plantés ont été détectés dans ces scénarios **à forts effets**, mais ce plafond de puissance reflète la facilité artificielle de l'expérience, pas une garantie pour les études réelles. Pour DESeq2, les couvertures approximatives moyennes des IC Wald non rétrécis sur les vrais signaux sont respectivement 100 %, 92,5 % et 93,65 % selon ces trois scénarios ; il s'agit d'estimations très peu précises (12 jeux, 4/20/80 signaux). **L'observation de FDR au-dessus de 5 % dans plusieurs cellules demeure un signal de revue méthodologique**, même si le seuil logiciel d'alerte à 35 % n'a pas été dépassé. Le test ne permet **ni** de conclure à une inflation réelle du FDR **ni** d'établir sa maîtrise nominale.

**Action requise avant qualification confirmatoire :** campagne indépendante bien plus grande avec précision Monte-Carlo préspecifiée, effets faibles et hétérogènes, tailles d'échantillon variées, plan confondu partiellement et facteurs techniques réalistes ; conserver les résultats négatifs et figer les versions logicielles et critères de décision avant de réexaminer les performances.

### Résultats observés de l'audit de filtrage (12 octobre 2026)

L'[exécution Bioconductor 38214629245](https://github.com/rberrah/rberrah.github.io/actions/runs/38214629245) a réussi **8 études simulées indépendantes par scénario** avec les **trois moteurs R** et deux politiques de préfiltrage. Les sorties détaillées sont archivées sous `rnaseq-paired-prefilter-reference-validation`. Ces chiffres ne sont **pas** des critères de validation externe ni des propriétés universelles de DESeq2, edgeR ou limma-voom.

| Scénario (400 gènes, 8 cohortes) | Gènes écartés par le filtre plus strict | Vrais signaux écartés parmi ceux plantés | Variation de puissance DESeq2 / edgeR / voom (points de pourcentage) |
| --- | ---: | ---: | --- |
| `small_cohort_6v6` | 190,1 | 9,6 / 20 | 0 / 0 / 0 |
| `sparse_signals_20v20` | 196,3 | 1,9 / 4 | −12,5 / −12,5 / −3,1 |
| `unequal_batch_12v28` | 194,8 | 9,6 / 20 | −0,6 / −2,5 / −1,3 |

La variation rapportée est **filtre strict − analyse de référence**, sur les **mêmes cohortes**. L'absence de variation de puissance dans la petite cohorte ne signifie pas que les prétraitements sont équivalents : la puissance y est déjà extrêmement basse. Sous signaux rares, le filtrage plus strict peut faire disparaître des gènes réellement modifiés et ainsi réduire sensiblement le nombre de découvertes vraies.

Le protocole a également quantifié la variation moyenne du facteur d'échelle effectif sur les échantillons, après centrage en log2 ; elle dépend du scénario et du moteur. **Ce changement de normalisation fait partie du contraste**, il n'a pas été artificiellement tenu fixe. Certains FDP et FDR diminuent sous filtrage strict, mais une diminution peut simplement résulter d'une absence de découverte (FDP = 0 par convention) et s'accompagner d'une perte de puissance. Il ne faut donc ni sélectionner un filtre parce que son FDR observé est bas ni classer les méthodes sur huit simulations.

**À approfondir avant qualification scientifique :** effets de composition, gènes à zéro dans tous les échantillons, sérialisation des paramètres exacts de normalisation, tests de sous-groupes, grand nombre de répliques Monte-Carlo et cohortes publiques indépendantes. La sensibilité documentée sert à **détecter les fragilités de l'analyse**, pas à optimiser son résultat après observation.

### Premières mesures de normalisation gelée — 12 octobre 2026

L'[exécution R 38255238228](https://github.com/rberrah/rberrah.github.io/actions/runs/38255238228) a réussi les **3 scénarios × 8 jeux indépendants × 3 moteurs R**, en comparant **trois états de traitement par jeu**. Les CSV par jeu, les corrections d'échelle et les versions R sont archivés dans l'artefact `bioconductor-rnaseq-frozen-normalization-audit`.

Différences de **puissance en points de pourcentage** dans le scénario à signaux rares (`sparse_signals_20v20`) :

| Méthode | Filtrage strict et réajustement des modèles aux facteurs gelés (B − A) | Changement dû à la normalisation recalculée sous filtre strict, modèles réajustés (C − B) | Différence totale (C − A) |
| --- | ---: | ---: | ---: |
| DESeq2 | −6,25 | −6,25 | −12,50 |
| edgeR QL | −12,50 | 0 | −12,50 |
| limma-voom | −3,13 | 0 | −3,13 |

Ces résultats montrent que, **pour ces données simulées**, le changement de normalisation contribue à la variation de puissance de DESeq2 après application du filtre strict, alors que les différences observées pour edgeR/voom sont concentrées dans la première étape du chemin. Les différences sont exactement additives par construction. Le FDR empirique est également exporté pour chaque bras et chaque étude ; le rapport ne prétend pas identifier une contribution universelle au FDR nominal.

**Aucune attribution généralisable ou hiérarchisation des méthodes** n'est possible avec huit cohortes artificielles indépendantes par scénario. Le découpage dépend de l'ordre des opérations, les dispersions et les autres paramètres internes sont réajustés, et la matrice est artificiellement enrichie en faibles comptages. Un véritable audit de robustesse requiert une simulation beaucoup plus vaste et des cohortes publiques indépendantes, avec les mêmes critères fixés avant inspection.

### Cinquième campagne : normalisation gelée pour analyser les variations de FDR et de puissance

La sensibilité au préfiltrage met en évidence des changements combinés de **la famille des gènes testés**, de **la normalisation** et des **modèles statistiques**. Il ne serait pas scientifiquement correct d'attribuer l'intégralité de ces changements au seul filtre.

Le nouveau protocole `multiomics-engine/rnaseq_normalization_lock.R`, avec exécution de référence `multiomics-engine/tests/bioconductor_rnaseq_normalization_lock.R`, compare trois bras sur **exactement les mêmes cohortes, échantillons, lots techniques et vérités biologiques plantées** :

| Bras | Préfiltrage externe | Facteurs de normalisation des échantillons | Modèles de référence |
| --- | --- | --- | --- |
| `reference_native` | Référence, pas de filtre additionnel | Estimés par la méthode R | DESeq2 / edgeR QL / limma-voom réajustés |
| `strict_frozen` | Filtre plus strict, indépendant des labels | **Gelés aux facteurs du bras référence** | Méthodes R réajustées, paramètres internes/dispersion réestimés |
| `strict_native` | Même filtre strict | **Réestimés sur les données filtrées** | Méthodes R réajustées |

La décomposition est strictement additive pour chaque jeu et chaque méthode :

- **Changement total** = `strict_native − reference_native`.
- **Chemin « filtre + réajustement aux facteurs gelés »** = `strict_frozen − reference_native`.
- **Chemin « normalisation réestimée avec filtre strict »** = `strict_native − strict_frozen`.

Le bras aux facteurs gelés emploie `sizeFactors` DESeq2, ou des bibliothèques effectives (offsets : `lib.size × norm.factors`) pour edgeR QL/limma-voom. Les mêmes gènes retirés sont identifiés et les facteurs réellement utilisés doivent correspondre aux facteurs gelés, sous peine d'échec explicite.

**Précaution d'interprétation majeure :** ce protocole produit une **décomposition conditionnelle au chemin suivi** ; la première différence englobe le changement de famille d'hypothèses, le filtrage natif et le réajustement de dispersion, la seconde englobe les conséquences aval d'un changement de normalisation sous le même préfiltrage strict. Ce n'est donc pas une attribution causale pure ou symétrique du « poids » du filtre versus celui de la normalisation. Une campagne plus importante, des simulations réalistes indépendantes et des études externes restent nécessaires avant toute revendication de calibration nominale du FDR à 5 %.

La CI dédiée `.github/workflows/multiomics-rnaseq-normalization-lock.yml` lance une courte campagne sous PR et des études plus nombreuses lors des exécutions programmées ou manuelles ; toutes les sorties et versions R sont archivées. **Les premières mesures sont consignées ci-dessus avec leur provenance**, mais huit jeux ne suffisent pas à établir une calibration nominale.

### Quatrième campagne : sensibilité au filtrage indépendant et à la normalisation

Dans les analyses RNA-seq, retirer les gènes faiblement exprimés modifie la **famille d'hypothèses corrigées**, et peut changer l'estimation des facteurs de normalisation, des dispersions et de la puissance. Comparer uniquement des listes de gènes significatifs sans expliciter les filtres peut donner une impression trompeuse de reproductibilité. Le manuel officiel [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) rappelle par ailleurs que son `results()` effectue par défaut un filtrage indépendant interne, et que les q-values écartées peuvent être `NA`.

Le protocole `multiomics-engine/rnaseq_filter_sensitivity.R` et `multiomics-engine/tests/bioconductor_rnaseq_filter_sensitivity.R` utilise **les mêmes cohortes simulées** et les mêmes vrais gènes différentiels dans deux bras, avec graines indépendantes d'une cohorte à l'autre :

1. **Référence :** données de comptage entières puis procédures internes normales de DESeq2, edgeR QL ou limma-voom.
2. **Filtre de sensibilité préspécifié :** conserver seulement les gènes avec au moins 20 comptages dans au moins 40 % des échantillons, puis exécuter les mêmes procédures internes. Le filtre utilise uniquement les comptages (pas l'appartenance aux groupes, pas les p-values, pas les q-values).

Pour rendre les limites du filtrage observables, les jeux synthétiques subissent au préalable un **amincissement binomial indépendant des groupes** et gelé par gène (30 % à 0,025×, 30 % à 0,2× et 40 % sans amincissement). Cette manipulation fait varier les abondances et les probabilités de détection sans inventer a posteriori de différences biologiques. Trois contextes difficiles sont réutilisés : petite cohorte (6 + 6), groupe déséquilibré avec batch (12 + 28), et signaux rares (1 %).

Les métriques restent basées sur la **vérité biologique plantée avant tout filtrage** : gènes éliminés, vrais effets éliminés, faux positifs, FDP/FDR, puissance calculée sur tous les vrais gènes, erreur log2FC. Chaque modèle recalcule **sa propre normalisation** après changement de la matrice ; les variations des facteurs d'échelle effectifs sont consignées, ainsi que les différences *appariées* de FDR et de puissance et leurs intervalles bootstrap calculés par cohorte indépendante.

**Limites et règle de décision :** le contraste « filtre plus strict vs procédure de référence » mesure l'effet **global du changement de préfiltrage et du réajustement normalisation/modèle**, pas l'effet causal isolé d'un seul paramètre. Ni un meilleur FDR observé, ni un gain de puissance ne doivent déclencher un choix automatique du filtre après lecture des résultats. Un test logiciel réussi ne vaut pas contrôle démontré du FDR nominal ni qualification confirmatoire. Le workflow `.github/workflows/multiomics-rnaseq-filter-sensitivity.yml` exécute une courte campagne en PR et prévoit des campagnes indépendantes plus grandes, avec résultats négatifs et versions R archivés.

### Troisième campagne : faible effet, petite cohorte, batch et signaux rares

Les expériences précédentes détectaient tous les vrais signaux, mais avec des effets artificiels **±1,5 log2FC** et généralement 40 sujets ; ce résultat ne permettait pas d'évaluer la puissance en conditions difficiles.

Le plan de simulation reproductible `multiomics-engine/rnaseq_stress_designs.R` prescrit quatre cas (400 gènes, exactement une étude simulée indépendante par graine, même étude évaluée par les trois méthodes R) :

| Scénario | Taille / proportion de vrais signaux | Effets artificiels | Risque mis à l'épreuve |
| --- | --- | --- | --- |
| `small_effect_20v20` | 20 + 20 sujets, 5 % | ±0,35 log2FC | Faible effet, faible puissance |
| `small_cohort_6v6` | 6 + 6 sujets, 5 % | ±0,75 log2FC | Manque de précision et de puissance |
| `unequal_batch_12v28` | 12 + 28 sujets, 5 % | ±0,75 log2FC | Lot technique corrélé à la condition, mais contraste identifiable |
| `sparse_signals_20v20` | 20 + 20 sujets, 1 % | ±0,75 log2FC | Signaux rares, nombreuses hypothèses réellement nulles |

La génération inclut dispersion négative binomiale hétérogène, facteurs de profondeur de séquençage et effets de lot variant selon les gènes. Les **sujets biologiques sont indépendants**, les groupes/batchs sont explicitement vérifiés et toutes les estimations sont comparées aux effets plantés connus.

Les scripts `multiomics-engine/tests/bioconductor_rnaseq_stress_grid.R` et `.github/workflows/multiomics-rnaseq-stress-grid.yml` distinguent :
- **PR :** 8 jeux par scénario, uniquement pour déceler les régressions logicielles.
- **Campagne mensuelle :** 120 jeux par scénario, avec possibilité de 60/120/200 exécutions sur demande ; 4 scénarios × 3 méthodes, avec archives de toutes les mesures par jeu, diagnostics, versions R, FDR moyenne, intervalle bootstrap et erreur Monte-Carlo.
- **Interprétation :** si la borne inférieure de l'intervalle de FDR dépasse 5 %, la cellule est **signalée pour revue**. L'absence d'alerte, ou le succès de CI, n'est **jamais** un certificat de contrôle du FDR, de puissance suffisante, de stabilité hors simulation ou de validité clinique.

### Résultats de la première campagne de stress R — vérifiés

Le [workflow GitHub Actions 38177062789](https://github.com/rberrah/rberrah.github.io/actions/runs/38177062789), exécuté le **11 octobre 2026**, a terminé les quatre scénarios avec **8 jeux indépendants par scénario**, chacun traité par DESeq2, edgeR QL et limma-voom (contrôle de logiciel, **pas** une campagne de calibration nominale). Les métriques et l'environnement R sont présents dans l'artefact `bioconductor-rnaseq-stress-calibration`.

**Puissance moyenne estimée** = proportion de tous les gènes réellement modifiés qui ont été retenus à BH q ≤ 0,05 :

| Scénario (8 jeux) | DESeq2 | edgeR QL | limma-voom |
| --- | ---: | ---: | ---: |
| Faible effet 20 + 20 | 2,50 % | 1,88 % | 0 % |
| Petite cohorte 6 + 6 | 2,50 % | 1,25 % | 0 % |
| Déséquilibré 12 + 28 | 25,00 % | 21,88 % | 10,63 % |
| Signaux rares 20 + 20 | 34,38 % | 28,13 % | 28,13 % |

**FDR empirique moyen** = moyenne du nombre de faux positifs divisé par le nombre total de découvertes, avec FDP = 0 quand aucun gène n'est sélectionné :

| Scénario (8 jeux) | DESeq2 | edgeR QL | limma-voom |
| --- | ---: | ---: | ---: |
| Faible effet 20 + 20 | 6,25 % | 12,50 % | 25,00 % |
| Petite cohorte 6 + 6 | 20,83 % | 0 % | 0 % |
| Déséquilibré 12 + 28 | 5,40 % | 4,17 % | 3,13 % |
| Signaux rares 20 + 20 | 13,75 % | 6,25 % | 10,42 % |

**Interprétation prudente :** les faibles puissances sont cohérentes avec des signaux difficiles et/ou un faible effectif, mais ces estimations dépendent de la distribution simulée. Les FDR sont particulièrement volatils lorsque presque aucun gène n'est détecté ; certaines cellules à FDR nul sont aussi celles où la puissance est nulle. Sur **8** études simulées, les intervalles bootstrap et les erreurs Monte-Carlo sont trop larges pour démontrer une calibration à 5 % ou classer les méthodes. Un pipeline exécuté sans erreur peut produire une analyse **non informative par manque de puissance**. Le nombre minimal de trois sujets par groupe requis par le garde-fou d'entrée ne doit jamais être interprété comme une validation de la puissance d'une vraie étude.

**Important :** cette section associe un protocole implémenté à une exécution réellement vérifiée, sans extrapoler au patient, à des études omiques externes ou au contrôle nominal du FDR. Une couverture exhaustive de la variabilité biologique, de la dépendance entre gènes, des compositions cellulaires, des tailles d'effets réalistes et des conditions de collecte nécessitera des plans complémentaires. Le critère de puissance est calculé sur **tous les vrais gènes plantés**, y compris ceux éliminés lors du filtrage.

### Deuxième campagne : mélanges contrôlés de vrais et de faux signaux

Le fichier `multiomics-engine/tests/bioconductor_mixed_truth_reference.R` génère des données de comptages négatifs binomiaux à **400 gènes**, **40 sujets indépendants** et lot technique mesuré. Trois scénarios préspécifiés introduisent **4 (1 %), 20 (5 %) ou 80 (20 %) vrais gènes différentiels**. La moitié de ces gènes est sur-exprimée et l'autre moitié sous-exprimée, à effet nominal ±1,5 sur l'échelle log2. Les mêmes simulations servent aux trois moteurs réels **DESeq2, edgeR QL et limma-voom**.

Les métriques calculées sans connaître le résultat à l'avance sont : **FDP** = faux gènes sélectionnés / tous les gènes sélectionnés dans chaque jeu, **FDR estimée** = moyenne des FDP entre jeux indépendants, **puissance** = vrais signaux détectés / tous les vrais signaux plantés (y compris ceux éliminés par filtrage), erreurs de signe, biais et erreur absolue du log2FC. Des intervalles bootstrap à 95 % sur la FDR sont produits en rééchantillonnant les **jeux indépendants**, jamais les gènes comme s'il s'agissait de nouvelles études. Les seules couvertures d'IC calculées dans cette campagne concernent les intervalles de Wald approximatifs, **non rétrécis**, de DESeq2 ; elles ne doivent pas être assimilées aux intervalles des autres moteurs.

**Ce que ce test ne démontre pas :** la présence de vrais gènes modifiés ne garantit pas la calibration. Avec 12 simulations par scénario, les estimations restent fragiles ; la vérité simulée n'englobe pas les vrais batchs biologiques, les biais d'annotation, les études longitudinales ni l'indépendance des prélèvements d'une cohorte réelle. Les fichiers CSV et `sessionInfo()` sont archivés et les échecs grossiers restent bloquants. Une calibration nominale exploitable nécessitera une simulation beaucoup plus importante, indépendante de la mise au point du pipeline, ainsi qu'une réplication externe avec protocole gelé.

### Que mesure réellement ce nouvel audit ?

Le benchmark `multiomics-engine/tests/bioconductor_null_calibration.R` exécute **trois méthodes R réelles** (DESeq2, edgeR QL, limma-voom) sur trois plans simulés, chacun avec un **effet biologique nul connu** :

- **Standard :** deux groupes de 20 sujets, séries techniques équilibrées et correction du batch explicite (12 jeux indépendants).
- **Qualité complexe :** deux groupes de 20 sujets, fortes différences de profondeur et forte surdispersion négative binomiale (8 jeux indépendants).
- **Plan déséquilibré :** 12 témoins et 28 traités, série technique partiellement liée au groupe ; rang de la matrice de design vérifié (8 jeux indépendants).

On rapporte **par méthode ET scénario** la proportion de simulations où au moins un gène nul présente q BH ≤ 0,05, avec un **intervalle binomial exact à 95 %** calculé par `binom.test`. Sous nul global, cette proportion estime le FDR et le FWER. Ces intervalles quantifient seulement l'erreur Monte-Carlo. **Zéro simulation positive sur douze donne malgré tout une borne supérieure exacte proche de 26,5 % ; sur huit, proche de 36,9 %.** Le seuil d'alerte automatisé de 35 % ne détecte que des régressions grossières.

Les sorties réelles (`bioconductor_null_replicates.csv`, `bioconductor_null_summary.csv`, `bioconductor_null_scenarios.csv`, environnement R) sont archivées dans les artefacts de la CI. L'ajout du scénario dans le code **n'est pas une exécution probante** : les valeurs mesurées doivent être lues dans un workflow réussi, jamais inventées dans cette documentation. Les résultats ne sont pas assimilables à une calibration sous alternatives non nulles.

**Critères de validation statistique :** définir le niveau alpha, le
nombre de répétitions et la précision Monte-Carlo, les scénarios, les
tests de puissance, les tolérances de calibration et les échecs
autorisés **avant de lire le benchmark**. Les bornes grossières utilisées
par la CI historique ne constituent pas une certification nominale
du contrôle de FDR à 5 %.

## Autres preuves nécessaires avant un parcours confirmatoire complet

1. Validation du protocole analytique et de la provenance (lectures,
   mapping, comptage, annotation, QC, contaminations, batch et
   indépendance réelle des sujets).
2. Plan statistique préenregistré, taille d'effet pertinente et
   justification de la puissance ; procédures d'exclusion et
   d'analyses de sensibilité préspécifiées.
3. Vérifications de diagnostics/dispersion et confrontation à une méthode
   de référence indépendante ; non-concordance explicitement conservée.
4. Réplication externe de résultats biologiques sur des données non
   utilisées pour améliorer ou sélectionner la méthode.
5. Révision humaine indépendante des diagnostics et des conclusions.
   La validation de l'interface par utilisateurs fait l'objet d'un
   protocole distinct.

**Règle de communication :** même lorsque les contrôles R1 passent,
écrire « exécution de la méthode RNA-seq de référence et préflight
technique satisfaisants, revue méthodologique requise », et non
« résultat confirmé », « validé universellement » ou « prêt pour
décision clinique ».
