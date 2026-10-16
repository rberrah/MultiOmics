Pipeline omique ou multi-omique PMx Explain — guide simple des fichiers et méthodes
=======================================================================

FRANÇAIS
========

LECTURE RAPIDE — LES MOTS UTILISÉS DANS L’OUTIL
=================================================

L’interface affiche d’abord le sens pratique. Le terme technique reste disponible entre parenthèses ou dans l’aide « ? » lorsque cela est utile pour la reproductibilité.

- Tableau des échantillons = metadata / sample sheet. Il indique qui a été étudié, quel prélèvement a été analysé, dans quel groupe et à quel moment.
- Matrice = tableau contenant les valeurs mesurées. Les lignes sont des gènes, protéines ou métabolites ; les colonnes correspondent aux mesures des échantillons.
- Critère étudié = outcome / endpoint. Exemple : réponse au traitement, concentration, score, décès ou délai jusqu’à un événement.
- Série technique = batch. Exemple : plaque, jour, run ou lot technique.
- Facteur d’ajustement = covariable. Exemple : âge, sexe, centre ou traitement concomitant.
- Variable biologique = feature. Selon la couche, il s’agit d’un gène, d’une protéine ou d’un métabolite.
- Résultat corrigé = q BH / FDR. Il tient compte du grand nombre de variables testées simultanément.
- Validation croisée = évaluation d’une prédiction sur des sujets qui n’ont pas servi à construire le modèle.

PRINCIPE GÉNÉRAL
================

Une seule matrice RNA, protéomique ou métabolomique suffit pour une exploration, une comparaison de groupes, une analyse longitudinale ou une association à un critère lorsque les données et le plan le permettent. Deux couches ou davantage sont nécessaires aux méthodes d'intégration inter-omique (MOFA2, DIABLO, corrélations croisées). Le tableau des échantillons peut contenir d'autres omiques non importées : seules les omiques effectivement chargées comptent dans les statistiques.

Pour explorer une seule matrice sans tableau d'échantillons, il est possible de générer un tableau minimal APRÈS confirmation explicite que chaque colonne représente un sujet biologiquement indépendant avec une seule mesure. L'outil ne devine ni les groupes, ni les temps, ni les réplicats. Cette facilité est descriptive seulement et ne permet pas de prétendre à une inférence sur un design inconnu.

Le logiciel automatise ce qui peut l’être sans masquer les hypothèses. Il ne transforme pas un plan expérimental non identifiable en analyse valide, ne déduit pas la puissance du seul nombre de sujets, ne fusionne pas artificiellement les q-values de plusieurs méthodes et ne présente pas une validation croisée interne comme une validation externe.

Le parcours conseillé est :
1. poser la question biologique ;
2. décrire correctement sujets, prélèvements et mesures ;
3. vérifier le plan et la qualité ;
4. exécuter la méthode compatible ;
5. lire amplitude + incertitude + q/FDR ;
6. si plusieurs omiques sont chargées : examiner leur convergence sans inventer de corrélations en cas d'omique unique ;
7. replacer le signal dans son contexte biologique ;
8. conserver les limites d’interprétation et le manifeste reproductible.

1. TABLEAU DES ÉCHANTILLONS
===========================

Le tableau des échantillons ne contient pas les milliers de valeurs moléculaires. Il explique simplement à l’outil à quoi correspondent les colonnes des matrices.

Colonnes essentielles :
- subject_id : participant, animal, culture ou unité biologique indépendante ;
- sample_id : prélèvement biologique ;
- assay_id : mesure technique correspondant exactement à une colonne de matrice ;
- omic : transcriptomics, proteomics ou metabolomics.

Colonnes utilisées selon l’étude :
- condition : groupe ou traitement ;
- timepoint : visite ou moment expérimental ;
- batch : série technique ;
- technical_replicate : répétition technique ;
- outcome : critère étudié ;
- survival_time : durée de suivi ;
- survival_event : événement observé/censuré ;
- sample_type : biological / pooled_qc / qc / blank pour les workflows MS ;
- injection_order : ordre d’injection MS ;
- toute autre colonne peut être sélectionnée explicitement comme facteur d’ajustement.

Exemple :
P001 / P001_T0 / RNA001  / transcriptomics
P001 / P001_T0 / PROT001 / proteomics
P001 / P001_T0 / MET001  / metabolomics

P001 est la même personne ; P001_T0 est le même prélèvement ; RNA001, PROT001 et MET001 sont trois mesures différentes. Un réplicat technique n’est jamais compté comme un nouveau sujet biologique.

VÉRIFICATION AUTOMATIQUE DES IDENTIFIANTS ET CRITÈRES
- Le même subject_id peut apparaître plusieurs fois (plusieurs omiques, visites ou répétitions), mais un même participant ne doit pas recevoir deux valeurs contradictoires du critère clinique (outcome).
- Pour une analyse de survie, renseigner survival_time et survival_event ensemble. L’événement doit être codé 1 (observé) ou 0 (censuré), avec une durée strictement positive. Ces valeurs doivent rester identiques entre les omiques d’un même sujet.
- Pour un comptage (outcome de type count), fournir un entier positif ou nul.
- Dans une comparaison de groupes indépendants ou une étude longitudinale avec bras fixes, un sujet appartient à un seul groupe condition ; une différence de groupe entre lignes du même sujet est refusée.
- En cas d'incohérence, le calcul est interrompu au lieu de choisir silencieusement la dernière ligne importée. Corriger le tableau des échantillons, pas les valeurs des matrices.
- Ces contrôles n'établissent pas que les annotations sont vraies ; ils préviennent uniquement certaines erreurs de saisie et d'appariement.

2. CONTRÔLE QUALITÉ ET FAISABILITÉ
===================================

Avant une interprétation forte, l’outil vérifie notamment :
- nombre de sujets biologiquement indépendants dans chaque groupe ;
- déséquilibre entre groupes ;
- chevauchement réel des sujets entre omiques ;
- complétude des mesures répétées ;
- nombre de sujets disposant réellement de plusieurs temps ;
- données manquantes ;
- mesures atypiques ;
- accord entre répétitions techniques ;
- série technique et risque de confusion avec groupe/temps ;
- pour la survie, nombre de sujets exploitables et nombre d’événements observés ;
- éligibilité des méthodes avancées MOFA2 et DIABLO.
- pour les plans longitudinaux à deux groupes : interaction groupe × temps avec modèle à intercept aléatoire ;
- pour les plans longitudinaux à trois groupes ou plus : test omnibus multi-ddl de l’interaction groupe × temps, puis classement corrigé BH ; les comparaisons post-hoc spécifiques ne sont pas inventées automatiquement.

`preAnalysisDiagnostics` possède quatre niveaux :
- ready : aucun problème déterministe détecté ;
- usable_with_cautions : analyse possible avec précautions documentées ;
- review_required : une limite de qualité/effectif/design doit être explicitement revue ;
- blocked_or_redesign_required : au moins un problème déterministe empêche une interprétation forte défendable ; il faut corriger le plan/l’annotation ou limiter l’analyse à une description appropriée.

Exemples de situations considérées comme bloquantes ou non défendables pour une inférence forte :
- comparaison de groupes avec moins de deux sujets indépendants dans un groupe ;
- série technique complètement confondue avec la condition biologique ;
- route longitudinale sans véritables mesures répétées chez plusieurs sujets ;
- survie sans suffisamment de sujets exploitables/événements pour estimer le modèle demandé ;
- crossover sans période et séquence explicitement modélisées.

Un groupe contenant moins de cinq sujets indépendants est signalé comme hautement exploratoire, même si un logiciel peut techniquement produire une p-value.

La puissance statistique n’est volontairement PAS déduite du seul nombre de sujets. Une estimation défendable de puissance nécessite un effet attendu, une variabilité et le modèle prévu ; une simulation est souvent préférable pour les designs multi-omiques complexes.

3. LC-MS / GC-MS
=================

ENTRÉE BRUTE STANDARDISÉE
-------------------------
Pour les données non ciblées, la voie brute de référence accepte mzML / mzXML et autres formats ouverts lisibles par xcms. Les formats propriétaires Thermo / Waters / Agilent / Bruker doivent être convertis en mzML en amont, typiquement avec ProteoWizard/msconvert lorsque le lecteur vendor correspondant est disponible.

Pipeline local de référence :
1. manifeste explicite raw_file ↔ sample_id / assay_id / condition / sample_type / injection_order ;
2. détection des pics chromatographiques par centWave ;
3. alignement des temps de rétention par obiwarp lorsque plusieurs injections sont présentes ;
4. groupement/correspondance des pics par peak density ;
5. gap filling par ChromPeakAreaParam ;
6. export d’une matrice `metabolomics_peak_area.csv`, des définitions de features et d’un manifeste de traitement.

Commande :
Rscript multiomics-engine/run_raw_ms.R raw_ms_manifest.csv raw_ms_output raw_ms_parameters.json

Templates :
- /multiomics/templates/raw_ms_manifest.csv
- /multiomics/templates/raw_ms_parameters.json

Les fichiers CSV/TSV MS de type long ou large séparent les aires des pics des
colonnes techniques m/z, temps de rétention, adduits, scores et annotations
présentes dans l'export. Le fichier d'origine reste inchangé, son SHA-256 et
le SHA-256 de la matrice dérivée sont inscrits au rapport quand Web Crypto est
disponible. La table technique séparée est téléchargeable. Un identifiant
ChEBI donné par l'utilisateur reste une déclaration et non une identification
validée chimiquement. Plusieurs pics d'un même ChEBI ne sont jamais additionnés
ni choisis arbitrairement.

Important : cette étape crée des features analytiques m/z–RT. Elle ne transforme pas automatiquement ces features en identités métabolites. L’annotation ChEBI/HMDB/KEGG/PubChem/InChIKey reste une étape distincte, conservatrice et traçable.

Lorsque sample_type et injection_order sont fournis :
- blanks et pooled-QC sont exclus des échantillons biologiques ;
- un ratio blank/biologique est par défaut un SIGNAL DE QC et ne supprime pas automatiquement une feature ;
- un mode de suppression explicite reste disponible pour un SOP de laboratoire validé ;
- une dérive instrumentale peut être corrigée à partir des pooled-QC ordonnés ;
- la stabilité des pooled-QC peut être évaluée par RSD ;
- aucune imputation MNAR n’est effectuée par défaut ;
- une hypothèse left-censored doit être déclarée et accompagnée d’une analyse de sensibilité sans imputation.

4. DIFFÉRENCES TECHNIQUES ET FACTEURS D’AJUSTEMENT
===================================================

Une série technique totalement confondue avec le groupe ou le temps ne peut pas être « corrigée » statistiquement : les deux effets sont non séparables. L’outil bloque alors l’interprétation forte.

Quand le design le permet, les séries techniques et facteurs d’ajustement choisis sont inclus dans le modèle. Aucun facteur d’ajustement n’est inventé automatiquement.

5. QUESTIONS ET ROUTES STATISTIQUES
====================================

Explorer :
- une seule omique : ACP descriptive sur les variables de cette omique, avec au moins trois sujets indépendants et sans prétendre mesurer un couplage inter-omique ;
- plusieurs omiques : ACP multi-blocs équilibrée et prise en compte documentée des blocs partiellement manquants ;
- backend de référence : MOFA2 seulement pour au moins deux couches et si les garde-fous d’éligibilité sont compatibles.
- une analyse exploratoire ne prouve pas une différence entre groupes et les loadings ne sont pas des tailles d'effet de groupe.

Comparer des groupes :
- estimation de l’amplitude et du sens ;
- tests adaptés à la nature de la couche ;
- correction Benjamini-Hochberg ;
- comparaisons inter-omiques sur le chevauchement réel des sujets.

Temps / mesures répétées :
- le sujet reste l’unité de corrélation ;
- lmerTest dans le backend de référence ;
- pente aléatoire du temps essayée automatiquement seulement si le plan la supporte, sinon repli sur intercept aléatoire ;
- interaction condition × temps lorsque le plan l’autorise.

Critère clinique ou expérimental :
- continu : régression linéaire ;
- oui/non : logistique ;
- comptage : Poisson ;
- multicatégoriel : modèle adapté ;
- survie : Cox ;
- la performance prédictive est évaluée séparément de l’association statistique.

6. RNA-SEQ ET MATRICES NORMALISÉES
===================================

Backend R de référence : une seule omique est acceptée ; les moteurs de comparaison sont choisis selon le design et le type de données, sans exiger une couche artificielle. Le site GitHub Pages n'héberge pas R et indique explicitement si le calcul R a été exécuté, indisponible ou non demandé.
- comptes RNA-seq : DESeq2 comme analyse principale et edgeR + limma-voom comme analyse de sensibilité automatisée lorsque disponibles ;
- edgeR quasi-likelihood est également implémenté comme adaptateur de référence ;
- matrices normalisées/log : limma ;
- matrices de design non pleines-rang : rejetées plutôt que fitted avec des coefficients non identifiables.

La concordance DESeq2/voom peut rapporter :
- corrélation des statistiques signées ;
- nombre de résultats significatifs par méthode ;
- chevauchement/Jaccard ;
- accord du sens de l’effet.

Les q-values des différentes méthodes ne sont jamais moyennées ou transformées en « q consensus ». Le désaccord reste visible comme analyse de sensibilité.

7. MOFA2
========

La route de référence :
- nécessite au moins 16 sujets réellement communs aux blocs analysés ;
- élimine les features constantes ou insuffisamment observées ;
- conserve une seed déterministe ;
- enregistre facteurs, poids et variance expliquée lorsqu’elle est disponible ;
- considère explicitement la présence d’effets techniques connus avant l’interprétation.

Un facteur MOFA2 est une structure de covariance latente. Ce n’est ni une preuve causale, ni automatiquement un biomarqueur, ni une validation externe.

8. DIABLO
=========

La route DIABLO est réservée à un critère catégoriel supervisé et nécessite une classe suffisamment représentée pour permettre une validation interne.

Lorsque le plan le permet :
- tuning de la sparsité keepX par `tune.block.splsda()` ;
- validation M-fold répétée ;
- BER comme critère de tuning ;
- nombre de composantes et keepX retenus enregistrés ;
- `mixOmics::perf()` pour l’évaluation interne ;
- seed enregistrée.

Une signature DIABLO est un signal multivarié supervisé. Une validation externe indépendante reste nécessaire avant toute revendication de biomarqueur généralisable.

9. ASSOCIATION, PRÉDICTION ET CAUSALITÉ
========================================

Ces trois affirmations sont séparées dans l’interface :

Association : « cette variable est associée au groupe/critère après les ajustements déclarés ».
Prédiction : « ce modèle prédit un sujet non utilisé pour sa construction avec telle performance ».
Causalité : nécessite un design et des hypothèses causales supplémentaires ; elle n’est jamais déduite d’une simple association ou d’une bonne prédiction.

La validation prédictive disponible évite les fuites : sélection, ajustements, imputation/centrage/échelle et tuning sont appris dans l’entraînement. La performance finale est calculée sur les sujets laissés de côté.

Métriques :
- binaire : AUC, accuracy, log-loss ;
- multiclasse : accuracy ;
- continu/comptage : RMSE et R² hors-échantillon ;
- survie : Cox ridge avec tuning interne et C-index de Harrell dans les folds externes.

Une validation croisée interne ne remplace jamais une validation externe.

VALIDATION EXTERNE FIGÉE
------------------------
Le dépôt contient un validateur séparé qui consomme uniquement des prédictions déjà générées par un modèle figé. Il ne sélectionne aucune variable, ne retune aucun hyperparamètre et ne réentraîne aucun modèle sur la cohorte de validation.

Types pris en charge :
- binaire : AUC, Brier, log-loss, calibration ;
- continu : RMSE, MAE, R², calibration ;
- comptage : RMSE/MAE/R² + déviance de Poisson ;
- survie : C-index de Harrell ;
- multiclasse : accuracy, balanced accuracy, log-loss et macro-AUC one-vs-rest.

Une incertitude bootstrap déterministe est calculée lorsqu’elle est estimable. Les ré-échantillonnages non estimables restent manquants ; ils ne déclenchent jamais un ré-entraînement.

Utilisation :
- depuis le site : charger le CSV de prédictions figées et le JSON de configuration dans le panneau « Validation externe depuis le navigateur » ; le calcul est envoyé au endpoint local R `/external-validation` ;
- en ligne de commande :
  Rscript multiomics-engine/run_external_validation.R external_validation_predictions_binary.csv external_validation_binary.json external_validation_output

Dans les deux voies, `independent_cohort=true` est requis pour utiliser l’étiquette « validation externe ». Aucune variable n’est sélectionnée et aucun modèle n’est réentraîné pendant cette étape.

Templates :
- /multiomics/templates/external_validation_predictions_binary.csv
- /multiomics/templates/external_validation_binary.json

Le statut `external_validation` n’est utilisé que si l’analyste affirme explicitement que la cohorte est indépendante du développement. Le logiciel peut empêcher un mauvais étiquetage, mais il ne peut pas prouver l’indépendance historique de la cohorte : celle-ci doit être documentée par la provenance de l’étude.

10. INTERPRÉTATION BIOLOGIQUE
=============================

Reactome replace les résultats dans des voies biologiques connues. L’univers de référence doit correspondre aux variables effectivement mesurées et mappables. L’outil peut également montrer combien de couches soutiennent une même voie.

Un enrichissement de voie ne prouve pas que la voie est activée, causale ou cliniquement pertinente. La couverture des identifiants et les molécules effectivement contributrices doivent rester visibles.

11. IDENTIFIANTS
================

Résolution conservatrice :
- transcriptomique : Ensembl et symboles de gènes ;
- protéomique : UniProt / Ensembl ;
- métabolomique : ChEBI et noms/alias déterministes ;
- InChIKey : UniChem vers ChEBI lorsque possible ;
- HMDB : HMDB → ChEBI via UniChem ;
- KEGG Compound : recherche exacte des cross-références ChEBI ;
- PubChem CID : PubChem → ChEBI via UniChem quand le type a été déclaré ;
- une correspondance externe n’est acceptée que si elle conduit à un ChEBI unique ;
- plusieurs correspondances = ambiguous ;
- aucune correspondance = unresolved.

Le pipeline préfère l’absence de mapping à un mapping biologiquement incertain.

12. RAPPORT REPRODUCTIBLE
=========================

La page de résultats comporte désormais une synthèse scientifique expliquant :
- ce que les résultats permettent de conclure ;
- ce qu’ils ne permettent pas de conclure ;
- les preuves q/FDR visibles ;
- les avertissements qui doivent rester dans le rapport ;
- la frontière entre association, prédiction et causalité.

Deux exports supplémentaires sont disponibles :
- `multiomics-methods-report.md` : rapport humainement lisible des choix et limites ;
- `multiomics-reproducibility-manifest.json` : manifeste machine-readable contenant question scientifique, design, critère, paramètres de l’interface, métadonnées locales des fichiers, preuves corrigées visibles, avertissements et limites d’interprétation.

Le manifeste enregistre également des empreintes SHA-256 calculées localement sur les fichiers d’entrée et sur une représentation canonique de la table d’échantillons. Elles permettent de vérifier qu’un rapport correspond exactement aux données importées. Un fingerprint rapide éventuellement conservé dans l’interface n’est pas considéré comme un substitut cryptographique au SHA-256.

Le manifeste n’embarque pas les données de recherche elles-mêmes. Pour une analyse destinée à un manuscrit, conserver séparément :
- fichiers d’entrée immuables ;
- checksums ;
- manifeste ;
- sorties complètes ;
- versions des packages/session R ;
- code/commit utilisé.

13. COMMENT LIRE LES PRINCIPAUX RÉSULTATS
==========================================

- Amplitude / fold ratio : taille et sens du changement.
- p : test avant correction multiple.
- q BH / FDR : résultat corrigé pour le grand nombre de variables testées.
- Pattern : forme générale du changement observé.
- r référence / r comparaison : corrélation dans chaque groupe.
- Δr : différence de corrélation entre groupes.
- FDR de voie : enrichissement corrigé dans l’univers mesuré/mappable.
- Nombre de couches concordantes : nombre d’omics soutenant indépendamment la même voie.

Une petite p-value seule n’est jamais présentée comme une conclusion biologique suffisante.

14. VALIDATION PUBLIQUE ET LIMITES
===================================

La suite automatisée comprend notamment des jeux publics tels que Nutrimouse, TCGA breast, IntLIM NCI-60/BRCA, AgingHFCD, STATegra, LRRK2 G2019S, PaintOmics avec signal multi-omique planté et missRows NCI-60.

L’automatisation peut détecter de nombreuses incompatibilités et réduire les erreurs statistiques courantes. Elle ne peut pas certifier automatiquement :
- absence de biais de sélection ;
- absence de confondeurs non mesurés ;
- pertinence clinique d’une taille d’effet ;
- puissance suffisante sans effet/variance attendus ;
- causalité sans design causal ;
- validité externe sans nouvelle population indépendante.

Ces limites doivent rester visibles même lorsqu’une analyse s’exécute sans erreur.

BACKEND R
=========

Installation :
Rscript multiomics-engine/install_backend_dependencies.R

Démarrage local :
Rscript multiomics-engine/run_backend.R

Adresse par défaut :
http://127.0.0.1:8787

Un serveur distant recevant des données de recherche doit ajouter au minimum TLS, authentification, limites de taille/requêtes et allow-list d’origine.

15. INTERFACE GUIDÉE ET LECTURE DES FIGURES
===========================================

Le bouton « Voir une analyse complète » charge les matrices de démonstration
(3 omiques) et exécute la même analyse que pour les données personnelles.
Dans les résultats, les vues proposées dépendent des données : la heatmap fonctionne avec une seule omique ; les vues de voies et de métabolites sont informatives seulement si des effets et des identités exploitables existent. La démo peut charger RNA seul, protéines seules, métabolites seuls ou les trois couches. Les données d'exemple sont synthétiques, non une preuve biologique.
Trois familles de vues sont accessibles :

1. Heatmap : chaque ligne montre une variable centrée-réduite (z-score).
   Un bleu ou rouge ne signifie pas un log2 fold change ni une différence
   statistiquement significative. Les échantillons sont ordonnés selon
   les annotations de l'étude ; absence de clustering automatique.

2. Voies biologiques : métabologramme circulaire inspiré de la figure 4C
   de Guyon et al. (EMBO Molecular Medicine 2022 ; 14:e15343).
   À gauche figurent les métabolites, à droite les transcrits enzymatiques.
   L'utilisateur peut sélectionner glycolyse, cycle TCA, acides aminés ou
   pentoses phosphates. La sélection est une liste explicite de membres,
   PAS un test d'enrichissement ni un score d'activité. Une variable
   non reconnue ou absente ne doit pas être interprétée comme inchangée.
   Le centre indique une moyenne descriptive d'effets sélectionnés sur leur
   échelle déclarée (log2 ou unité native). Il ne s'agit pas d'un score
   d'activité et les couleurs ne sont pas comparables entre omiques.

3. Carte métabolique : représentation pédagogique, inspirée de la figure 4D
   de la même publication. Les ovales sont des métabolites ; les rectangles,
   des transcrits d'enzymes. Les effets mesurés sur une échelle exploitable
   (log2 ou échelle native déclarée) sont colorés avec une étiquette d'échelle.
   Les données absentes, les identités inconnues, les signaux MS multiples
   non agrégés et les familles discordantes restent séparés ou grisés.
   Le réseau cible les régions documentées et leurs voisins, mais les liens
   sont schématiques, sans flux, sans stœchiométrie Rhea vérifiée ni causalité.

IMPORTANT : l'article utilise des traceurs isotopiques au carbone 13.
   Le logiciel ne calcule PAS de flux 13C à partir de simples
   abondances métabolomiques et mesures d'ARN. Une interprétation
   de flux exige les isotopologues, leur correction et un modèle de
   traçage spécifique.

La mise en page n'altère aucun calcul, seuil q-BH ni modèle R/JS.
Les paramètres détaillés et workflows experts restent accessibles
à la demande ; la démo sert à comprendre, non à valider une méthode.
