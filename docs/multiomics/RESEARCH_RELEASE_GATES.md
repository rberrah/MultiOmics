# MultiOmics — critères de sortie pour un logiciel de recherche (13 octobre 2026)

## Research v1 freeze status — 19 October 2026

The Research v1 freeze decision is defined by [RESEARCH_V1_FREEZE_20261019.md](RESEARCH_V1_FREEZE_20261019.md) and `research-v1-manifest.json`. The frozen baseline includes the real external NKI survival and GSE191328 longitudinal qualifications, plus the exact 171-package R reference closure. This freeze preserves unfavorable evidence and does **not** convert issues #44, #47 or #20 into resolved scientific claims.


**Objectif** : produire une analyse mono- ou multi-omique exploitable, transparente et reproductible pour un article scientifique *du chercheur*. Ce document n'annonce ni une certification universelle des analyses, ni un article sur l'outil.

## Décision de diffusion

**Ne pas labelliser « confirmatoire validé » sur la base des seuls tests internes ou du succès CI.** Autoriser une analyse confirmatoire seulement pour les combinaisons de données, plan, modèle et moteur dont les critères ci-dessous sont documentés. Dans les autres cas : exploration explicitement marquée, ou blocage avec explication courte.

## État d'intégration vérifié — 14 octobre 2026

Les changements MultiOmics suivants ont été fusionnés dans `main` après des CI de branche réussies au premier passage :

- #73 : environnement R/Bioconductor figé et rejeu indépendant ; #77 : acceptation navigateur, FDR/puissance et jeux TCGA ; #79 : ZIP du navigateur ; #89 : import MS peak-AUC avec gestion d'état explicite.
- #90 : intention confirmatoire et contraste RNA-seq choisi par le chercheur, transmis à R et contrôlé dans les résultats.
- #91 : export intégral des variables ajustées DESeq2/voom/limma ; limites de transfert et archivage privé local explicites.
- #92 : séparation des tables exploratoires et des CSV R exhaustifs, vérification d'intégrité et ZIP.
- #93 : RNA-seq apparié avec blocage sujet, matrice/metadata concordants, erreurs bloquantes et tests réels de modèles R.
- #94 : aide débutant, template 4 paires/8 échantillons, obligation de R apparié, sans second sélecteur de contraste contradictoire.
- #95 : vraie session Chromium → R/DESeq2 → CSV R complet → ZIP → vérificateur Python ; cohorte publique Airway GSE52778 et confrontation au DESeq2 natif. Tests sur le code assemblé, non sur un ancien backend épinglé.

L'ancienne liste de dépendances #54/#81/#82/#83/#84/#86/#87 correspond aux branches de développement, **pas** à leur intégration effective. Les validations de branche ne remplacent pas la recette du dernier commit déployé. Suivre le SHA de `main`, ses workflows de déploiement, puis le smoke test public. **Ne jamais convertir des CI vertes en certification statistique ou clinique.**

## Décomposition des filtres Cook et BH — expérience native appariée (18 octobre 2026)

La [PR #57](https://github.com/rberrah/MultiOmics/pull/57)
a ajusté **192 véritables modèles DESeq2 sous des simulations nulles**
et comparé **768 résultats BH** obtenus par les quatre combinaisons
filtre Cook actif/inactif × filtrage indépendant actif/inactif.
L'authenticité du calcul natif, l'identité des cohortes dans chaque
groupe de quatre appels, les p-values et le recalcul BH ont été
contrôlés indépendamment.

**Aucun des six scénarios (32 études par scénario) n'a changé de
nombre d'études avec une fausse découverte selon ces quatre bras.**
L'independent filtering n'a écarté aucun gène supplémentaire
évaluable, et Cook a retiré seulement deux p-values dans toute
la campagne (aucun faux appel au niveau étude). Le remplacement
des valeurs aberrantes était désactivé **partout** afin d'isoler
le filtrage : cette expérience ne reproduit donc pas la totalité
du comportement DESeq2 par défaut.

Les [effectifs réels, le protocole, la provenance des 768 tables et
les limites d'inférence](DESEQ2_COOK_FILTER_NATIVE_FINDINGS_20261018.md)
sont conservés dans Git. **Il serait incorrect de présenter ce
résultat comme l'inutilité générale des filtres Cook/BH sur des
données biologiques réelles**. L'issue #44 reste ouverte, et aucun
comportement de production n'est changé.

## Contre-épreuve du générateur Gamma–Poisson (18 octobre 2026)

La [PR #56](https://github.com/rberrah/MultiOmics/pull/56)
a rejoué **384 ajustements DESeq2 natifs sur 96 nouvelles études
entières**, construites par un **nouveau code Gamma–Poisson**,
avec dispersions par gène substantiellement **hors de la tendance
moyenne–dispersion**. Deux programmes scientifiques distincts ont
produit des preuves numériques identiques sur **384/384 lignes**,
16 mesures communes, avec **zéro divergence**.

Les taux d'études ayant au moins un faux appel BH sont
**15/32, 25/32 et 25/32** pour le bras `dispFit`
contre **2/32, 3/32 et 5/32** pour le MAP final
dans les trois scénarios. Ce résultat **renverse l'avantage**
de la tendance trouvé avec un générateur où la dispersion vraie
suivait presque parfaitement la moyenne. Voir
[les chiffres et limites de cette falsification](INDEPENDENT_GAMMA_POISSON_NATIVE_FINDINGS_20261018.md).

**Conséquence Research v1 :** une méthode ne peut être choisie,
déclassée ou présentée comme FDR-calibrée sur la seule base de
simulations dont les hypothèses favorisent sa structure. Le nouveau
générateur reste artificiel, et un contrôle FDR à 5 % sur des
données RNA-seq biologiques n'est **pas** établi.
L'issue #44 reste ouverte et aucun paramètre de production
n'a changé.

## Diagnostic statistique complémentaire — dispersion et test de Wald (18 octobre 2026)

La [PR #53](https://github.com/rberrah/MultiOmics/pull/53)
a exécuté **288 vrais ajustements DESeq2** à partir de **96 jeux nuls
synthétiques indépendants**, en comparant les résultats avec
dispersion estimée et avec **vraie dispersion du générateur NB**
sur les mêmes matrices. Les taux observés de jeux présentant au
moins un faux appel sont passés de **4/32 à 1/32** (plan équilibré),
**5/32 à 1/32** (forte dispersion) et **3/32 à 1/32** (lot/groupes
déséquilibrés) lorsque la vraie dispersion était fournie.
Le moment moyen de la statistique de Wald au carré s'est rapproché de
1 dans les trois plans. Ces indices suggèrent un rôle des
estimations de dispersion **dans ces simulateurs**, mais ne
fournissent ni identification causale complète ni garantie de FDR
à 5 % dans une étude réelle, les IC binomiaux restant larges.

Voir les [résultats natifs, mesures de dispersion et limites](DESEQ2_DISPERSION_WALD_NATIVE_FINDINGS_20261018.md).
La **revue statistique humaine reste obligatoire** ; aucun
paramètre de production n'a été changé.

## Preuves complémentaires indépendantes — LC-MS et mécanisme FDR

La [sensibilité sur les douze chromatogrammes originaux Sacurine](REAL_EXTERNAL_SENSITIVITY_AND_NULL_ORACLE_FINDINGS_20261018.md)
(PR #50, cinq contrôles réussis) a donné **2 727/3 135 features
réciproquement retrouvées** à tolérance fixée 20 ppm/30 s sous le
couple de perturbation xcms ppm25/SNR8 → ppm20/SNR10 ; 97,36 %
de concordance du drapeau RSD QC sur **2 497** paires où les deux
valeurs étaient calculables. Ce résultat ne valide pas l'identité
chimique, ne masque pas les mesures de blanc non évaluables et
ne remplace pas une comparaison de pipelines de prétraitement
**indépendants**.

Le diagnostic [normalisation DESeq2, 96 études nulles appariées](REAL_EXTERNAL_SENSITIVITY_AND_NULL_ORACLE_FINDINGS_20261018.md)
(PR #51, vrais ajustements R) a montré que fournir les véritables
facteurs de profondeur du simulateur **ne supprime pas l'excès des
fausses découvertes**. Le mécanisme initial reste **non résolu**,
l'issue #44 est ouverte et aucun seuil statistique de production
n'a été modifié.

## Mise à jour de qualification — 18 octobre 2026 (1 728 fits R)

La campagne à haute répétition préspécifiée de [PR #43](https://github.com/rberrah/MultiOmics/pull/43) a fait tourner les **vrais** DESeq2, edgeR-QL et limma-voom sur **128 jeux nuls indépendants × 3 scénarios**, et **64 jeux à vérité mixte × 3 scénarios**, soit 1 728 ajustements R au total. Les données, graines, versions et contrôles BH indépendants sont archivés dans le [workflow scientifique réel](https://github.com/rberrah/MultiOmics/actions/runs/39143010669). Au nul global, DESeq2 a rendu au moins un faux positif dans **22/128 (17,19 %), 23/128 (17,97 %) et 18/128 (14,06 %)** jeux des trois scénarios, alors que la cible BH nominale est de 5 %. Les **bornes inférieures des IC exacts simultanés à 95 %** sur les neuf cellules méthode × scénario étaient **9,11 %, 9,70 % et 6,82 %**, toutes supérieures à 5 % : **une inflation est donc démontrée pour ces générateurs**, sans pouvoir être généralisée automatiquement aux données biologiques réelles. Les intervalles edgeR-QL et limma-voom restent inconclusifs. Les q-values BH ont passé un audit numérique indépendant ; cela ne démontre pas la calibration des p-values en amont.

**Conséquence opérationnelle :** le statut `PASS` de la CI atteste seulement l'exécution et la reproductibilité informatique. Il n'autorise ni l'étiquette « FDR universellement validé à 5 % » ni le choix automatique d'une autre méthode a posteriori. Revue statistique indépendante maintenue, analyse de la cause dans [issue #44](https://github.com/rberrah/MultiOmics/issues/44) ; voir les [résultats détaillés et leurs limites](RNASEQ_HIGH_REP_FDR_FINDINGS_20261018.md). Le moteur de production n'est pas changé par ce résultat.

## Alerte scientifique ouverte — inflation FDR sous simulations

La PR diagnostique #75 reste **ouverte et non fusionnée comme correction inférentielle**. Son exécution réelle : [GitHub Actions #38369835987](https://github.com/rberrah/rberrah.github.io/actions/runs/38369835987). Cinq nuls globaux préspécifiés × 60 réalisations, trois variantes DESeq2 par jeu = 900 fits. Pour la variante DESeq2 Wald par défaut, le nombre de jeux avec au moins une fausse découverte est respectivement **10/60, 11/60, 13/60, 7/60 et 6/60**, soit des fréquences observées de **16,7 %, 18,3 %, 21,7 %, 11,7 % et 10,0 %**, dans cet ordre de scénarios (dispersion fixée, tailles de librairies fixées, effet de lot absent, nul initial équilibré, nul initial déséquilibré).

En nul global, cet indicateur est à la fois la fréquence empirique FWER et l'estimateur FDR sur ces simulations ; sa hausse ne prouve pas que DESeq2 soit généralement mal calibré sur des données biologiques réelles. Avec seulement 60 réalisations/scénario, les intervalles de confiance sont larges. Les variantes de test évaluées ne doivent **jamais être choisies a posteriori** d'après leurs faux positifs. L'attribution du signal (génération NB, estimation de dispersion, formulation des tailles de librairie, indépendance, lot, filtrage, tests de référence) exige des contrôles supplémentaires.

Conséquence produit : l'exécution du moteur R et un contrat d'entrée validé ne signifient **pas** `FDR universellement contrôlé à 5 %` ; l'interface et les rapports doivent conserver le statut de revue indépendante requise. Ne pas masquer ce résultat négatif. Distinguer toujours validité du câblage, calibration fréquentiste, validité biologique externe et généralisation prédictive.

## Complément empirique — attribution mécanistique des faux appels (14 octobre 2026)

Deux campagnes **préspécifiées** et exécutées avec les vrais packages R/Bioconductor ont été intégrées, sans modification du moteur DESeq2 de production :

- [PR #97 — facteurs de taille estimés ou exacts du générateur](https://github.com/rberrah/rberrah.github.io/pull/97), [CI et artefact R](https://github.com/rberrah/rberrah.github.io/actions/runs/38565888024) : deux tailles de cohorte (40/80 sujets), 40 jeux nuls par design, 160 fits DESeq2. Fréquences de jeux avec au moins un faux appel : **17,5 % / 17,5 %** (estimé/oracle, n=40), **10 % / 7,5 %** (n=80). Les facteurs de taille connus ne font pas disparaître le problème.
- [PR #100 — dispersions et référence de Wald](https://github.com/rberrah/rberrah.github.io/pull/100), [CI et artefacts R](https://github.com/rberrah/rberrah.github.io/actions/runs/38602501639) : deux nouveaux échantillons de 40 jeux nuls, n=40 et n=80, 320 vrais ajustements. Tableau ci-dessous.

| Paramètres DESeq2 | N=40 : études avec faux appels/40 | N=80 : études avec faux appels/40 |
| --- | ---: | ---: |
| Facteurs de taille et dispersions estimés, Wald normal | 3 (7,5 %) | 4 (10 %) |
| Facteurs de taille vrais, dispersions estimées, Wald normal | 3 (7,5 %) | 4 (10 %) |
| Facteurs de taille et dispersions vrais, Wald normal | 2 (5 %) | 2 (5 %) |
| Facteurs de taille et dispersions vrais, Wald avec référence t | 0 (0 %) | 1 (2,5 %) |

**Interprétation bornée** : l'effet observé de la normalisation paraît faible dans ces générateurs ; une partie de l'écart peut être liée à l'estimation de dispersion, à la référence du test de Wald ou à leur interaction. Il ne s'agit pas d'une décomposition causale indépendante puisque la modification d'un paramètre affecte la variance, les coefficients et les tests réajustés. Les campagnes #97 et #100 ont des graines différentes : ne pas confondre l'écart *entre campagnes* avec un contraste apparié. Le filtrage indépendant n'a pas changé les fréquences de découverte dans #100, mais cela ne valide pas son innocuité universelle.

**Précision Monte-Carlo insuffisante** : 40 réalisations par scénario donnent des IC binomiaux très larges ; **0/40** faux appel est compatible avec un taux allant jusqu'à **8,81 %** (borne exacte supérieure 95 %). Les critères « FDR contrôlé à 5 % », « Wald t supérieur » et « validité clinique » restent **non établis**. Aucun changement de référence statistique ne doit être appliqué sur la seule base de ces données.

Prochaine qualification à financer/calculer : générateurs hors famille NB, variations de dispersion et de corrélation, designs appariés réalistes, effets mixtes pour la puissance, cohortes biologiques indépendantes, puis réplication Monte-Carlo avec précision définie avant exécution. La reproductibilité informatique des fits, la calibration fréquentiste et la généralisation aux patients restent trois objectifs distincts.

## Matrice d'acceptation minimale

| ID | Cas | Vérification observable | Critère bloquant |
| --- | --- | --- | --- |
| R1 | RNA-seq indépendant, comptages entiers | Intention + contraste explicites, matrice/metadata concordants, covariables et lots identifiables | Pas d'inférence confirmatoire par repli sur navigateur |
| R2 | RNA-seq apparié | `subject_id` relie exactement deux conditions, `sample_id` identifie chaque prélèvement, modèle `~ subject_id + condition` ajustable | Une paire incomplète/doublon, modèle déficient ou contraste non défini bloque |
| R3 | Données transformées/non entières | Déclaration du type de données et branche de méthode adéquate | Interdire DESeq2 sur quantification non entière sans justification/méthode adaptée |
| R4 | Protéomique / métabolomique | QC manquants, blanks, dérive, filtres, imputations et limites explicités selon plate-forme | Ne pas inférer de validation biologique générale depuis un QC réussi |
| R5 | Mono-omique et trois omiques | Import réel → exécution moteur réel → rapport → ZIP → réimport | Aucun résultat synthétique maquillé en résultat expérimental |
| R6 | Sorties R | Tous les gènes/features ajustés, NA conservés, IDs et contraste inchangés ; `reference_r/` séparé des exports navigateur | Aperçu 100 lignes / sortie `too_large` ≠ export complet |
| R7 | Archivage | `methods.md`, paramètres, versions, checksums SHA-256, identifiants originaux, scripts de vérification | La corruption d'un octet et un fichier manquant doivent être détectés |
| R8 | Protection | Aucun fichier patient envoyé sans opt-in; ZIP privé explicite ; aucune fuite vers logs/télémétrie | Les données privées ne figurent jamais dans les traces publiques |
| R9 | Accessibilité débutant | Parcours FR/EN sur plan, appariement, groupe de référence, q-value et statut d'analyse | Éviter les contrôles avancés avant le parcours guidé |
| R10 | Reproductibilité externe | Jeu public figé, provenance SHA/versions, comparaison numérique ligne à ligne avec méthode indépendante ou appel de référence | Un benchmark démonstratif seul ne constitue pas validation indépendante |
| R11 | Déploiement | Test depuis l'URL Pages, base path, hydration, import CSV MS AUC, export ZIP et réimport | CI de branche verte mais application publique cassée = échec |

## Qualification statistique — règles non négociables

1. Séparer **erreur de câblage** (wrapper identique à l'appel R natif), **calibration statistique** (type I / FDR sous scénarios), **validité externe biologique** (jeu indépendant) et **validité prédictive** (split hors-échantillon / cohorte externe). Une preuve d'une catégorie ne vaut pas pour les autres.
2. Signaler pour chaque scénario N, taille du groupe, taux de filtrage, dispersion, batch/paired/correlation, pourcentage de NA, méthode et version des packages, seed, nombre de répétitions, faux positifs et incertitudes Monte Carlo.
3. Évaluer FDR réalisé et fréquence de découvertes séparément, avec IC binomiaux ou incertitude de simulation; éviter de conclure sur la moyenne de FDP seule dans les scénarios sans découvertes. Distinguer contrôle FDR et puissance.
4. Ne pas « corriger » une inflation FDR en changeant arbitrairement vers une loi t. #75 doit rester un diagnostic jusqu'à comparaison avec appels natifs, hypothèses, calibrage et erreurs de mise en œuvre.
5. Pour MOFA2/DIABLO et tout gain multi-omique : inclure chaque omique seule, fusion naïve, permutation des appariements, split au niveau sujet/cohorte, métriques de calibration, analyses par sous-groupe et contrôles d'informations fuyantes. Un résultat multi-omique inférieur doit rester visible.
6. Le jeu Airway (GSE52778) teste un design apparié connu. Confronter toutes les lignes, pas seulement les 100 premières, à un appel de référence; confirmer les attentes biologiques sans qualifier un effet isolé de sensibilité universelle.

## Artefacts requis pour chaque scénario

- `data_manifest.json` : provenance, téléchargements, empreintes, conditions de partage, structure des échantillons.
- `analysis_config.json` : plan, intention, contraste directionnel, filtres, méthodes, versions, seed.
- `qc.json` : tailles entrantes/sortantes, avertissements et blocages, décisions reproductibles.
- `reference_r/<method>.csv` : lignes exhaustives *si et seulement si* le backend les a effectivement produites.
- `results_browser/` : résultats d'exploration clairement séparés.
- `methods.md` : texte compréhensible, limites et détails suffisants pour Methods.
- `manifest.json` et vérification hors ligne de l'intégrité de l'archive.
- `validation_summary.json` : statut de chaque gate, version du code, commandement CI, preuves et échecs.

Ces noms constituent un **contrat cible** : ne pas prétendre qu'ils existent tous déjà dans la version publiée. La structure implémentée peut utiliser des chemins équivalents si le manifeste les décrit sans ambiguïté.

## Définition de « prêt pour les chercheurs »

- Un nouvel utilisateur peut réaliser le chemin complet avec une démo documentée puis ses propres fichiers, sans connaître `metadata` ni modèle linéaire à l'avance.
- Chaque analyse indique ce qui a été réellement exécuté (navigateur versus R de référence), le nombre d'échantillons indépendants, les exclusions et les limitations de validité.
- Les erreurs du design sont bloquées avant inférence, jamais converties silencieusement en une autre méthode.
- Le dossier d'export se vérifie, s'ouvre et conserve toutes les variables réellement analysées.
- Des jeux publics **externes** et des contrôles négatifs révèlent aussi les échecs, y compris les cas où intégrer plusieurs omiques n'améliore pas les résultats.
- Un test humain (débutant et expert) reste requis avant de revendiquer une accessibilité démontrée.

**Pas de fusion / déploiement « confirmatoire » avant passage documenté de cette matrice sur le SHA intégré.**
