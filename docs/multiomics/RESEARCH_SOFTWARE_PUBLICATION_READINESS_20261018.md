# MultiOmics — état de maturité et décision de publication scientifique

**Revue initiale du 18 octobre 2026 ; décision de gel mise à jour le 19 octobre 2026. Périmètre : dépôt autonome MultiOmics uniquement.**

## Trois décisions différentes

| Usage / décision | État vérifiable | Formulation scientifiquement autorisée |
| --- | --- | --- |
| Utilisation pour produire des résultats publiables dans un article de biologie | **Possible sous conditions** | Des jeux de données réels peuvent être traités ; archiver la version, le plan, les paramètres, les fichiers d'entrée, les sorties R réellement calculées, le rapport et la revue scientifique du cas |
| Première version citée et figée d'un logiciel de recherche | **Research v1 / 1.0.0 : gel logiciel accepté sous gate CI** | Le code autonome, le ZIP reproductible, les contrats natifs R, les qualifications publiques Airway/Nutrimouse/Sacurine/NKI/GSE191328 et la fermeture transitive R sont intégrés. L'identité définitive est le SHA `main` issu du merge de freeze ; un GitHub Release/DOI pérenne reste une étape de diffusion, pas une condition de validité du code gelé. |
| Publication méthodologique du logiciel dans une revue à évaluation par les pairs | **Dossier incomplet** | Il faut démontrer indépendance, intérêt scientifique comparatif, utilité externe et limites avec résultats quantitatifs et reproductibles ; ne pas confondre la réussite CI et une validité statistique universelle |

## Fonctions réelles disponibles

- Analyse transcriptomique, protéomique, métabolomique et multi-omique, interface novice FR/EN et sombre/clair, import de matrices/metadata et démos clairement synthétiques.
- Moteur R de méthodes statistiques de référence **local ou autorisé séparément**. La version statique du site GitHub Pages **n'exécute pas R**.
- Modes exploratoires et *confirmatory intent* explicitement séparés. Comparaisons de groupes, interaction groupe × temps et survie ne sont étiquetées natives que lorsque **toutes les couches** disposant de données ont leurs tableaux R complets vérifiés. Les autres outcomes feature-wise, la prédiction et l'intégration multibloc restent exploratoires pour leurs conclusions d'inférence.
- Export ZIP, sources, paramètres, provenance, méthodes, contrôles d'intégrité et rapport R Markdown quand le moteur compatible est disponible.
- Cohortes publiques externes et validations numériques réelles Airway, Nutrimouse, GEO et QC LC-MS, avec preuves de source et limites détaillées.

## Points non négociables avant une version « Research v1 »

**P0 — Fiabilité scientifique**

1. L'inflation **DESeq2 FDR** en nul global est significative dans trois générateurs indépendants (22/128, 23/128, 18/128 études faussement significatives) : [résultats persistants](RNASEQ_HIGH_REP_FDR_FINDINGS_20261018.md) et [issue #44](https://github.com/rberrah/MultiOmics/issues/44). Le moteur, l'UI et les rapports ne doivent pas promettre un FDR universel contrôlé à 5 %. Revue humaine et justificatif des méthodes restent obligatoires.
2. La qualification native **12 mzML réels Sacurine v3** est désormais exécutée : 12/12 sources MD5 vérifiées, 6 urines/3 pooled-QC/3 blancs, ordres source conservés, xcms natif, audit numérique indépendant et sensibilité préspécifiée aux paramètres (PR #48/#50). Les résultats défavorables sont conservés : 64,67 % des 2 850 features à RSD QC évaluable dépassent 30 %, 72,71 % des 535 ratios blanc/urine évaluables dépassent 20 %, et 2 600/3 135 ratios blanc restent non évaluables. **Le P0 résiduel** de l'[issue #47](https://github.com/rberrah/MultiOmics/issues/47) est une comparaison de prétraitement réellement indépendante/standards analytiques lorsque l'identification ou la précision métabolomique est revendiquée, ainsi qu'une politique de QC aval préspécifiée.
3. Séparer expressément : validité d'exécution R, calibration fréquentiste, plausibilité biologique, validation externe et utilisation clinique. Pas de revendication de validation universelle ni médicale.
   Les voies longitudinales et survie projetées en R natif conservent leurs workflows numériques indépendants et disposent désormais d'une qualification externe réelle sur NKI/DMFS et GSE191328. Cette qualification établit la parité d'exécution pour ces contrats précis, pas la calibration fréquentiste ni la validité biologique universelle.
4. Faire relire au moins quelques analyses réelles par un statisticien/biologiste indépendant, y compris échecs et données manquantes, et effectuer une validation d'utilisabilité par des novices. L'utilisateur prend en charge les tests humains.

**P1 — Reproductibilité d'une version citée**

5. **Réalisé pour le code Research v1** : version 1.0.0, manifeste, citation, versions npm/R/Bioconductor, conteneur immuable, inventaire complet de 171 packages R, ZIP autonome et gate de réexécution. Le SHA exact du merge final devient l'identité du gel. La création d'un GitHub Release/tag/DOI est une étape de diffusion ultérieure ; ne jamais le nommer « scientifiquement validé ».
6. Séparer la version public navigateur, les options de déploiement du moteur R et le niveau de sécurité/infrastructure ; fournir une installation local R sans données patients envoyées à un serveur non autorisé.
7. Fournir un exemple biologique complet documenté montrant le trajet **input original → QC → design → vraies sorties R → figures → ZIP/Rmd → reproduction indépendante**, avec toutes les lignes de résultats et comparaisons à l'appel package natif.

**P2 — Publication d'un article sur MultiOmics**

8. Positionner l'apport scientifique **usability + composabilité multi-omique + transparence d'inférence + reproductibilité** et non un moteur statistique nouveau ou une calibration 5 % non démontrée ; confronter à un pipeline R/Bioconductor adapté dans plusieurs cohortes publiques indépendantes.
9. Démontrer un cas d'usage externe véritable, les limites des données manquantes et le taux de concordance numérique, puis documenter l'usage par des chercheurs autres que le développeur.
10. **Choisir la revue selon la licence.** Le dépôt est actuellement sous licence **PolyForm Noncommercial**, publique/source-available mais **pas OSI open-source**. JOSS exige une licence OSI et une histoire publique de développement d'au moins six mois, en plus d'un logiciel complet et d'une preuve d'adoption/recherche. **Ne pas remplacer la licence sans décision explicite du détenteur des droits.** Une publication scientifique *utilisant* MultiOmics n'est pas la même soumission qu'un article *sur* le logiciel.

## Critère de rédaction d'un article issu d'une analyse MultiOmics

La section Methods doit expliciter versions/SHA, origine et taille du jeu, contrôles indépendants des matrice et groupes, modèle statistique exactement exécuté, contraste, méthode de correction multiple, résultats QC, détails des exclusions/NA, métriques retenues et revue de sensibilité. Les figures/p-values « exploratoires browser-only » ne doivent pas être décrites comme une référence R ou comme certifiées.

**Décision au 19 octobre 2026 :** le **gel Research v1 / 1.0.0 du logiciel et de ses contrats de preuve est acceptable** si le SHA final passe toutes les CI de release. Cela ne rend pas MultiOmics universellement confirmatoire. Les issues #44 (calibration DESeq2), #47 (validation analytique LC-MS au-delà de l'exécution xcms réelle) et #20 (transférabilité biologique des scores vers Omics2PK) restent des limites scientifiques explicites et non des fonctionnalités cachées comme validées. Les revues humaines et tests utilisateurs restent externes au gel logiciel.
