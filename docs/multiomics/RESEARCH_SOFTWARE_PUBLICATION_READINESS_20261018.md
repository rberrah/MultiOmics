# MultiOmics — état de maturité et décision de publication scientifique

**Revue du 18 octobre 2026. Périmètre : dépôt autonome MultiOmics uniquement.**

## Trois décisions différentes

| Usage / décision | État vérifiable | Formulation scientifiquement autorisée |
| --- | --- | --- |
| Utilisation pour produire des résultats publiables dans un article de biologie | **Possible sous conditions** | Des jeux de données réels peuvent être traités ; archiver la version, le plan, les paramètres, les fichiers d'entrée, les sorties R réellement calculées, le rapport et la revue scientifique du cas |
| Première version citée et figée d'un logiciel de recherche | **Préparation avancée, non encore figée** | Le code autonome, les jeux de référence, le ZIP reproductible et les vérifications numériques existent ; un tag immuable, un dépôt pérenne avec DOI et des scénarios utilisateurs indépendants restent à finaliser |
| Publication méthodologique du logiciel dans une revue à évaluation par les pairs | **Dossier incomplet** | Il faut démontrer indépendance, intérêt scientifique comparatif, utilité externe et limites avec résultats quantitatifs et reproductibles ; ne pas confondre la réussite CI et une validité statistique universelle |

## Fonctions réelles disponibles

- Analyse transcriptomique, protéomique, métabolomique et multi-omique, interface novice FR/EN et sombre/clair, import de matrices/metadata et démos clairement synthétiques.
- Moteur R de méthodes statistiques de référence **local ou autorisé séparément**. La version statique du site GitHub Pages **n'exécute pas R**.
- Modes exploratoires et *confirmatory intent* explicitement séparés. Les comparaisons de groupes ne sont étiquetées natives que lorsque **toutes les couches** disposant de données ont des tableaux R complets vérifiés. Longitudinal/survie restent exploratoires dans leur projection UI.
- Export ZIP, sources, paramètres, provenance, méthodes, contrôles d'intégrité et rapport R Markdown quand le moteur compatible est disponible.
- Cohortes publiques externes et validations numériques réelles Airway, Nutrimouse, GEO et QC LC-MS, avec preuves de source et limites détaillées.

## Points non négociables avant une version « Research v1 »

**P0 — Fiabilité scientifique**

1. L'inflation **DESeq2 FDR** en nul global est significative dans trois générateurs indépendants (22/128, 23/128, 18/128 études faussement significatives) : [résultats persistants](RNASEQ_HIGH_REP_FDR_FINDINGS_20261018.md) et [issue #44](https://github.com/rberrah/MultiOmics/issues/44). Le moteur, l'UI et les rapports ne doivent pas promettre un FDR universel contrôlé à 5 %. Revue humaine et justificatif des méthodes restent obligatoires.
2. Terminer et documenter la qualification **12 mzML réels Sacurine v3** (6 biologiques, 3 QC, 3 blancs, vrais ordres et empreintes source) : [issue #47](https://github.com/rberrah/MultiOmics/issues/47) et [protocole](SACURINE_ZENODO_V3_12RAW_XCMS_QA_20261018.md). La source seule ne prouve ni extraction correcte ni validité biologique.
3. Séparer expressément : validité d'exécution R, calibration fréquentiste, plausibilité biologique, validation externe et utilisation clinique. Pas de revendication de validation universelle ni médicale.
4. Faire relire au moins quelques analyses réelles par un statisticien/biologiste indépendant, y compris échecs et données manquantes, et effectuer une validation d'utilisabilité par des novices. L'utilisateur prend en charge les tests humains.

**P1 — Reproductibilité d'une version citée**

5. Figer une version avec tag GitHub Release, SHA exact et versions npm/R/Bioconductor, plus ZIP autonome, test de réexécution de bout en bout et empreintes des jeux publics. Ne pas publier un tag « scientifiquement validé » si le dossier de qualification du SHA n'est pas vert.
6. Séparer la version public navigateur, les options de déploiement du moteur R et le niveau de sécurité/infrastructure ; fournir une installation local R sans données patients envoyées à un serveur non autorisé.
7. Fournir un exemple biologique complet documenté montrant le trajet **input original → QC → design → vraies sorties R → figures → ZIP/Rmd → reproduction indépendante**, avec toutes les lignes de résultats et comparaisons à l'appel package natif.

**P2 — Publication d'un article sur MultiOmics**

8. Positionner l'apport scientifique **usability + composabilité multi-omique + transparence d'inférence + reproductibilité** et non un moteur statistique nouveau ou une calibration 5 % non démontrée ; confronter à un pipeline R/Bioconductor adapté dans plusieurs cohortes publiques indépendantes.
9. Démontrer un cas d'usage externe véritable, les limites des données manquantes et le taux de concordance numérique, puis documenter l'usage par des chercheurs autres que le développeur.
10. **Choisir la revue selon la licence.** Le dépôt est actuellement sous licence **PolyForm Noncommercial**, publique/source-available mais **pas OSI open-source**. JOSS exige une licence OSI et une histoire publique de développement d'au moins six mois, en plus d'un logiciel complet et d'une preuve d'adoption/recherche. **Ne pas remplacer la licence sans décision explicite du détenteur des droits.** Une publication scientifique *utilisant* MultiOmics n'est pas la même soumission qu'un article *sur* le logiciel.

## Critère de rédaction d'un article issu d'une analyse MultiOmics

La section Methods doit expliciter versions/SHA, origine et taille du jeu, contrôles indépendants des matrice et groupes, modèle statistique exactement exécuté, contraste, méthode de correction multiple, résultats QC, détails des exclusions/NA, métriques retenues et revue de sensibilité. Les figures/p-values « exploratoires browser-only » ne doivent pas être décrites comme une référence R ou comme certifiées.

**Décision actuelle :** utilisable pour la recherche avec qualification au cas par cas. Pas encore légitime de vendre le site, son backend et tous ses scénarios comme un logiciel *universellement confirmatoire*. Le parcours vers Research v1 est limité principalement par la calibration inférentielle, la validation complète LC-MS, l'usage indépendant réel et la publication d'une version figée reproductible.
