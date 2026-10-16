# MultiOmics — qualification d'usage scientifique et de résultats exportables

**13 octobre 2026 — protocole déclaré avant résultats des nouvelles campagnes.**
Ces contrôles sont indépendants du manuscrit préliminaire (PR #76) et des
branches de Codex #52/#54. Les calculs restent déterministes ou utilisent des
bibliothèques R de référence ; aucune IA générative n'interprète les données.

## Objectif

Une démonstration qui produit une figure n'est pas une preuve suffisante pour
une analyse destinée à un article. Un utilisateur doit pouvoir **analyser
des fichiers, obtenir des tableaux structurés, documenter les méthodes,
archiver les entrées exactes, rouvrir le projet et recalculer les sorties**
sans transformation silencieuse.

## P0 — recette navigateur complète

`tests/e2e/multiomics-research-output-acceptance.spec.js` exécute une vraie
session SvelteKit dans Chromium, avec les données de démonstration
synthétiques et sans moteur R simulé :
- RNA seul, protéines seules, métabolites seuls : chacune fournit des lignes
  et des p/q numériques finies dans [0,1] quand applicables, et une
  déclaration honnête de l'absence de moteur R.
- Trois omiques : tables séparées pour chaque omique, rapport de méthodes,
  aucune certification de validité biologique/clinique.
- Un export complet de résultats, un export de méthodes sans matrices brutes
  et un projet portable contenant **exactement les entrées et leurs SHA-256**.
- Import du projet dans une session vierge, aucune injection de résultats
  archivés en lieu et place d'une nouvelle analyse, nouvel ajustement
  déterministe et comparaison stricte feature/effet/p/q.
- Blocage d'appels inattendus à Reactome, ChEBI et à un R local dans le
  scénario autonome (une consultation externe doit être volontaire).
- L'analyse sait produire les données d'un article **exploratoire** ;
  cette recette n'établit pas un test statistique confirmatoire validé.

## P0 — réelle donnée biologique avec comparateurs honnêtes

Source officielle `mixOmics-org/mixOmics` au commit
`ef3e760526623d9e91e945e4b50be48d764efc40`, jeu
`breast.TCGA$data.train`, échantillons Her2 et LumA avec transcriptome
et protéome. **Une étude biologique publique réelle différente de Nutrimouse,
mais pas une nouvelle cohorte externe indépendante** : l'échantillon est
déjà traité par le package source, et cinq sous-échantillonnages répétés
réutilisent le même ensemble initial.

L'outil de référence `mixOmics` calcule quatre bras sur les **mêmes
sujets tenus à l'écart** :
1. DIABLO ARN + protéines ;
2. sPLS-DA ARN seul ;
3. sPLS-DA protéines seules ;
4. DIABLO ARN + protéines **mal appariées volontairement** (contrôle négatif).

5 partitions stratifiées préspécifiées, un composant et 10 features
par omique, normalisation, imputation et filtrage ajustés sur les
sujets d'entraînement seulement. L'hyperparamétrage fixé n'est
pas présenté comme optimal. Un bras à deux vues sélectionne jusqu'à
20 features contre 10 pour chacun des bras mono-omiques :
cela crée un avantage de capacité potentiel, explicitement signalé.
**Les performances sont descriptives et aucun p-value de supériorité
n'est fabriqué en traitant les observations de folds répétés comme
indépendantes.** Les résultats négatifs restent archivés.

Les données réellement exécutées mesurent la valeur de l'intégration
dans cette étude ; elles ne valident pas un nouveau biomarqueur, ne
permettent pas de généraliser entre hôpitaux ni de conclure à une
supériorité universelle des méthodes multi-omiques.

**Premier résultat de l'exécution réelle (13 octobre 2026)** —
sur les cinq partitions TCGA, DIABLO a une BER descriptive de 17,5 %,
ARN seul de 10,9 %, protéines seules de 17,5 %, et le contrôle
désapparié de 10,3 %. **Aucun bénéfice multi-omique n'est démontré**.
La [question de qualification biologique #78](https://github.com/rberrah/rberrah.github.io/issues/78)
prescrit l'investigation de l'appariement, une meilleure harmonisation
de la complexité des bras et une réplication externe. Une CI verte ne
supprime pas cet échec de preuve de valeur ajoutée.

## P0 — validité statistique : puissance ET faux positifs

Second diagnostic indépendant des 400 études / 1200 fits (#52) et
des 300 nulles / 900 fits (#75). Pour 60 études par cas :
- deux groupes indépendants 20/20 ou 12/28 ;
- nul global **ou** 40 vrais signaux préspécifiés parmi 400 gènes ;
- matrice NB avec dispersions hétérogènes, taille de bibliothèque
  variable et lot mesuré non confondu complètement ;
- test DESeq2 de Wald normal, puis test avec référence *t* sur la même
  matrice ajustée ; identiques filtres et contrastes ;
- FDR empirique, puissance sur effets plantés, MAE des effets retrouvés,
  proportions de rejet et intervalles (Clopper–Pearson nul complet ;
  bootstrap par étude pour vérités mixtes).
- Erreurs de convergence écrites explicitement, pas remplacées par zéro
  ni exclues silencieusement du dénominateur.

**Règle décisionnelle préalable :** AUCUNE substitution automatique
`useT=TRUE` n'est autorisée sur ce seul diagnostic. Un meilleur
contrôle des faux positifs peut s'accompagner de perte de puissance.
Un autre générateur avec mécanisme différent et une cohorte biologique
indépendante restera nécessaire.

## Évaluation publique réelle et déploiement

Le site public GitHub Pages n'héberge pas à lui seul un moteur R.
L'exécution R réelle dépend du backend disponible et du mode choisi.
En complément de la validation de la branche, le workflow existant
`multiomics-public-route-smoke.yml` vérifie déjà hydratation
et démo publique desktop/mobile **après** le déploiement Pages :
il ne faut pas le confondre avec la validation d'un code encore
présent uniquement dans une PR.

Pour que le logiciel produise des résultats utilisables dans un article
tiers, chaque auteur devra conserver : vrais fichiers d'entrée,
métadonnées, traitements/exclusions, modèle, contraste, famille de
correction multiple, versions exactes, résultats complets, contrôles
négatifs, figures, et provenance des API externes.

Les exports JSON permettent d'archiver les entrées et méthodes,
mais les flux R, les données externes et les formats bruts MS binaires
demandent toujours une archive additionnelle.

## Conditions avant de revendiquer « professionnel et publiable »

- Démontrer qu'une **nouvelle équipe** peut utiliser le site pour
  une étude réelle, avec suivi des erreurs et vérification du
  comportement de la plateforme.
- Publier les résultats défavorables et les métriques, sans
  modifier les seuils selon leur apparence.
- Fournir un conteneur R complet pour les méthodes supportées ;
  les tests environnementaux de #73 portent pour l'instant sur
  trois méthodes RNA-seq et une simulation synthétique.
- Obtenir l'amélioration multi-omique, si elle existe, sur un
  jeu de données externe avec un *pipeline gelé* et **des patients
  indépendants de ceux utilisés pour tout réglage**.
- Réaliser l'essai de terrain avec des débutants, tests clavier,
  lecteurs d'écran, export/réimport, et données réelles de formats
  hétérogènes. Un succès Playwright ne remplace pas ces observations.

## Revue scientifique et prépublication

La présentation sous forme de logiciel de recherche peut viser
une *Application Note* à Bioinformatics Advances **si** l'on documente
une utilité pratique démontrée, le comparatif aux autres outils et
un résultat biologique reproductible, ou JOSS lorsque l'installation,
les tests, l'historique public et l'impact logiciel sont convaincants.
Un préprint bioRxiv catégorie Bioinformatics viendra **après**
l'analyse des preuves et l'accord des auteurs, pas avant.

La revue JOSS vérifie aussi un historique public de développement
**d'au moins six mois** avec de vrais signes de développement continu.
L'ancienneté du dépôt d'hébergement de PMx Explain ne prouve
pas celle de l'outil MultiOmics lui-même.

Aucun algorithme n'a été ajouté au moteur pour fabriquer une
performance et aucune modification des travaux en cours de Codex.
