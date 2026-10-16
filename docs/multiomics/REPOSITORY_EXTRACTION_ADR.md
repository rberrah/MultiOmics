# ADR — Extraire MultiOmics du dépôt PMx Explain avant un article

**Statut : proposition technique, migration non exécutée.**
**Périmètre : uniquement MultiOmics.** Aucun changement à l'outil PK/PD pédagogique ni aux autres projets.

## Pourquoi séparer les dépôts

Le dépôt `rberrah.github.io` regroupe actuellement le site pédagogique PMx Explain,
la route SvelteKit de MultiOmics, son moteur JavaScript, ses outils de
validation R/Bioconductor, ses jeux de références publics, ses workflows et
ses scripts d'exports. Cela allonge les CI, complique le versionnement
scientifique et rend difficile de citer uniquement MultiOmics.

**Avis : oui, créer un dépôt autonome MultiOmics avant la publication de
l'article logiciel ou méthodologique.** En revanche, ne pas déplacer les fichiers
sur `main` en une opération ni casser les URL existantes.

## Architecture cible

Un dépôt autonome (nom de travail `multiomics`, nom public non fixé) :

- `apps/web/` : interface SvelteKit mono-/multi-omique, démos,
  aide débutant et visualisations interactives ;
- `packages/analysis-js/` : statistiques exploratoires
  déterministes, visualisations et schémas de données ;
- `packages/research-export/` : contrat d'échange, SHA-256,
  ZIP, Rmd et sources de rejeu ;
- `r-engine/` : backend Plumber et méthodes de référence
  DESeq2, edgeR, limma-voom, fgsea, MOFA2, DIABLO ;
- `validation/` : jeux publics épinglés, cas négatifs,
  scripts R/JS et rapports FAIR des vérifications ;
- `docs/` : protocoles scientifiques, tutoriels, limites,
  sauvegardes de versions et futur lien omiques→PK/PD.

Le portail PMx Explain ne conserverait à terme qu'une entrée/menu
vers MultiOmics. Les fonctionnalités omiques sont développées et publiées
dans le nouveau dépôt, et non dupliquées.

## Contrats à préserver absolument

1. **URL existante** `/multiomics` : conserver une page d'accueil ou
   une redirection explicite vers le nouvel outil. Ne pas casser les
   anciens liens de cours et les démonstrations.
2. **Contrat des données** : manifeste d'échange versionné,
   identifiants sujets/prélèvements/assays, tissu/temps/units, type de
   mesure, NA/zéro, provenance et SHA.
3. **Reproductibilité** : ZIP contenant le code R réellement exécuté
   lorsqu'il existe, le notebook `ANALYSE_REPRODUCTIBLE.Rmd`,
   le lancement en une commande et les figures R comparables.
4. **Fiabilité** : ne jamais remplacer une méthode statistique native R
   par une approximation JS à l'insu du chercheur.
5. **Sécurité** : le site ne doit pas envoyer les données patients vers
   un service tiers par défaut ; backend R local/autorisé clairement défini.
6. **Versions** : chaque publication liée à un tag immuable du moteur R,
   du web, du schéma de données et des dépendances R/Bioconductor.
7. **Intégration future** : SWEET/Pilier A/Pilier B, ou autre
   consommateur PK/PD, utilise une sortie versionnée ; aucune dépendance
   directe au routage du portail.

## Inventaire de migration (à confirmer par graphe d'import automatisé)

| Origine actuelle | Destination provisoire |
| --- | --- |
| `src/routes/multiomics/**` | `apps/web/src/routes/**` |
| `src/lib/multiomics/**` | `packages/analysis-js/` et `packages/research-export/` |
| `src/lib/components/MultiomicsVisualizations.svelte` | `apps/web/src/lib/components/` |
| `src/lib/components/ExternalValidationPanel.svelte` | `apps/web/src/lib/components/` |
| `multiomics-engine/**` | `r-engine/` |
| `static/multiomics/**` | `apps/web/static/` |
| `scripts/*multiomics*` et tests associés | `validation/` |
| `.github/workflows/multiomics-*.yml` | `.github/workflows/` du nouveau dépôt |
| `docs/multiomics/**` | `docs/` |

Les imports vers `$lib/stores/language`, la gestion du `base` SvelteKit,
les styles globaux et les utilitaires partagés doivent être recensés
puis convertis en petits adaptateurs explicites. Les scripts ou
composants du portail non omiques ne doivent pas être copiés par défaut.

## Migration par étapes, sans interruption

1. **Geler une baseline validée** du site actuel : SHA `main`,
   golden datasets Airway, Nutrimouse, TCGA, jeux négatifs,
   parcours ZIP extrait et tests FR/EN.
2. **Créer le dépôt** et copier l'historique utile ou conserver une
   traçabilité Git vers le commit origine, avec licence, politique de
   contribution, fichier CITATION.cff une fois auteurs confirmés,
   documentation d'installation et version du schéma.
3. **Reconstruire les dépendances** : imports, CSS, chemins d'assets,
   variables R et JavaScript, CI, contenus et routes.
4. **Exécuter en parallèle** site actuel et nouveau dépôt sur les mêmes
   données publiques et comparer les résultats complets, les graphiques,
   la provenance du ZIP, les erreurs et le temps d'accès.
5. **Publier une première release indépendante** seulement si ces
   comparaisons passent et qu'un utilisateur peut lancer un notebook Rmd
   sans éditer le code.
6. **Modifier ensuite le portail** : carte et lien vers le site
   MultiOmics autonome ; maintenir une page de compatibilité
   `/multiomics` et consigner les anciennes URL.
7. **Retirer la copie historique** de l'ancien dépôt lors d'une PR distincte
   après validation du déploiement et des liens publics. Garder un
   pont de migration pendant une release.

## Critères bloquants pour la migration

- Pas de publication/renommage massif avant comparaison numérique
  des résultats R avec l'ancienne version.
- Pas de changement des URLs sans redirection testée desktop/mobile.
- Pas de suppression du backend R ou des test suites indépendantes.
- Pas de CI verte obtenue en désactivant/ignorant les tests lents ou
  en remplaçant des packages réels par des doubles synthétiques.
- Les limitations empiriques connues (calibration FDR, sensibilité
  normalisation de MOFA2, absence de supériorité DIABLO prouvée)
  restent publiques dans le nouveau dépôt.

**Décision proposée :** extraire après l'adoption du carnet Rmd et
une dernière validation du ZIP actuel. Rédiger l'article en parallèle,
mais citer la première release indépendante seulement après sa recette.

## Nouveau dépôt créé et déploiement (vérification du 16 octobre 2026)

- **Dépôt disponible :** [`rberrah/MultiOmics`](https://github.com/rberrah/MultiOmics), privé, accessible avec droits administrateur depuis le connecteur GitHub, actuellement vide.
- **Décision :** ne pas déplacer le code ni les données avant acceptation des tests du notebook `.Rmd` sur jeux biologiques publics et revue du périmètre scientifique. Le site `rberrah.github.io` continue à servir MultiOmics comme aujourd'hui.
- **Privé n'implique pas inaccessible aux visiteurs :** le navigateur reçoit les ressources statiques compilées publiées à partir du dépôt public du portail, jamais les sources privées directement.
- **Point bloquant à la future extraction :** le `GITHUB_TOKEN` d'une Action du dépôt public ne permet pas d'obtenir automatiquement le code de l'autre dépôt privé. Utiliser une GitHub App avec permissions minimales, ou construire et publier un artefact statique depuis le dépôt privé au moyen d'un workflow explicitement autorisé. Ne pas mettre de clé dans du JS client, ni exposer le code/données sensibles dans les artefacts ou logs publics.
- **Pages depuis le dépôt privé :** l'éligibilité dépend de l'offre GitHub. Pour éviter cette dépendance et conserver les URL, la voie recommandée est un build testable de MultiOmics suivi d'une publication statique dans `rberrah.github.io`, ou un service d'hébergement distinct avec URL/redirect contrôlés.
- **Validation de migration :** mêmes tables R natives, mêmes valeurs JS tracées, mêmes données NA/couleurs, ZIP/Rmd/HTML, tests navigateurs FR/EN et redirections; aucune suppression avant cette recette.

