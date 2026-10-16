# MultiOmics — audit de migration autonome (16 octobre 2026)

## Périmètre exact et intégrité de copie

- Dépôt source : `rberrah/rberrah.github.io`.
- **Référence figée** : `dbe503f98cfde460f95a60617f3d743adb6b915b` (PR #112 déjà fusionnée).
- Dépôt cible : `rberrah/MultiOmics` ; branche `migration/import-validated-multiomics-20261016`.
- **244 fichiers liés aux omiques transférés**, comprenant 206 fichiers actifs (moteur JS, pages, visualisations, R, données publiques, documentation, tests) et 38 workflows historiques.
- Vérification Git blob : **242/244** identiques à la source avant ajustements autonomes ; les deux écarts initialement intentionnels portaient sur les commandes de clonage R dans la page d'analyse et le test supplémentaire de disposition UI. Le lien R de cette page a ensuite été réorienté vers le nouveau dépôt.
- Les anciens workflows sont conservés **à titre documentaire** dans `migration/legacy-workflows/`. Aucun de ces 38 workflows n'est activé directement : plusieurs reposaient sur les pages et outils PMx Explain absents ici.

## Travail d'autonomisation réalisé

- `src/routes/+layout.svelte`, `src/routes/+page.svelte` et `src/app.html` : shell indépendant, page d'accueil, navigation et structure HTML.
- `src/lib/stores/language.js` : état FR/EN indépendant des traductions PMx Explain.
- `src/lib/styles/theme.css`, ressources statiques et fontes locales : styles et scripts de guidage conservés.
- `package.json`, `package-lock.json`, `svelte.config.js`, `vite.config.js` et configuration Playwright : point d'entrée et build indépendant des générateurs TDM/chapitres PMx.
- `svelte.config.js` utilise des URL à dossier (`trailingSlash: 'always'`) et `BASE_PATH` paramétrable, nécessaire à l'hébergement projet GitHub Pages `/MultiOmics/`.
- `scripts/test_standalone_migration.mjs` : audit statique local des routes, scripts de guidage, boutons de démonstration, Rmd et sources R. Exécuter `npm run test:migration`.

## Commandes d'acceptation à exécuter avant fusion

```sh
npm ci
npm run test:migration
npm run check
npm run test:multiomics
npm run build
BASE_PATH=/MultiOmics npm run build
npx playwright install chromium
npm run test:e2e
Rscript multiomics-engine/tests/backend_smoke.R
Rscript multiomics-engine/tests/browser_reference_agreement.R
```

**Statut : les commandes ne sont pas encore certifiées sur le nouveau dépôt.**
L'ancien `package-lock.json` conserve des métadonnées de nom/version `pmx-explain` : synchroniser avec `multiomics` lors de la recette npm, sans changer les versions scientifiques silencieusement.

Un workflow GitHub Actions natif à ce nouveau dépôt reste à activer après vérification des droits d'écriture GitHub Actions. Les tests historiques du dépôt source ne se substituent pas aux tests du nouveau site.

## Critères de validation scientifique et de sécurité

1. **RNA-seq Airway/GSE52778** : rejouer le cas public apparié et comparer les 17 199 gènes au DESeq2 natif, y compris NA, contrastes, ordre, filtrage et versions. Le résultat initial était vérifié **dans le dépôt source** ; **pas encore rejoué ici**.
2. **Figures R** : comparer les valeurs, échelles z-score, NA, légendes et ordres d'échantillons entre le nouveau site et R ; figures comparables scientifiquement, non identiques pixel par pixel.
3. **R Markdown/ZIP** : construire un ZIP depuis ce dépôt, le décompresser sans l'éditer, exécuter `Rscript run_all.R` et `Rscript render_notebook.R`, vérifier HTML, CSV exhaustifs, SHA et rapports de comparaison.
4. **Modes analytiques** : JS = exploration, parfois simplifiée/approximative ; moteur R local = méthodes de référence sous hypothèses vérifiées. Aucun basculement automatique vers JS en mode confirmatoire.
5. **Cas négatifs** : contrôler indépendance des sujets, batch, appariement, missingness et échelle des données. Respecter les limitations publiées de calibration FDR, de MOFA2 et DIABLO ; aucune certification universelle.
6. **Déploiement public** : tester réellement la page GitHub Pages autonome (FR/EN, ordinateur/mobile, hydratation, 4 boutons, démo et CSV MS-AUC). Les sites statiques ne peuvent pas exécuter R côté serveur : backend local/autorisé explicitement.
7. **Confidentialité** : vérifier l'absence d'envoi de matrices patient à des tiers sans consentement, de journaux sensibles et de données cliniques dans les artefacts CI.

## Continuité du site

Le site déployé `https://rberrah.github.io/pharmacometrie/multiomics/tool/` **reste propriétaire de la route publique existante** pendant le développement ; aucune suppression ni redirection du portail vers le nouveau dépôt tant que la qualification autonome et le déploiement n'ont pas réussi.

## Résultat actuel

**Copie du périmètre MultiOmics terminée, première autonomie du projet réalisée ; validation de build/CI et publication autonome restent ouvertes.** L'article logiciel et l'extension omiques → PK/PD attendent une release figée et validée, pas seulement une copie des sources.
