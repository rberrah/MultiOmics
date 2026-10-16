# Audit technique et scientifique — MultiOmics autonome
**Date : 16 octobre 2026 · Périmètre : rberrah/MultiOmics**

Statut : audit de code, des contrats de CI, des artefacts et des preuves expérimentales accessibles. Les qualifications déclarées dans le tableau ci-dessous concernent les designs et commits explicitement testés, **pas une certification statistique universelle**. Le rendu public définitif dépend aussi du déploiement PMx, hors de ce dépôt.

## Synthèse

MultiOmics est une application autonome SvelteKit avec moteur d'exploration déterministe en navigateur et moteur de référence R/Bioconductor **optionnel et distinct**. Elle accepte 1 à 3 couches transcriptomiques / protéomiques / métabolomiques et quatre démonstrations synthétiques, avec garde-fous de design, export d'un projet scientifique, dossiers de provenance et contrôle de reproductibilité.

Le découplage du repo principal est techniquement atteint. Le navigateur **n'exécute pas DESeq2/xcms/MOFA2/mixOmics** à distance et ne doit pas présenter ses calculs JS comme leurs équivalents scientifiques.

## Inventaire observé

| Domaine | Évidence | État et limites |
|---|---|---|
| Source autonome | `package.json`, build SvelteKit 5, composants/modules propres | Fonctionnel, pas d'import obligatoire de PMx |
| Applications dans le même dépôt | `/multiomics/`, `/multiomics/tool/`, `/multiomics/learn/`, `/multiomics/disciplines/`, `/multiomics/resources/` | Deux overlay GitHub Pages, mais **déploiement final PMx à vérifier séparément** |
| Identité visuelle | Exactement même SHA blob `src/lib/styles/theme.css` que PMx (`6a8450d...`) | Les tokens sont partagés ; l'ancien header n'avait ni monogramme ni nav mobile PMx |
| Header PR #19 | Monogramme `Om` de 34px, terracotta, Source Serif, 4 liens internes, retour PMx/site principal | Modifié en branche, Playwright + double URL à qualifier avant fusion |
| FR/EN, clair/sombre | Stores partagés `pharmacometrie-language` et `pk-theme`, anciens préférences rétrocompatibles | Vérifier après assemblage final, à 320/390 px et clavier |
| Tests navigateur | ~30 spécifications e2e dans `tests/e2e` | **Précédemment seuls deux tests e2e** dans le workflow principal ; PR #19 élargit la couverture à l'interface, un seul omic, la provenance et la prudence scientifique |
| JavaScript scientifique | Tests déterministes, design checks, null models, nested CV, FDR, identifiants | Exploratoire uniquement ; la sélection de méthodes et les hypothèses dépendent du design |
| R natif | DESeq2, edgeR QL, limma-voom, limma, modèles mixtes, fgsea, MOFA2, DIABLO, xcms | Implémentations réelles et CI ; **ne pas inférer une calibration générale** des p/q-values |
| Qualification transcriptomique | Airway/GSE52778, vraie exécution R/Rmd, contrôle numérique | Réussite pour le design/source/test fixés, pas toutes les expériences |
| Qualification intégration | Nutrimouse réel, MOFA2/DIABLO, permutations négatives, ablations/single-view | Pas de gain multi-omique systématique démontré ; annotations/facteurs potentiellement techniques |
| Stabilité d'environnement | R 4.5.2 + Bioc 3.22, 171/206 transitive package closure, Python MOFA2 3.12.12/18 distributions | Versions et empreintes observées ; **pas d'installateur R offline à bibliothèques figées** |
| MS brute | faahKO réel (4 fichiers), source publiée Sacurine/Zenodo 1346742, 2 urines + 2 blancs, vrai xcms | 4,799 groupes dont 310 rapports blanc/urine calculables ; 158 ≥ 20 % (screen, non FDR) |
| MS QC externe | Deux expériences drosophiles, 8 panneaux, 429 features mesurées | QC descriptifs externes, pas le dépistage des contaminants sur des fichiers bruts |
| Étude Sacurine QC | PR #18, 9 mzML originaux Zenodo 3244991, 965.5 Mo, MD5 tous authentifiés | L'exécution R a détecté 101,526 pics / 3,347 features ; **résultat final QC en correction sur l'appariement assay_id**. Ne pas déclarer sa validation terminée |
| MetaboLights MTBLS404 complet | Source W4M expirations connexion, EBI 404 sur chemins testés | Non acquis ; tests d'archive offline verts mais ne valent pas un test brut |
| Reproductibilité de partage | ZIP portable, empreintes SHA-256, licences, source de citation, R Markdown scientifique | Distinguer ZIP autonome browser-only et vrai moteur R ; bien publier un artefact de release pérenne |
| Licence | PolyForm Noncommercial 1.0.0 + contenus CC BY-NC-SA ; anciennes licences préservées | Source disponible ; **ne pas écrire « open source OSI »** |
| Validation humaine | Parcours pédagogiques en Playwright mais absence de protocole novices/experts indépendants | Non démontrée |

## Actions ordonnées

### P0 — release et intégrité
1. Finir PR #19 : header PMx cohérent, logo omics reconnaissable, navigation mobile accessible, liens PMx/site principal, Responsive 320/390/1280 px, langue et thème, y compris ZIP portable et deux overlays URL. **Ne pas modifier le repo PMx depuis ici.**
2. Tester explicitement le shell dans `.github/workflows/standalone-ci.yml`; ce n'était pas inclus auparavant. Ajouter un smoke de single-omic, provenance, avertissements de scientific summary pour réduire le risque d'exposer du contenu incorrect.
3. Fusionner uniquement après dernier SHA vert ; **puis** communiquer le SHA immuable au PR d'intégration PMx. PMx doit vérifier le site assemblé (sitemap, 404, CSS, hydratation) et conserver les URLs héritées.
4. Lever le blocage scientifique PR #18 : vérifier l'ordre source de 9 chromatogrammes, aligner `assay_id` et le tableau `xcms::featureValues`, puis confirmer les QC RSD réellement estimables ; conserver les manquants sans imputation ni p-value inventée.

### P1 — méthode et valorisation
5. Réaliser une comparaison **à design égal** des résultats de référence contre Workflow4Metabolomics sur une cohorte LC-MS externe avec QC, blancs, acquisition et métadonnées vérifiées. Définir critères de concordance avant l'analyse.
6. Renforcer la matrice de qualifications statistiques (simulations nulles, faible effectif, batch confondu, profils incomplets, splits par sujet ; IC Monte Carlo, power et FDR par méthode). Conserver les échecs.
7. Ajouter au protocole de validation une étude multi-laboratoire indépendante. Une seule source publique/famille Sacurine ne démontre pas la transférabilité universelle.
8. Tests utilisateurs externes : minimum tâches mono-/multi-omics, import CSV, compréhension des limites, export/rejeu R, cas d'échec ; temps, taux de succès, erreurs de manipulation, version FR/EN.

### P2 — durcissement
9. Travailler vers un `renv.lock` ou un miroir binaire réellement reconstructible en plus de la vérification SHA de l'environnement observé ; cela n'est pas déjà acquis.
10. Limiter l'endettement lié aux warnings Svelte (sélecteurs inutilisés, attributs non accessibles, tooltips non interactifs), préserver le contraste WCAG clair/sombre et accélérer la volumétrie de graphiques.
11. Maintenir un manifeste de diffusion indiquant SHA du shell, SHA du moteur JS/R, artefacts référentiels, limitations et versions du navigateur.

## Dossier de publication — état

**Possible à préparer :** article *software/methods* présentant architecture, expérience novice, politiques de garde-fous, distinctions browser/reference R, tests reproductibles, comparateurs numériques et études de cas publiques.

**Prématuré à affirmer :** « pipeline confirmatoire universel », supériorité multi-omique, contrôle FDR global en toute situation, application clinique validée, robustesse interlaboratoire démontrée. Pour un article de **validation** plus fort : tests utilisateurs indépendants, comparateurs de référence, vrai jeu externe avec QC + blancs + designs complets, validation de transport et rapports de toutes les analyses négatives.

## Coordination de dépôts

- Source unique du code MultiOmics : `rberrah/MultiOmics`.
- Code PMx Explain et déploiement Pages : agent dédié, **aucune mutation depuis cette branche**.
- Séparer PR science xcms (#18) et PR UI (#19) : pas de mélange de données ni de commits non qualifiés dans le pin portail.
- Preuve du site réellement servi non vérifiable par une CI locale de source seule. Ne pas confondre `success` build et production validée.
