# MultiOmics RNA-seq — vérification des chemins de normalisation (15 octobre 2026)

## Résultat effectivement obtenu, avec références reproductibles

Source : [exécution GitHub Actions du 15 octobre 2026](https://github.com/rberrah/rberrah.github.io/actions/runs/38704282110) (PASS, 1 tentative), [artefact brut avec 5 fichiers CSV/texte](https://github.com/rberrah/rberrah.github.io/actions/runs/38704282110/artifacts/11958683082).

**Plan déjà exécuté** : trois scénarios (6 contre 6, 12 contre 28 avec batch inégal, 20 contre 20 avec signaux rares), 80 **jeux de données indépendants par scénario**, trois moteurs R réellement installés (DESeq2, edgeR QL, limma-voom), trois états sur chaque jeu : référence/native, filtrage strict avec facteurs d'échelle figés sur la référence, puis filtrage strict avec renormalisation native. Contraste \`treated - control\`, mêmes matrices, sujets, vérité artificiellement injectée et batch dans le triplet. **Les trois états sont des refits réels ; ils n'isolent pas un effet causal pur du filtrage.**

Les chiffres ci-dessous sont les **différences moyennes appariées** « strict + renormalisation » moins « référence ». L'augmentation du FDR empirique sous vérité mixte est défavorable ; une réduction de FDR seule n'est **pas** une réussite si la puissance recule. Pourcentage en **points**, et non en pourcentage relatif.

| Scénario | Méthode | Δ puissance (points) | Δ FDP moyen / FDR empirique (points) |
| --- | --- | ---: | ---: |
| 6 vs 6 | DESeq2 | −0,125 | −3,438 |
| 6 vs 6 | edgeR QL | −0,250 | −2,500 |
| 6 vs 6 | limma-voom | −0,125 | −0,625 |
| 20 vs 20, signaux rares | DESeq2 | −4,063 | −2,708 |
| 20 vs 20, signaux rares | edgeR QL | −2,813 | +0,729 |
| 20 vs 20, signaux rares | limma-voom | −1,250 | +2,917 |
| 12 vs 28, batch inégal | DESeq2 | −2,188 | +0,188 |
| 12 vs 28, batch inégal | edgeR QL | −1,563 | −0,028 |
| 12 vs 28, batch inégal | limma-voom | −0,375 | −0,729 |

Le résultat **n'autorise pas** un pré-filtrage strict automatique, ni une recommandation globale fondée sur les p-values produites : l'effet dépend du scénario et du moteur. **Important :** ce FDR est la moyenne des proportions de fausses découvertes sur simulations *avec signaux réels et nuls* ; ce n'est ni le FWER sous nul global du protocole #100, ni une estimation de FDR en patients.

## Nouvel instrument de qualification : intervalles appariés, sans surinterprétation

\`rnaseq_scale_path_paired_uncertainty()\` calcule maintenant, **avant toute décision de production**, les intervalles percentile bootstrap à 95 % des différences moyennes :

- Unité d'échantillonnage : **jeu de données indépendant**, jamais les gènes ; au sein d'un couple scénario–méthode, le même tirage d'indices est partagé par les trois bras ; les rééchantillonnages entre méthodes sont indépendants et ne constituent pas un contraste interméthode.
- Les trois chemins sont comparés : total (strict-native − référence), filtrage/refit à échelles figées (strict-frozen − référence), renormalisation additionnelle (strict-native − strict-frozen).
- Deux critères distincts : puissance et FDP moyen sur vérité mixte. Tous les contrastes vérifient l'identité arithmétique exacte et la concordance des graines/IDs de simulation.
- 2 000 tirages de bootstrap déterministes, seed 20261015, restauration de l'état aléatoire appelant, sortie \`rnaseq_norm_lock_paired_uncertainty.csv\`, résultats détaillés par jeu conservés.
- Intervalles descriptifs **conditionnels au générateur simulé**, non des bornes de conformité du FDR à 5 %, et non une validation externe biologique ou clinique. La CI réussie implique uniquement cohérence d'exécution.

Le statut de certification reste \`FALSE\` ; aucun changement du moteur de référence n'est autorisé par ce diagnostic. Comparer avec les bornes de Hoeffding existantes, les nuls globaux des PR #75/#97/#100 et des jeux biologiques publics *avant* toute décision d'inférence confirmatoire.

## Audit de déploiement

Le site MultiOmics publié au commit \`a93b1ce4\` a réussi son [déploiement GitHub Pages](https://github.com/rberrah/rberrah.github.io/actions/runs/38705569059) et le [smoke test public](https://github.com/rberrah/rberrah.github.io/actions/runs/38707124691) dès le premier passage : exploration mobile et bureau, import LC-MS AUC synthétique 25 features × 12 échantillons, ZIP complet et empreintes SHA-256, aucune requête HTTP(S) tierce après import. Ces succès testent l'interface et les exports, **pas** la calibration statistique.

L'évaluation indépendante d'accessibilité par des chercheurs débutants reste à organiser ; les tests automatiques ne remplacent pas cette étape.
