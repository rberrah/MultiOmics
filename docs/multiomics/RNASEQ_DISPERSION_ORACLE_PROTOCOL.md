# MultiOmics P0 — attribution du FDR par dispersion connue et référence de Wald

**Protocole fixé le 14 octobre 2026 avant la nouvelle simulation.** Aucun moteur de production, aucune q-value patient et aucun seuil du site ne sont modifiés.

## Question motivant cette analyse

Les simulations indépendantes #52/#75/#77 montrent, sous plusieurs modèles binomiaux négatifs sans effet de groupe, un risque de fausse découverte supérieur à 5 % ; l'expérience #97 estime respectivement 17,5 % et 17,5 % avec normalisation estimée/vraie (n=40), puis 10 % contre 7,5 % (n=80) sur 40 jeux par plan. La normalisation n'explique donc pas, isolément, l'anomalie. Il reste à tester l'estimation des dispersions et la référence de Wald, **sans sélectionner un vainqueur à posteriori**.

## Générateur et plan figés

400 gènes tous nuls, moyennes basales log-uniformes entre 45 et 300, dispersions vraies U(0,08 ; 0,32), facteur de bibliothèque log-normal (écart-type 0,22) et effet de lot log(1,35). Deux plans équilibrés n=40 et n=80 avec conditions indépendantes, batch équilibré et design identifiable \`~ batch + group\`. 40 nouvelles simulations par plan, avec graines non utilisées dans #97.

Toutes les méthodes utilisent **strictement les mêmes matrices et la même pré-sélection sans labels** : au moins 10 comptages dans 20 % des sujets. Quatre bras par jeu :

1. Facteurs de taille et dispersion estimés, Wald normal (référence).
2. Facteurs de taille vrais, dispersion estimée, Wald normal.
3. Facteurs de taille vrais, dispersion vraie par gène, Wald normal.
4. Facteurs de taille vrais, dispersion vraie par gène, Wald avec référence t.

Un total de 80 jeux × 4 bras = **320 ajustements DESeq2 réels**. La substitution de vrais paramètres, possible uniquement en simulation, ne doit jamais être présentée comme une solution applicable à des données de patients. Le remplacement automatique des valeurs aberrantes est neutralisé de manière cohérente entre les bras afin de ne pas introduire une variation supplémentaire.

## Critères prédéfinis et gestion des résultats négatifs

Principaux : proportion de simulations avec au moins une découverte fausse (FWER = FDR sous le nul global), intervalle binomial exact 95 %, différence appariée entre bras (bootstrap exploratoire 2 000 rééchantillonnages), proportion de p-valeurs brutes sous 0,05.

Secondaires : correction BH avec/sans filtrage indépendant pour **le même fit**, ratio dispersion ajustée/vraie, erreur de normalisation, nombre de tests et p/q non finis. CSV complets par simulation, versions R/DESeq2, graines, manifeste, diagnostics et erreurs conservés en artefacts ; si un fit manque, le workflow échoue et préserve les données partielles.

Une fréquence élevée même avec vraie dispersion peut indiquer que des hypothèses supplémentaires sont problématiques (Wald en petit échantillon, dépendance/ajustement, autres éléments du générateur), **mais ne permet pas d'en attribuer la cause unique**. Une fréquence basse avec des paramètres vrais n'autorise pas à remplacer le modèle en production. Les 40 simulations laissent des IC larges, l'extrapolation biologique est interdite.

## Exécution

\`PMX_RNASEQ_DISPERSION_ORACLE_REPS=40 Rscript multiomics-engine/validation/rnaseq_deseq2_oracle_dispersion.R tmp/dispersion-oracle\`

Artefacts : \`per_fit.csv\`, \`summary.csv\`, \`paired_differences.csv\`, \`manifest.json\` et \`errors.csv\` lorsque nécessaire. Statut clinique ou scientifique automatiquement certifié : **jamais**.
