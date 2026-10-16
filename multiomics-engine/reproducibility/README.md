# MultiOmics — preuve exécutable de reproductibilité R (RNA-seq R1)

**Date : 13 octobre 2026. Statut : en validation CI, pas de certification
statistique ni clinique.** Ce chantier indépendant ne modifie ni l'analyseur
web, ni les formules R de production, ni les travaux de l'autre agent sur la
calibration à 400 jeux (PR #52) et le passage des contrastes (PR #54).

## Objectif et portée

L'environnement de référence doit pouvoir être **capturé puis reconstruit**
ailleurs. Le code réellement testé doit être relié à ses versions R, Bioconductor
et packages, à ses entrées, ses sorties et aux éventuelles erreurs.

La CI de `.github/workflows/multiomics-r-frozen-replay.yml` comporte deux
jobs distincts :

1. **Capture :** conteneur officiel Bioconductor 3.22 / R 4.5.2 fixé à son
   digest d'image `sha256:e34177896c1579f66d48b756e4229b58237499183a089da084af72d7aba1d40a`
   (pas un tag `latest`) ; installation de véritables méthodes DESeq2,
   edgeR et limma-voom ; analyse d'un seul jeu de comptages RNA-seq
   **totalement synthétique**, semence fixe 20261013, 24 échantillons et
   160 variables. Entrée, sortie et packages sont archivés.
2. **Restauration indépendante :** autre job, conteneur neuf de **même
   image exacte**, téléchargement des preuves, restauration depuis
   `renv.lock` dans une bibliothèque privée puis **réexécution** sur le
   même jeu synthétique. Le test refuse toute dérive des versions requises,
   des noms de gènes ou des valeurs de sortie (tolérance 1e-9).

Une exécution réussie démontre uniquement un **rejeu logiciel numérique
sur ce test**. Ce n'est ni une mesure de contrôle de FDR, ni une validation
biologique, ni une garantie que toutes les méthodes de l'outil, notamment
MOFA2, DIABLO, MS ou la survie, peuvent déjà être reconstruits.

## Artefacts disponibles après une CI effectivement réussie

- `synthetic_input.rds` — 24 sujets **virtuels**, gènes et vérité artificielle ;
- `reference_results.rds` — colonnes d'effet, p-valeur et q-valeur des
  trois implémentations R, triées par identifiant ;
- `runtime_manifest.json` — version de R/Bioconductor, versions
  individuelles des packages, limites, statut non certifié ;
- `renv.lock` — versionnement des packages directement nécessaires et
  de leurs dépendances transitives, créé **après installation** ;
- `sessionInfo.txt`, `source-commit.txt`, `SHA256SUMS` — versions,
  commit analysé et empreintes des fichiers capturés ;
- `replay-result.txt` — présent seulement en cas de réussite du rejeu.

Les artefacts sont conservés 30 jours par GitHub Actions : les télécharger
et les déposer dans une archive scientifique permanente avant toute
publication. Le ZIP source du logiciel et son propre commit, produits par
la PR #71, sont une preuve complémentaire.

## Limites et critères de refus

- **Ne pas confondre** image de base fixée et versions des packages. C'est
  précisément le `renv.lock`, suivi d'une restauration distincte, qui
  apporte la vérification supplémentaire.
- Même si les deux jobs réussissent, seules DESeq2, edgeR et limma-voom
  sont exercées ici, avec un jeu de données et une configuration fixés.
- Les contrôles des identifiants, des comptages entiers, du design et
  du contraste sont exécutés avant la régression ; une entrée fractionnaire
  intentionnellement invalide doit être refusée.
- Les versions de bases biologiques externes ne sont pas couvertes par
  ce jeu synthétique. Le contenu d'un dépôt tiers peut changer sans
  changer la version de l'algorithme R ; figer aussi les bases et les
  annotations dans une étude biologique réelle.
- Une restauration qui échoue, un package manquant ou une différence
  de résultat doit bloquer la déclaration « reproductible ».
- Les comparaisons R1 sur **400 jeux préspécifiés** de la PR #52 sont
  nécessaires mais distinctes. Elles devront être analysées avec
  incertitude et les échecs de convergence, sans interpréter une suite
  de CI vertes comme une preuve de FDR contrôlée à 5 %.
- L'accessibilité universelle nécessite encore des sessions réelles
  débutants et des vérifications NVDA/VoiceOver. Les tests Playwright
  seuls n'en apportent pas la preuve.

## Exécution par un autre chercheur

Lancer via GitHub Actions, ou sur un poste doté de Docker, en montant le
dépôt puis en réutilisant exactement l'image référencée dans le workflow.
Ne pas mettre d'éléments patients dans la commande de démonstration.

La première étape installe seulement les paquets indispensables et lance :

```sh
Rscript multiomics-engine/reproducibility/rnaseq_replay.R capture tmp/omics-r-frozen
```

La deuxième étape **doit être un environnement neuf**, restauré depuis le
lockfile. Après cette restauration, lancer :

```sh
Rscript multiomics-engine/reproducibility/rnaseq_replay.R replay tmp/omics-r-frozen
```

La documentation institutionnelle de Bioconductor décrit les images
officielles : <https://bioconductor.org/help/docker/>. Le mécanisme
`renv::snapshot()` / `renv::restore()` est décrit dans la documentation
R : <https://search.r-project.org/CRAN/refmans/renv/html/snapshot.html>.

**Niveau de preuve visé pour cette PR :**
« exécution reproductible de trois moteurs R de référence sur un jeu de
données synthétique préétabli, avec restauration de l'environnement
informatique », sans inférence scientifique au-delà.
