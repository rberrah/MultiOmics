# MultiOmics — résultats RNA-seq à haute répétition : inflation du FDR sous nul global

**Observation figée le 18 octobre 2026. À conserver même si des travaux ultérieurs expliquent ou corrigent le mécanisme.**

## Provenance vérifiable

- Protocole préspécifié : [RNASEQ_HIGH_REP_FDR_PRECISION_20261018.md](RNASEQ_HIGH_REP_FDR_PRECISION_20261018.md), fusionné par [PR #43](https://github.com/rberrah/MultiOmics/pull/43).
- Exécution des calculs scientifiques **réels** : [GitHub Actions #39143010669](https://github.com/rberrah/MultiOmics/actions/runs/39143010669), passée après les 1 728 fits R.
- [Artefact original de 17 fichiers, dont 576 jeux simulés, q-values auditées et rapports](https://github.com/rberrah/MultiOmics/actions/runs/39143010669/artifacts/12142955814) (rétention Actions 90 jours, à archiver durablement avant expiration).
- SHA des sources **réellement exécutées**, avant fusion : `796ba2c89cb65852016380b61615ace7aff50655`. Fusion squash PR #43 : `0845c88d957bc6bea91bd1e7e5877e0564c8ca14`.
- ZIP GitHub Actions : SHA-256 observé `ea6ee39580999bb94523d0d41772d497d16efefeeef05fd2712d5b63f0bb2fb4` ; empreinte du fichier compressé, non des CSV individuellement.
- R 4.5.2 / Bioconductor 3.22 dans le conteneur épinglé, 400 gènes, contrastes natifs `~ batch + condition` ; packages et graines exacts conservés dans l'artefact.
- Trois méthodes sur **les mêmes données simulées**, et non trois groupes d'études indépendants. Nouvelles répétitions 1001 à 1128 (nul global) et 1001 à 1064 (vérité mixte).

## Résultat primaire : nul global, BH q ≤ 0,05

Au nul global, tout gène détecté est un faux positif. La fréquence des études avec ≥1 découverte est donc exactement l'estimateur empirique du FDR (égal au FWER **dans ces scénarios**). Chaque cellule repose sur **128 jeux indépendants**, avec ajustements natifs DESeq2, edgeR-QL et limma-voom. Intervalles de Clopper–Pearson bilatéraux simultanés à 95 %, avec correction de Bonferroni sur les **neuf** cellules méthode × scénario.

| Générateur nul | Méthode R native | Jeux avec ≥1 faux appel / 128 | FDR observé | IC simultané 95 % | Statut vs 5 % |
|---|---|---:|---:|---:|---|
| Groupes équilibrés, profondeur modérée, lot équilibré | DESeq2 | **22** | **17,19 %** | **[9,11 ; 28,16] %** | **Inflation démontrée pour ce générateur** |
| Idem | edgeR-QL | 14 | 10,94 % | [4,68 ; 20,67] % | Indéterminé |
| Idem | limma-voom | 9 | 7,03 % | [2,28 ; 15,64] % | Indéterminé |
| Groupes équilibrés, forte dispersion, profondeur hétérogène | DESeq2 | **23** | **17,97 %** | **[9,70 ; 29,06] %** | **Inflation démontrée pour ce générateur** |
| Idem | edgeR-QL | 4 | 3,12 % | [0,45 ; 10,10] % | Indéterminé |
| Idem | limma-voom | 2 | 1,56 % | [0,06 ; 7,57] % | Indéterminé |
| Groupes 12/28, lot partiellement confondu mais identifiable | DESeq2 | **18** | **14,06 %** | **[6,82 ; 24,48] %** | **Inflation démontrée pour ce générateur** |
| Idem | edgeR-QL | 6 | 4,69 % | [1,08 ; 12,41] % | Indéterminé |
| Idem | limma-voom | 9 | 7,03 % | [2,28 ; 15,64] % | Indéterminé |

La borne inférieure **simultanée** pour DESeq2 dépasse 5 % dans les trois plans. C'est un résultat statistique négatif substantiel, pas une erreur d'intégrité du calcul BH : les q-values renvoyées par les packages ont fait l'objet d'un rejeu algorithmique indépendant réussi. **La distribution des p-values sous le nul en amont de BH n'est pas certifiée.**

## Résultat secondaire : vérité mixte

Dans des simulations de 64 études indépendantes par taux de gènes réellement modifiés (effet planté ±1,5 log2FC), FDR = moyenne des fractions de découvertes fausses `V/max(R,1)` calculée **par étude entière**, et non par pseudo-réplication des gènes.

| Fraction de vrais signaux | DESeq2 — FDP moyen | edgeR-QL — FDP moyen | limma-voom — FDP moyen |
|---|---:|---:|---:|
| 1 % | 9,17 % | 3,48 % | 4,84 % |
| 5 % | 6,82 % | 3,95 % | 3,92 % |
| 20 % | 5,06 % | 3,67 % | 3,78 % |

Les bornes supérieures simultanées Hoeffding–Bonferroni sur ces 64 études sont **très conservatrices** (environ 25–31 %) ; ces résultats n'établissent donc pas une garantie formelle de FDR nominal par méthode. Les puissances quasi parfaites sont attendues ici avec des effets plantés très grands et des intensités élevées : elles ne prouvent pas une sensibilité équivalente sur cohortes biologiques réelles.

## Implications pour l'outil et le protocole scientifique

1. **Conserver le moteur de production inchangé** : aucune modification des p-values, de BH, du seuil de 0,05, du modèle DESeq2 ou sélection automatique de edgeR-QL/limma d'après ce benchmark.
2. **Maintenir la revue statistique indépendante obligatoire** pour toute analyse confirmatoire, en particulier lorsque l'utilisateur croit que `q≤0,05` prouve automatiquement un taux de faux appels contrôlé à 5 % dans son étude.
3. Ne pas confondre **validité du code de calcul BH**, **calibration du test Wald**, **hypothèses du générateur**, **validité biologique externe** et **usage clinique**.
4. Diagnostiquer séparément les facteurs de taille, l'estimation des dispersions, les statistiques/test de Wald, l'influence éventuelle du filtrage indépendant, la distribution empirique des p-values nulles et les hypothèses de génération des comptes, par un protocole **préspécifié** et une implémentation externe indépendante.
5. Étendre ensuite à des générateurs hors famille binomiale négative, à des échantillons biologiques appariés et à de vraies cohortes multi-instruments. Aucun de ces résultats ne permet une certification scientifique universelle ou clinique.
6. Reporter le travail d'attribution causale et de reproductibilité indépendante dans [issue #44](https://github.com/rberrah/MultiOmics/issues/44). Cette issue reste **ouverte** ; les seuls tests logiciels PASS ne la résolvent pas.

### Rejouer la campagne à l'identique

Depuis le commit source précisé plus haut, lancer le workflow
`.github/workflows/rnaseq-monte-carlo-precision.yml` et archiver son
`RNASEQ_MONTE_CARLO_PRECISION_REPORT.md`, les deux CSV de synthèse, les
deux CSV des études individuelles, les deux audits BH natifs ainsi que
`rnaseq_precision_run_manifest.txt`. Aucun résultat de patient n'a été
utilisé dans cette campagne.

**Conclusion permise :** inflation statistiquement démontrée pour DESeq2
dans ces **trois générateurs nuls préspécifiés**. **Conclusion non
permise :** « DESeq2 est universellement mal calibré » ou « edgeR/limma
sont universellement meilleurs ».
