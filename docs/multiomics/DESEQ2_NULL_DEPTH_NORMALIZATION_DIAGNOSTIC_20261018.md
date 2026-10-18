# MultiOmics — diagnostic ciblé de l'inflation FDR DESeq2 : facteurs de profondeur simulés

**Préspécifié le 18 octobre 2026 avant nouvelle campagne scientifique.**

## Motivation

Une campagne sous le nul global négatif-binomial a produit, avec les
méthodes R réelles, des fréquences de fausses découvertes **DESeq2**
22/128, 23/128 et 18/128 (trois scénarios). Les bornes inférieures des
intervalles simultanés exacts dépassaient 5 % dans chacun de ces
générateurs. [Résultats vérifiés](RNASEQ_HIGH_REP_FDR_FINDINGS_20261018.md)
et [issue #44](https://github.com/rberrah/MultiOmics/issues/44).

Les q-values avaient été recalculées correctement par l'algorithme BH :
cela ne prouve pas la validité de la distribution des p-values du test
Wald. Une explication possible est l'estimation des facteurs de taille,
mais elle n'est **ni démontrée ni présumée vraie** par ce protocole.

## Intervention purement diagnostique, sans modifier l'outil

Dans le simulateur existant, chaque compte RNA-seq suit une binomiale
négative avec les **vrais facteurs de profondeur par échantillon connus**,
le batch et la condition définis avant les mesures. Le contraste de
condition est **exactement nul**. Le pipeline original refait son
estimation de taille par la méthode usuelle DESeq2.

Cette expérience applique, sur les **mêmes jeux simulés**, les deux bras :

| Bras | Origine des facteurs de taille | Identique |
| --- | --- | --- |
| `DESeq2_default` | Estimés normalement par DESeq2 | |
| `DESeq2_oracle_true_depth` | Vrais facteurs de profondeur générateurs, normalisés à moyenne géométrique 1, **inconnus en pratique** | |

Les deux bras gardent le même modèle
`~ batch + condition`, le même préfiltrage de gènes (≥10 comptes dans
≥20 % des échantillons, règle exacte de l'adaptateur existant), la
même estimation de dispersion et le même test de Wald,
le même contraste, le même filtre indépendant DESeq2
et Benjamini–Hochberg `q≤0,05`. **Aucun des deux bras n'est un
nouveau paramètre de production.**

## Horizon fixe et mesures

- Trois scénarios originaux : équilibré, profondeur/dispersion extrêmes,
  groupes 12/28 avec batch partiellement confondu.
- **32 nouvelles études indépendantes par scénario** (répétitions
  2001–2032, graines originales définies par la même formule,
  hors jeux 1–128 et 1001–1128 antérieurs).
- Chaque étude a 400 gènes, 40 échantillons et aucun vrai effet de
  condition. Deux vrais ajustements DESeq2 par étude :
  **96 études × 2 = 192 ajustements natifs**.
- Estimand primaire : fréquence des études contenant au moins une
  découverte BH (FDR = FWER **uniquement sous le nul global**).
  Intervalle binomial exact à 95 % par cellule méthode/scénario.
- Critères descriptifs de mécanisme : différence entre facteurs de
  taille estimés et connus (`RMSE(log(SF))`, biais médian en log),
  fraction des p-values brutes ≤0,05, proportion de q-values finies
  et différences de fausses découvertes sur les études appariées.
- L'intégrité BH est réexécutée indépendamment des packages depuis
  leurs p-values sélectionnées (les `padj=NA` de filtre indépendant
  ne sont pas converties en zéro).

**N=32 ne permet pas de certifier un FDR nominal de 5 %**. Une diminution
des faux appels avec l'oracle indiquerait qu'une erreur de normalisation
pourrait contribuer au résultat *dans ce générateur*, mais ne démontrerait
pas que DESeq2 serait mal normalisé sur des cohortes biologiques réelles.
Une inflation persistante malgré la taille connue orienterait plutôt vers
la dispersion, le test Wald, les estimations ou le générateur. Les
hypothèses doivent être discriminées dans un benchmark indépendant ensuite.

La réussite GitHub CI doit indiquer **seulement** qu'un diagnostic natif
avec contraintes, données et prétraitements documentés s'est terminé ;
toute inflation ou valeur défavorable doit être conservée dans
l'artefact. **Aucun changement automatique du moteur ni sélection
de méthode sur la base de ce test**.

## Reproductibilité

Scripts :

- `multiomics-engine/rnaseq_deseq2_null_oracle.R` — générateur
  explicite et contrôles de pairage.
- `multiomics-engine/validation/rnaseq_deseq2_null_normalization_oracle.R`
  — vrais ajustements R, export par étude.
- `multiomics-engine/tests/test_rnaseq_deseq2_null_oracle.R` — tests
  sans Bioconductor (graines, groupes, modèles, appariement).
- `.github/workflows/deseq2-null-normalization-oracle.yml` — R
  4.5.2 / Bioconductor 3.22, empreinte de 171 packages figés et
  artefacts per-study/manifest.

Le résultat est **une analyse de sensibilité mécanistique**
et non une validation universelle de DESeq2, du site ou d'un usage clinique.
