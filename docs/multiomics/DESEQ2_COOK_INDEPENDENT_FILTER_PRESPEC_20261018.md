# Étude préspécifiée : filtres de Cook et filtrage indépendant DESeq2 au nul global

**18 octobre 2026, protocole avant lecture des résultats de cette campagne.**

## Contexte

La calibration DESeq2 de MultiOmics sur deux générateurs négatifs
binomiaux *synthétiques* a fait apparaître une forte sensibilité
à l'estimation des dispersions. La [contre-épreuve Gamma–Poisson
indépendante](INDEPENDENT_GAMMA_POISSON_NATIVE_FINDINGS_20261018.md)
montre que le comportement d'une tendance seule peut se dégrader
drastiquement quand les vraies dispersions s'écartent de la
relation moyenne–dispersion.

Une question distincte reste ouverte : quelle fraction de
la variation des faux appels peut être due, dans ces simulations,
au filtrage des valeurs aberrantes **Cook** et au filtrage
indépendant de Benjamini–Hochberg implémentés dans `results()` ?
Le [manuel DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
distingue ces étapes.

## Cohortes et quatre bras appariés sans refaire l'ajustement

- Trois plans NB source : équilibré, profondeur/dispersion
  importantes, déséquilibre des groupes/lot.
  **Nouvelles répétitions 6001–6032**, hors précédents jeux.
- Trois plans Gamma–Poisson hétérogènes : les 96 cohortes
  5001–5032, **réutilisées et identifiées explicitement**
  pour une autre étude de mécanisme (pas une nouvelle
  validation statistique indépendante).
- **192 études complètes**, chacune 400 gènes et 40 sujets
  synthétiques, coefficient de condition strictement nul.
- Une seule estimation native `DESeq2::DESeq` par étude,
  avec **`minReplicatesForReplace=Inf` pour tous**. Cela
  désactive le **remplacement** des comptes aberrants, qui est
  une étape distincte du filtrage des p-values dans `results()`.
  Les mêmes dispersions, même design `~ batch + condition`,
  même normalisation, mêmes gènes, même test Wald et même
  préfiltrage de comptes sont alors partagés par les quatre
  appels `results` suivants.

| Bras | `cooksCutoff` dans results | `independentFiltering` |
| --- | --- | --- |
| `cooks_on_if_on` | Valeur native par défaut | TRUE |
| `cooks_on_if_off` | Valeur native par défaut | FALSE |
| `cooks_off_if_on` | FALSE | TRUE |
| `cooks_off_if_off` | FALSE | FALSE |

**Ces bras ne correspondent pas à quatre ré-estimations
de modèle.** Ils comparent uniquement le post-traitement
de *la même estimation native*. La désactivation commune
du remplacement de valeurs aberrantes est une décision
diagnostique ; elle peut modifier les taux absolus par
rapport au comportement usuel du moteur de production.
Aucune méthode ni aucun seuil utilisateur n'est modifié.

## Hypothèses et sorties déterminées avant les résultats

- **Unité d'observation : l'étude entière**, et non les gènes.
  Un nul global implique FDR = FWER pour chaque stratégie.
- Pour chaque cellule, nombre d'études comportant au moins une
  fausse découverte BH `q≤0,05`, parmi exactement 32 études,
  avec un intervalle individuel Clopper–Pearson à 95 %.
  Aucun intervalle n'est présenté comme une preuve simultanée.
- Nombre de p-values encore finies après Cook ; nombre de
  q-values finies après le filtre indépendant ; différence
  d'effectifs entre les p-values brutes et filtrées ;
  nombre de faux appels par cohorte.
- Les p-values brutes doivent être **identiques** lorsque
  seul le filtrage indépendant varie. Activer Cook ne peut
  **pas ajouter** de p-values finies. La valeur de
  `alpha=.05` est imposée à `DESeq2::results()`.
- Dans chaque bras, toutes les q-values natives admissibles
  doivent être égales (tolérance `1e-8`) à un recalcul BH
  avec une implémentation indépendante qui ne réutilise
  pas `stats::p.adjust`.
- Toutes les **768 tables de résultats** attendues et les
  **192 ajustements natifs** sont exigés ; aucune étude
  manquante, graine copiée, changement de nombre de gènes
  ou pseudoréplication de gènes n'est admis.

La faible précision Monte-Carlo de 32 études par cellule
reste une limite substantielle. Les conséquences d'un
filtrage sur un générateur donné ne prouvent pas qu'il
s'agisse du mécanisme d'inflation du FDR dans des données
RNA-seq biologiques réelles. La réutilisation des cohortes
Gamma–Poisson interdit de présenter cette campagne comme
une validation externe totalement indépendante.

## Décision logicielle

Le workflow sous R 4.5.2 / Bioconductor 3.22 et son
environnement de 171 packages figés archives les résultats
par étude, la provenance et les limites. Le statut
**PASS** signifie que les calculs ont été correctement
exécutés et comparés, pas que le FDR nominal 5 % est
universellement vérifié.

L'[issue scientifique #44](https://github.com/rberrah/MultiOmics/issues/44)
reste ouverte. Seule une qualification scientifique
indépendante sur d'autres générateurs ou données pourrait
motiver une évolution méthodologique. Les seuils, l'interface
débutante, le moteur R de production et les autres dépôts
ne sont pas modifiés.
