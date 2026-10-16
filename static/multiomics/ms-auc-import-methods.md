# Import des aires de pics de spectrométrie de masse (MS AUC)

La mention **AUC MS** désigne l'intégrale du pic chromatographique (peak
area), c'est-à-dire une mesure de réponse instrumentale. Elle ne désigne
**ni l'AUC pharmacocinétique**, ni une concentration plasmatique.

## Fichiers acceptés

Exporter un CSV ou TSV depuis le logiciel de traitement LC-MS/GC-MS.

**Format long** (une ligne par couple identifiant de signal × injection) :

```csv
feature_id,assay_id,auc,subject_id,sample_id,condition,sample_type,injection_order
CHEBI:24996,MS1,102400,S1,PREL1,control,biological,1
CHEBI:24996,MS2,208300,S2,PREL2,treated,biological,2
CHEBI:30031,MS1,182000,S1,PREL1,control,biological,1
CHEBI:30031,MS2,91000,S2,PREL2,treated,biological,2
```

Les synonymes de colonnes courants sont reconnus : `metabolite`, `compound`,
`analyte`, `sample_id`, `sample_name`, `peak_area`, `area`,
`integrated_area`, `AUC`. Les identifiants sont conservés tels que fournis :
un nom d'analyte ou un identifiant d'ion de masse inconnu ne devient jamais
un métabolite validé.

Le fichier peut contenir `subject_id`, `sample_id`, `condition`,
`sample_type` (biological / qc / blank), `injection_order`,
`batch`, `timepoint`, `technical_replicate`.
Si `condition` est présent, les métadonnées sont extraites automatiquement ;
sinon, il faut fournir un fichier d'échantillons distinct.

**Format large** : première colonne `feature_id`, colonnes suivantes =
identifiants exacts des injections. Il exige les métadonnées distinctes.
Les fichiers `.xlsx` propriétaires ne sont pas importés directement ; exporter
en CSV/TSV pour éviter de perdre des identifiants ou des règles de séparation
décimale implicites.

### Conservation du fichier source et des annotations

Le fichier d'origine reste conservé localement, sans modification des octets.
L'import crée une matrice quantitative **et** une table d'annotations techniques
téléchargeable si les colonnes m/z, RT, adduits, formules, scores, etc. existent.
Dans un fichier long, les annotations sont conservées par couple signal ×
injection : aucune mesure potentiellement discordante n'est fusionnée. Le rapport
inclut le nom, la taille et, si disponible, le SHA-256 du fichier original ainsi
que le SHA-256 de la matrice transformée. Il faut archiver les fichiers à côté
du rapport, car les empreintes ne constituent pas une copie des données.

Un pic MS ne correspond pas nécessairement à une molécule. L'import ne lui
attribue jamais une identification chimique non déclarée. Une correspondance
`feature_id → CHEBI` fournie par l'utilisateur n'est qu'une annotation
déclarée, utilisée pour les figures mais pas pour modifier les statistiques.

## Règles scientifiques

- Une aire de pic doit être numérique, finie et non négative. Un signal
  absent peut rester vide. **Zéro n'est pas synonyme de donnée manquante.**
- Un doublon `(feature_id, assay_id)` est bloquant : l'outil ne décide pas
  arbitrairement de sommer ou moyenner des signaux.
- Les QC/blank, s'ils sont identifiés, peuvent servir à la détection
  des contaminants (blanks), au CV/RSD des pooled-QC et à la correction
  de dérive temporelle. Ils ne sont pas traités comme des individus.
- Le type `peak_area` entraîne un log2 avec pseudo-comptage documenté
  et un centrage médian par injection. Les rapports entre groupes ne sont
  interprétables qu'après lecture des étapes de prétraitement.
- Les aires de pics de **molécules différentes** ne représentent pas
  directement des concentrations comparables, car la réponse varie avec
  l'ionisation, la récupération, les standards et l'appareil.
- Aucune calibration en concentration, identification de métabolite, ni
  correction de batch non paramétrée n'est déduite automatiquement.
- Une seule couche MS permet une comparaison différentielle des groupes,
  temps ou outcomes ; **pas** de co-régulation inter-omique ni d'ACP
  multibloc. Une exploration intégrative exige au moins deux couches.

## Figures / effets

Une matrice `normalized` prétraitée est considérée telle que fournie :
l'effet est une **différence sur cette échelle**, jamais un log2FC présumé.
Les effets log2 sont affichés comme log2FC ; les effets en unité fournie
conservent leur étiquette explicite. Les couleurs sont normalisées
séparément par couche pour une présentation lisible, **sans comparaison
de l'amplitude des couleurs entre couches**.
Les q-values ajustées et les effets sont repris de l'analyse, sans nouveau
test par la visualisation.

Fichier exemple : `ms_peak_areas_example.csv` (valeurs **synthétiques**).
