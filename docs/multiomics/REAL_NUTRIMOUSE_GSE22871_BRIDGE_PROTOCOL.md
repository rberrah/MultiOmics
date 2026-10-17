# Pont de qualification Nutrimouse → GSE22871 : identités génétiques réelles, sans transfert d'échelle

**Statut scientifique : audit de faisabilité du mapping, pas validation externe d'une signature.**

## Question

Le modèle déjà gelé sur **30 souris Nutrimouse** (40 souris au total, 10 de
test) a sélectionné `Cyp3a11` et trois autres gènes **exclusivement selon leur
variance dans l'entraînement**. Ses scores sont des moyennes de quatre
valeurs centrées-réduites par les moyennes et écarts-types macroarray du
**même entraînement**. Ces paramètres ne sont **pas applicables** aux
intensités Affymetrix RMA log2 de GSE22871 sans calibration interplateforme.

L'objectif du présent audit est de déterminer, sur les sources originales,
si **ces quatre gènes précis** peuvent même être reconnus de manière
univoque dans le profil génique des 30 autres souris. Les conditions
biologiques (régime/génotype Nutrimouse ; dose PFOS/génotype GSE22871)
ne jouent aucun rôle dans la sélection des gènes ou des sondes.

## Références et protocole de reproduction

1. Récupérer le vrai `nutrimouse.rda` à la révision immuable mixOmics
   `ef3e760526623d9e91e945e4b50be48d764efc40` et rejouer le véritable
   programme R `test_real_nutrimouse_frozen_individual_score.R` sous
   R 4.5.2/Bioconductor 3.22, sans jamais réestimer sur GSE22871.
2. Contrôler les **SHA-256 des modèles source**, les 30 sujets d'apprentissage
   et les 10 sujets de test, l'absence de fuite des valeurs externes et les
   quatre noms de gènes figés. Résoudre chaque nom dans
   `org.Mm.eg.db` avec tous les identifiants Entrez (ne jamais accepter
   arbitrairement le premier en cas d'ambiguïté).
3. Reprendre les **45 101 sondes originales** de GSE22871/GPL1261 et les
   annotations NCBI originales comparées à `mouse4302.db`. Pour chaque
   gène source, exiger simultanément le même symbole dans les deux
   références, le même identifiant Entrez que dans `org.Mm.eg.db`, aucun
   multi-gène et exactement **une** sonde admissible. Les conflits et
   sondes multiples sont archivés et bloquent le mapping univoque ;
   aucune association PFOS ne sert à choisir une sonde « optimale ».
4. Archiver toute correspondance (y compris les échecs) sous
   `real_source_to_target_gene_bridge/locked_panel_to_external_probe_crosswalk.csv`,
   ainsi que le JSON avec les SHA, versions, fraction de correspondance,
   raisons de blocage et compte **zéro** de scores externes exportés.
5. Exécuter les tests de contrôle négatif : faux Entrez, références
   contradictoires, multi-sondes, symboles voisins et identifiants mal
   formés doivent empêcher la conversion silencieuse.

Les jeux diffèrent à la fois **par la technologie et par l'intervention
expérimentale**. L'identité d'une cible génétique ne prouve ni la même
spécificité de sonde, ni l'équivalence quantitative des signaux, ni
l'invariance des relations gene–phenotype.

## Conditions impératives avant une vraie projection hors cohorte

Une étape distincte doit pré-définir, sans regarder les résultats de PFOS,
une normalisation ou un modèle de transfert étalonné sur une **référence
indépendante appropriée** (mesures jumelées des plateformes, étalons
ou mesure biologique externe autorisant une calibration) ; vérifier
l'agrégation des sondes, les effets de plateforme, l'incertitude et la
stabilité/ablation du score. Le modèle, l'ensemble de référence et ses
paramètres devront être gelés et hachés avant d'appliquer la procédure aux
30 souris GSE22871.

Si les quatre gènes ne sont pas tous univoquement représentés, la
projection exacte du modèle gelé doit rester **bloquée**, sauf élaboration
d'un **nouveau** score explicitement identifié et validé séparément. Même
un mapping 4/4 ne permet **pas** la transformation directe
`(Affymetrix_RMA - moyenne_macroarray) / sd_macroarray`.

**Validation PK** : aucune des deux études ne contient de cinétiques
pharmacocinétiques par sujet. Toute amélioration PopPK/PD reste à tester
plus tard sur une cohorte réellement appariée omiques–concentrations.
L'issue [#20](https://github.com/rberrah/MultiOmics/issues/20) demeure
ouverte jusqu'à une démonstration adéquate.

Ce contrôle est autonome dans MultiOmics ; il ne touche pas les autres
dépôts et ne prétend pas que le contrat d'export Omics2PK suffit à établir
une validation biologique.


## Résolution de l'identité FAS par la source primaire — 17 octobre 2026

La colonne `FAS` dans le modèle original gelé pouvait désigner deux gènes
murins distincts selon un dictionnaire général : récepteur apoptotique
**Fas / Entrez 14102** ou synthase des acides gras
**Fasn / Entrez 14104**, dont l'alias historique est `FAS`.
La PR #37 a correctement bloqué une affectation **sur le seul symbole**.

**Une preuve primaire lève cette ambiguïté de nom**, sans regarder aucun
résultat de l'analyse : **Martin et al., Hepatology 2007, vol. 45, p. 772,
tableau 4**, utilisent pour la ligne **FAS** l'accession RefSeq
**`NM_007988`**, classée dans la **synthèse des acides gras**.
L'identifiant [NM_007988 est le transcript de Fasn, Entrez 14104](
https://www.ncbi.nlm.nih.gov/gene/14104), et non celui de Fas/CD95.
Sources : [article primaire, PDF, page 6, tableau 4](
https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/hep.21510) ;
[PubMed DOI 10.1002/hep.21510](
https://pubmed.ncbi.nlm.nih.gov/17326203/).

La nouvelle vérification impose simultanément :

- `FAS` demeure l'**étiquette et la valeur numérique d'origine** dans le
  modèle Nutrimouse ; aucun modèle ni score n'est réestimé.
- RefSeq **`NM_007988`** doit être retrouvé dans la version figée
  de `org.Mm.eg.db` avec l'unique triplet
  **`Fasn / Entrez 14104 / NM_007988`** ; la simple concordance de
  `Fas / 14102` est un **contrôle négatif**.
- La recherche de sondes sur GPL1261 utilise le symbole canonique
  `Fasn` et l'Entrez 14104, exige une annotation identique dans
  les deux sources officielles NCBI et Bioconductor, et n'accepte
  automatiquement que **une** sonde sans ambiguïté. Toute autre
  possibilité, y compris les sondes multiples, est archivée et rejetée.

**Ce qui reste incertain :** le tableau 4 établit l'identité du transcrit
**visé par l'auteur**, mais ne fournit pas, à lui seul, la séquence exacte
du clone nylon INRARRAY ; une validation d'hybridation des sondes ne peut
donc pas être déduite. De plus les signaux macroarray de Nutrimouse et
Affymetrix RMA log2 de GSE22871 n'ont **aucune calibration interplateforme**
gelée. Une correspondance correcte `FAS→Fasn` n'autorise donc pas
la projection du score à quatre gènes ni une conclusion CYP/PK.

### Deuxième obstacle indépendant — CYP4A10

Le candidat Affymetrix `1424853_s_at` est annoté simultanément
**Cyp4a10 (Entrez 13117)** et **Cyp4a31 (Entrez 666168)** dans
les deux références NCBI et Bioconductor. Même si cet identifiant
partage le gène d'intérêt, le signal n'est pas spécifique ; il reste
rejeté pour la projection **exacte** du modèle source. Le nombre
de gènes identifiés avec une sonde strictement unique reste à recalculer
après la correction FAS (CI de la PR de résolution primaire).

Le fichier `all_possible_source_gene_GEO_probe_matches_no_selection.csv`
conserve les alternatives sans les choisir selon le génotype,
l'exposition PFOS ou les valeurs d'expression. La définition du modèle
gelé et ses SHA ne changent pas. **Zéro score externe exporté.**
