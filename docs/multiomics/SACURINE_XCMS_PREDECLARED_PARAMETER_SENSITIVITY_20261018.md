# Sacurine v3 — analyse de sensibilité des paramètres xcms sur des spectres réels

**Protocole déclaré le 18 octobre 2026 avant le nouveau calcul.** Ce test complète la
[qualification instrumentale des douze fichiers](SACURINE_ZENODO_V3_REAL_INSTRUMENT_FINDINGS_20261018.md).
Il ne remplace pas les résultats négatifs du contrôle de qualité analytique.

## Hypothèse et données strictement identiques

Une table LC-MS de 3 135 features avait été extraite par `xcms` d'un unique
jeu externe [Zenodo Sacurine v3](https://zenodo.org/records/10130758) :
6 urines réelles, 3 pools QC, 3 blancs mesurés. Le même ensemble de douze
fichiers bruts et la **même chronologie source non contiguë (27–185)** sont
utilisés dans les deux traitements. Le téléchargement réel est vérifié par
l'empreinte MD5 officielle et les SHA-256 enregistrés dans le manifeste
d'acquisition. Aucun chromatogramme individuel n'est ajouté à Git.

## Deux réglages fixés, aucune optimisation sur les résultats

| Paramètre | Analyse primaire déjà qualifiée | Perturbation de sensibilité |
| --- | ---: | ---: |
| Détection `centWave` ppm | 25 | **20** |
| `centWave` snthresh | 8 | **10** |
| Peakwidth | 10–60 secondes | Identique |
| Prefilter | 3 scans, 100 unités | Identique |
| Alignement Obiwarp | `binSize=0.6` | Identique |
| PeakDensity correspondance | `bw=8`, `minFraction=0.5`, `minSamples=1`, `binSize=0.025`, `ppm=15` | Identique |
| Quantification `into` | Oui | Identique |
| Blancs exclus de la **définition** des groupes de features | Oui | Oui ; restent quantifiés |
| Définition biologique, origine et ordre des injections | 6/3/3, lot `ne1`, source vérifiée | **Strictement identiques** |

La perturbation est **combinée** (ppm et S/N changent en même temps) :
elle permet d'observer la sensibilité à ce couple de réglages, **pas**
d'attribuer causalement une variation à l'un ou l'autre.

## Comparaison préspécifiée sans réutiliser les identifiants xcms

Les ID de features peuvent changer lorsque xcms regroupe les pics différemment.
La correspondance est donc construite à partir des coordonnées natives
`mzmed` et `rtmed`, **sans supposer que FT0001 décrit le même signal
entre deux runs**.

1. Tolérance de coordonnées fixée a priori : **20 ppm** sur m/z
   et **30 secondes** sur le temps de rétention.
2. Pour être appariée, une feature doit avoir **exactement une candidate
   dans chaque sens** (référence→variante et variante→référence). Une
   ambiguïté, l'absence de candidate ou une relation non réciproque
   compte comme **non retrouvée**, jamais comme succès.
3. Le nombre de features retrouvées est rapporté sur **deux dénominateurs**
   (toutes les features de référence et toutes celles de la variante),
   sans sélectionner a posteriori le dénominateur le plus favorable.
4. Sur les seules features réciproques, le protocole rapporte le nombre
   de RSD de QC évaluables dans **les deux** traitements, et la concordance
   du signal descriptif **RSD >30 %** ; il fait de même pour le ratio
   blanc/urine **>20 %**. Les `NA` ne sont jamais assimilés à zéro.
5. Les différences signées de m/z (ppm) et de RT (secondes) et
   **chaque motif de non-appariement** sont conservés ligne par ligne.
   Aucun filtrage automatique de features pour publication n'est opéré.

Une feature appariée selon ces deux coordonnées n'est **pas** une
identification de métabolite. Les adducts, isotopes, fragmentation,
résolution chromatographique et signaux voisins restent à étudier
indépendamment. Une stabilité de l'ID, ou un chevauchement brut des noms,
n'est pas une preuve de robustesse.

## Critères informatiques et interprétation

Les tests doivent refuser identifiants dupliqués, m/z ou RT non finis,
alignement d'échantillons différent, erreurs de source, associations
plusieurs-à-un et faux scores de QC sur valeurs manquantes.

Le workflow requiert :

- présence des douze originaux et de leurs empreintes ;
- exécution `xcms` des deux réglages sur **les mêmes données** ;
- intégrité numérique des exports et rôles des échantillons ;
- tableau complet référence→variante, dénombrement exact des candidats
  univoques/ambigus/perdus et comparaison des écrans QC/blanc ;
- statut final **`scientific_sensitivity_review_required`**, sans
  seuil ad hoc de « réussite scientifique ».

Un workflow vert prouve que **le protocole a été exécuté et audité**,
pas que sa stabilité est favorable. Le taux de features retrouvées et
les mesures QC ne doivent être publiés comme observations scientifiques
qu'après inspection de l'artefact de la véritable exécution.

## Limites résiduelles de Research v1

Ce test reste **dans la même chaîne xcms**. Une validation réellement
indépendante exigera une autre implémentation ou méthode de prétraitement
(par exemple un pipeline scientifique comparateur explicite), un appariement
chimique annoté et/ou des standards analytiques, et plusieurs acquisitions
indépendantes. Il n'évalue ni la validité des hypothèses biologiques, ni
l'inflation du FDR DESeq2, ni une relation individus Omics→PK.

Les versions natives R et la bibliothèque des 206 dépendances sont
vérifiées avant le nouveau calcul. Tous les résultats défavorables
sont conservés ; aucune modification du moteur de production ou de
l'interface novice ne résulte automatiquement de cette sensibilité.
