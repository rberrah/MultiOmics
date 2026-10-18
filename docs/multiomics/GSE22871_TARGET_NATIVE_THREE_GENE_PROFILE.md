# GSE22871 : profil natif descriptif à trois gènes (18 octobre 2026)

## Question et limites

Le modèle **Nutrimouse** original contient quatre gènes sélectionnés dans les
30 souris d'apprentissage : `CYP3A11, FAS, CYP4A14, CYP4A10`. Le pont
indépendant avec les annotations des 45 101 sondes Affymetrix GPL1261 de
**GSE22871** permet de retrouver sans ambiguïté les trois premiers, dont
**FAS = Fasn (RefSeq NM_007988)**. La sonde candidate de CYP4A10 répond aussi
à Cyp4a31 : elle **reste exclue**, sans remplacement ni imputation.

Cette procédure supplémentaire calcule, pour chaque souris GSE22871 :

1. Les intensités **RMA log2 originales** des trois sondes strictement
   identifiées, sans renormalisation.
2. Le rang percentile avec ex æquo au rang moyen **au sein de la même
   puce et parmi les 45 101 sondes réelles** :
   `(nombre d'intensités inférieures + (nombre d'égales - 1)/2)/45100`.
3. Leur moyenne arithmétique **descriptive**, déclarée
   `target_native_three_probe_mean_midrank`, sans ajustement de paramètres
   inter-sujets ni utilisation du génotype ou de l'exposition PFOS.

Cela donne un exemple chiffré de **30 profils individuels réels, 90 mesures
de sonde**. Il ne s'agit **pas** d'une projection de la signature Nutrimouse
gelée à quatre gènes, ni d'un score CYP validé, ni d'un indicateur PK/PD.

## Exécution contrôlée

Le workflow
[`geo-gse22871-external-validation.yml`](../../.github/workflows/geo-gse22871-external-validation.yml)
récupère et hache les données originales, contrôle indépendamment les
annotations NCBI / Bioconductor, gèle les quatre caractéristiques Nutrimouse,
puis lance :

```bash
python3 scripts/qualify_geo_native_three_gene_profile.py --self-test
python3 scripts/qualify_geo_native_three_gene_profile.py tmp/multiomics/real-external-GSE22871
Rscript multiomics-engine/tests/test_geo_native_three_gene_profile.R --self-test
Rscript multiomics-engine/tests/test_geo_native_three_gene_profile.R tmp/multiomics/real-external-GSE22871
```

Les contrôles refont **indépendamment** les 90 intensités et percentiles et
les 30 moyennes depuis le fichier GEO compressé original en R, sans appeler
les calculs Python. Le workflow archive les valeurs individuelles, les SHA-256,
le rapport et le test R. L'échec d'une correspondance de sonde, d'une
empreinte ou de la reproduction numérique fait échouer la qualification.

Artefacts :
- `external_native_three_gene_descriptive_profile/individual_native_three_gene_probe_measurements.csv`
- `external_native_three_gene_descriptive_profile/individual_native_three_gene_descriptive_index.csv`
- `external_native_three_gene_descriptive_profile/native_three_gene_audit.json`
- `external_native_three_gene_descriptive_profile/independent_R_native_three_gene_replay.txt`

## Limites à conserver dans les articles

Le rang par puce permet de décrire la position relative de la sonde, **pas**
de déduire la quantité absolue de protéine, l'activité CYP ou une concentration
plasmatique. Une combinaison arithmétique de trois sondes, choisies parce
qu'elles sont techniquement identifiables, n'est pas nécessairement une
signature biologique fidèle aux quatre caractéristiques initiales.

Les données GSE22871 servent ici de source de mesures externes au modèle
Nutrimouse, **pas** de validation interplateforme de ce modèle. Les souris
GSE22871 n'ont pas de mesures PK appariées dans ce pipeline. Ces profils
descriptifs ne doivent pas être sérialisés en
`individual_signatures.csv` avec le `fit_scope`
`training_frozen_test_projection`.

L'issue [#20](https://github.com/rberrah/MultiOmics/issues/20) doit
**rester ouverte** : il manque une véritable calibration indépendante
interplateforme, la résolution du quatrième gène (ou un nouveau modèle
entièrement spécifié) et une validation biologique externe indépendante,
puis la qualification omiques–concentrations pour la PK.
