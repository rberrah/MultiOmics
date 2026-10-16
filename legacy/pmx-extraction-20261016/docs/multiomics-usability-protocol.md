# Multiomics — protocole d'évaluation réelle débutant / accessibilité

**Statut : protocole prêt, jamais encore exécuté avec des participants.**
Les tests Playwright, les vérifications de contraste et les CI scientifiques ne
constituent pas une étude d'utilisabilité humaine, ni un certificat WCAG.

**Périmètre :** [outil public mono- et multi-omique](https://rberrah.github.io/pharmacometrie/multiomics/tool/).
Ne tester que les fonctions réellement accessibles à la date de session. Ne pas
faire saisir de données patients, d'identifiants directs, de secrets ou de
résultats cliniques ; utiliser les démos intégrées ou des données synthétiques.

## Échantillon formatif initial

Huit sessions individuelles (deux personnes par profil) :
- **N1–N2 :** débutants sans pratique de bioinformatique ni des métadonnées ;
- **S1–S2 :** étudiants en sciences biomédicales, pharmacie ou biostatistiques ;
- **B1–B2 :** biologistes expérimentateurs connaissant leurs mesures mais non le pipeline ;
- **X1–X2 :** utilisateurs avancés (statistique, bioinformatique ou enseignant).

Les conclusions de huit sessions sont *des observations qualitatives*, pas des
proportions représentatives de l'ensemble des utilisateurs. Les participants
ne doivent pas avoir travaillé sur les écrans testés ni recevoir les solutions.

## Script indépendant (dans cet ordre)

Le modérateur lit uniquement l'objectif, sans indiquer où cliquer ni définir
les termes. Noter écran / appareil, heure, première action, erreurs, aides
sollicitées, interprétation verbalisée, résultat produit et export.

| ID | Consigne donnée à l'utilisateur | Observation factuelle attendue |
| --- | --- | --- |
| M1 | « Sans fichier, lancez une analyse ARN uniquement. Dites ce qu'on a mesuré. » | Trouve la démo mono-omique, voit les résultats sans backend R obligatoire et n'invente pas de résultats cliniques |
| M2 | « Lancez maintenant une analyse combinant plusieurs types de mesures. Qu'est-ce qui change ? » | Retrouve la démo multi-omique et distingue une mesure par omique d'une intégration prouvée |
| M3 | « Préparez deux mesures du même sujet et une mesure d'un autre sujet avec les modèles fournis. » | Comprend `subject_id` (individu), `sample_id` (prélèvement), `omic` (type de mesure), `condition` (groupe), réplicats biologiques/techniques ; n'assimile pas le nom de fichier à un sujet |
| M4 | « Un prélèvement a été mal nommé ou n'est pas apparié. Identifiez la cause du message d'erreur. » | Trouve le diagnostic, identifie l'identifiant incriminé et explique comment le corriger sans ignorer le garde-fou |
| M5 | « Que signifient fold ratio, p et q BH ? Peut-on dire qu'un q de 0,06 est une preuve ? » | Consulte l'aide, décrit un rapport de niveaux et une FDR contrôlée **dans la famille testée** ; n'affirme pas de découverte q≤0,05 |
| M6 | « Interprétez un résultat associé au génotype sur deux omiques. Prouve-t-il que les deux étaient nécessaires ? » | Reconnaît qu'une omique seule pourrait porter le signal ; sait que la normalisation/les contrôles négatifs peuvent changer l'interprétation |
| M7 | « Sauvegardez de quoi reproduire exactement l'analyse et expliquer ses limites à un collègue. » | Retrouve les résultats, paramètres, identifiants de jeu, éventuelle provenance/version, avertissements et limites ; distingue export et validation scientifique |
| M8 | « Reprenez M1 en anglais ou en français (langue opposée à celle de départ), sans souris. » | Navigation par Tab/Entrée, absence de piège clavier, traduction intelligible des principaux champs et résultats |

**Ne pas attribuer un succès à une fonctionnalité absente.** Si un scénario
n'est pas supporté par le site, le noter comme *non implémenté* (problème
produit), non comme un échec du participant.

## Questions d'interprétation indépendantes

Avant et après les tâches (sans corriger les réponses pendant la session) :
1. « Une métadonnée, c'est quoi ici ? » (information sur un échantillon).
2. « Deux omiques sur le même patient garantissent-elles une meilleure prédiction ? » (non).
3. « q BH = 0,03 signifie-t-il 3 % de chances que ce gène soit faux ? » (non ; interprétation de FDR sur un ensemble de rejets, sous hypothèses).
4. « Une AUC élevée garantit-elle de bonnes probabilités individuelles ? » (non ; calibration, Brier et cohorte externe à vérifier).
5. « Les résultats de cette démonstration suffisent-ils à soigner un patient ? » (non).

Ces réponses servent à comprendre les incompréhensions, **jamais à fabriquer
un pourcentage de « maîtrise » statistiquement infondé**.

## Passes d'accessibilité manuelles

Effectuer M1, M3, M5 et M7 sur navigateur desktop au zoom **200 %**,
largeur CSS **320 px**, clavier seul, puis lecteur d'écran **NVDA** sur
Windows et **VoiceOver** sur un appareil pris en charge. Utiliser au moins un
véritable téléphone pour M1/M2.

Contrôler le focus visible, l'ordre de tabulation, l'absence de défilement
horizontal, les annonces d'erreurs et de chargement, la lecture intégrale
des tableaux/figures, les bulles « ? » sans coupure, les contrastes,
les équivalents textuels des graphiques, la navigation FR/EN et l'option
réduction du mouvement. Ne jamais déduire la réussite NVDA/VoiceOver d'un
test DOM automatisé.

## Compte rendu / règle de correction

Remplir [le registre d'observation](multiomics-usability-results.md) **après
chaque session**. Pour toute tâche : statut
`autonome` / `avec aide` / `non terminée` / `non implémentée`,
temps observé (pas estimé), citation courte anonymisée et ID GitHub
de correction. Conserver aussi les cas où l'utilisateur se trompe
scientifiquement malgré un parcours de clics réussi.

En priorité corriger **P0** si la méthode peut conduire à une
interprétation statistique fausse, si une erreur de sujet passe inaperçue,
si l'analyse ne s'exécute pas ou si un chemin clavier indispensable
est bloqué ; **P1** si les termes/étapes restent incompris sans aide ;
**P2** pour ergonomie secondaire. Refaire exactement les tâches concernées
avec de nouveaux utilisateurs après correction.

**Critère de communication :** « testé avec utilisateurs » uniquement
après des sessions datées et documentées ; « accessible » uniquement
avec résultats de technologies d'assistance réelles et limitations
signalées. **« scientifiquement confirmatoire » demande séparément une
validation méthodologique et biologique externe** ; cette étude UX ne la remplace pas.
