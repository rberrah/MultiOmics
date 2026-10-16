# Réseau métabolique guidé par les données — périmètre scientifique

Le panneau **Réseau ciblé** exploite la sortie de l'analyse multi-omics sans modifier
les calculs, la correction des tests multiples ou les estimations d'effets.
Il affiche les régions comportant au moins une variable mesurée avec un effet
quantitatif exploitable sur une échelle explicitement déclarée : `log2`,
`as_supplied` ou `transformed_unknown`. L'étiquette d'échelle est conservée ;
un effet natif n'est jamais transformé fictivement en log2 fold change.
Le voisinage schématique comprend 0, 1 ou 2 pas. Une ACP descriptive sans
effets différentiels n'est pas une mesure d'activation de voie.

## Dictionnaire et confiance

- **Identifiant vérifié** : égalité stricte avec un identifiant ChEBI primaire
  explicitement référencé dans `src/lib/multiomics/metabolic-network.js`.
  Un identifiant ChEBI inconnu n'est jamais traité comme un nom chimique.
- **Référence externe vérifiée** : correspondance univoque retournée par la
  couche existante de résolution stricte (ChEBI/UniChem/HMDB/KEGG), lorsque le
  résultat pointe vers un ChEBI explicitement répertorié dans le dictionnaire.
  Les résultats `ambiguous`, `unresolved` et `api_error` sont refusés.
- **Nom seul** : correspondance textuelle exacte et unique dans le petit
  dictionnaire local ; présentée comme provisoire (contour pointillé), non
  assimilée à une preuve d'identité chimique.
- **Non reconnu ou ambigu** : aucune coloration ni raccordement artificiel.

Les différentes formes d'un métabolite (zwitterion, anion, acide, isotope,
stéréoisomère) ont souvent **des identifiants ChEBI distincts**. Ne pas
les fusionner automatiquement. Même principe pour des lipides isomères,
des métabolites dont la structure est incomplète et les mélanges.

Le dictionnaire actuellement inclus n'est **pas** le dictionnaire complet de
ChEBI. Seuls les identifiants explicitement renseignés sont considérés exacts.
Les noms sans numéro ChEBI renseigné restent provisoires ; le module n'invente
ni les identifiants manquants ni les correspondances entre bases.

Sources de référence pour contrôler et élargir le dictionnaire :
- ChEBI : https://www.ebi.ac.uk/chebi/ et
  https://www.ebi.ac.uk/chebi/downloads
- Rhea (réactions et participants ChEBI, distinction des microspecies) :
  https://www.rhea-db.org/help/reaction-participant
  https://www.rhea-db.org/help/download
- Reactome (événements et participants biologiques) :
  https://reactome.org/dev/content-service/

## Topologie et sélection de régions

Le catalogue intègre des régions de glycolyse, Krebs, PPP, acides aminés,
tryptophane-kynurénine, sérotonine, glutathion, urée, purines, lipides,
cholestérol/acides biliaires, choline, métabolisme à un carbone et hème.

Les liens sont **des adjacences schématiques d'orientation** ;
ils **ne sont pas** une reconstruction métabolique stœchiométrique,
une réaction Rhea vérifiée individuellement ni la preuve d'une réaction
physiquement possible dans le compartiment biologique de l'échantillon.
Les annotations des gènes représentent seulement leur *contexte*
enzymatique : elles ne sont ni des mesures de flux, ni des activités
enzymatiques. Le réseau n'utilise pas les associations statistiques
pour inventer de nouveaux arcs.

En mode automatique, les voies sont triées selon le nombre de points
mesurés et une voie qui n'ajoute aucun nouveau signal biologique est
masquée pour limiter la répétition des métabolites partagés (exemple :
tryptophane présent à la fois dans deux voies).

## Améliorations futures sans perte de reproductibilité

Pour une couverture quasi globale : importer à date fixe le ChEBI FULL
(identités et xrefs avec contrôles d'ambiguïté), les réactions Rhea
(identifiants des participants, direction et compartiments uniquement
quand disponibles) et les sous-voies Reactome propres à l'espèce.
Versionner le résultat dans le dépôt plutôt que télécharger et inférer
les liens à chaque analyse ; figer les versions des trois sources dans
le manifeste reproductible. La visualisation d'un réseau complet devra
continuer à fonctionner en mode focalisé.
