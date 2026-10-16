# PMx Explain — audit concurrentiel et critères de mise en production scientifique
_Date du bilan : 9 octobre 2026. Analyse fonctionnelle et documentaire, et non test comparatif sur des cohortes identiques._

## Positionnement

PMx Explain vise une plateforme mono- et multi-omique utilisable par une personne
sans programmation : question biologique → vérification de la structure des
données → méthode de référence adaptée → résultats interprétables → audit de
reproductibilité. L'exécution dans le navigateur reste une **exploration** ;
les méthodes de référence ne valent pas validation automatique d'une étude.

## Analyse concurrentielle sourcée

| Solution | Points forts documentés | Implication pour PMx Explain |
|---|---|---|
| [Galaxy](https://galaxyproject.org/eu/) | Historique, workflows, versions, partage, nombreuses ressources de calcul. | Conserver un historique complet et exécutable, pas seulement des tableaux de résultats ; améliorer l'installation et l'échelle. |
| [MetaboAnalyst 6.0](https://www.metaboanalyst.ca/MetaboAnalyst/faces/ModuleView.xhtml) | Spectres LC-MS, annotation, analyse statistique, voies, biomarqueurs. | Travailler la chaîne instrumentale MS, le QC et les conventions d'identifiants. |
| [Workflow4Metabolomics](https://workflow4metabolomics.org/) | Pipelines LC-MS/GC-MS/RMN de prétraitement, statistiques et annotations fondés sur Galaxy. | Prouver la qualité des fichiers issus de l'instrument, pas seulement celle des matrices importées. |
| [OmicsAnalyst 2.0](https://www.omicsanalyst.ca/home.xhtml) | Mono- et multi-omiques, réseaux, statistiques et exploration causale. | Étendre l'exploration réseau, sans présenter les corrélations comme de la causalité. |
| [PaintOmics 4](https://github.com/ConesaLab/paintomics4) | Mappage multi-omique sur KEGG/Reactome/MapMan. | Améliorer la couverture, les ambiguïtés d'identifiants et les annotations. |
| [mixOmics](https://mixomics.org/) | Méthodes documentées, notamment DIABLO supervisé, tuning et validation. | Valider les étapes multiblocs, la stabilité des sélections et la performance hors échantillon. |
| [iDEP](https://github.com/iDEP-SDSU/idep) | Analyse RNA-seq guidée et enrichissement ; recommande de vérifier indépendamment avant publication. | Construire un parcours plus court pour l'analyse RNA-seq courante et expliquer les hypothèses. |

**Aucun classement chiffré** n'est revendiqué sans tests utilisateurs et données identiques.

## Critères avant l'étiquette « professionnel »

**P0 — scientifique et sécurité**
- Méthodes de référence réellement exécutées, versions enregistrées ; aucune inférence confirmatoire sans plan identifiable.
- Vérification indépendante de l'erreur statistique selon le type de mesure et le plan ; contrôles négatifs/publics non présélectionnés ; cas MNAR, batch et dérive instrumentale.
- Aucune fuite d'information dans la validation croisée ; validation externe pour toute prétention de biomarqueur.
- Traçabilité des originaux et du protocole, archive locale vérifiée cryptographiquement ; confidentialité et consentement explicites pour les services externes.

**P1 — reproductibilité et opérations**
- Réexécution à environnement R/Bioconductor figé (R/Python/package lock et versions des bases biologiques) ; tests d'une réouverture du projet.
- Provenance de la chaîne LC-MS brute (mzML/mzXML, annotations, blanks/pooled-QC, drift correction), non couverte par le projet JSON texte.
- Gestion maîtrisée des jeux volumineux, erreurs récupérables, diagnostics et quotas ; version du modèle et des méthodes dans les exports.
- Support de dossiers publiables : scripts R rejouables, figure vectorielle avec légende et choix statistiques préspécifiés.

**P2 — usage novice**
- Études d'utilisabilité réalisées avec utilisateurs externes sans assistance : import, mapping, choix d'un plan, analyse, interprétation et reprise d'un projet.
- Tests WCAG 2.2 AA (clair/sombre, tableaux, graphiques, navigation au clavier, contraste), sans considérer quelques tests Playwright comme une certification complète.
- Manuel français et anglais, jeux de démonstration vérifiés, temps de réalisation et erreurs fréquentes mesurés.

## Éléments effectivement implémentés

- Analyse mono-omique, multi-omique, exploratoire et régression selon le plan ; connecteur R local optionnel ou requis.
- CI déterministe, contrôles méthodologiques spécialisés (HC3/MNAR/Cox) et benchmarks publics avec limites documentées.
- Projet portable v1 : fichiers sources **textuels** et sorties, SHA-256 vérifiés à l'import, **aucune réexécution automatique**, aucune connexion R/Reactome/ChEBI à l'import.
- États « exploratoire », « blocage », « revue indépendante requise » sans certificat scientifique automatique.

## Limites majeures non résolues

Le projet portable n'est pas une image exécutable du laboratoire : il n'inclut
pas les binaires spectraux ni un environnement R/Bioconductor figé. Les
versions externes des bases de connaissances peuvent changer. Les modèles
longitudinaux, la survie à risques concurrents, le MNAR et l'analyse
multiblocs ont des limites d'identifiabilité et de validation propres.
Ces écarts interdisent une revendication de « validation universelle ».
