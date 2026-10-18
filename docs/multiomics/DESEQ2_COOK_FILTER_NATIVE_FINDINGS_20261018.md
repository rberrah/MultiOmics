# MultiOmics — Cook et filtrage indépendant : résultats numériques DESeq2 natifs

**Étude exécutée le 18 octobre 2026**, sur des **données synthétiques
entièrement nulles**, dans le cadre de la qualification statistique
de MultiOmics. Cette étude n'est ni une validation externe sur patients,
ni un certificat général de contrôle du FDR.

## Sources et traçabilité

- [Protocole fixé avant le calcul](DESEQ2_COOK_INDEPENDENT_FILTER_PRESPEC_20261018.md).
- [PR #57 fusionnée](https://github.com/rberrah/MultiOmics/pull/57)
  dans `main` au commit `f43a8b8d822e9b2cab167a8cb781a544e52d1539`.
- [Exécution réelle R 4.5.2 / Bioconductor 3.22](https://github.com/rberrah/MultiOmics/actions/runs/39165660638)
  (192 vrais ajustements DESeq2, 768 appels `results()`,
  intégrité BH vérifiée dans les 768 tables).
- [Artefact complet par étude](https://github.com/rberrah/MultiOmics/actions/runs/39165660638/artifacts/12150732184),
  archivé temporairement par GitHub : 768 lignes `CSV` complètes,
  synthèse de 24 cellules et manifeste. Chaque tuple
  `(scenario, replicate, seed, arm)` a été vérifié unique, et
  chaque cellule comprend exactement 32 études complètes.

Les deux sources de simulation ne sont **pas** des cohortes réelles :
le premier générateur `rnbinom()` utilise 96 nouvelles études
(répétitions 6001–6032), et le second générateur Gamma–Poisson latent
reprend ses 96 études indépendantes des simulations précédentes
(répétitions 5001–5032). Cette réutilisation est déclarée et ne
permet pas de revendiquer deux nouvelles validations externes.

## Pourquoi les quatre bras sont strictement comparables

Un seul ajustement `DESeq2::DESeq` par étude, mêmes comptes
synthétiques, mêmes 40 sujets, mêmes 400 gènes d'origine,
même formule `~batch+condition`, même contraste,
mêmes dispersions et facteurs de taille. Seul l'appel
`DESeq2::results(alpha=.05)` varie :

| Variante | Filtre Cook | Filtrage indépendant BH |
|---|---|---|
| `cooks_on_if_on` | Par défaut | Activé |
| `cooks_on_if_off` | Par défaut | Désactivé |
| `cooks_off_if_on` | Désactivé | Activé |
| `cooks_off_if_off` | Désactivé | Désactivé |

**Contrainte importante :** le remplacement automatique des comptes
aberrants est désactivé **dans les quatre bras**
(`minReplicatesForReplace=Inf`). Nous isolons donc l'effet de
*filtrage* de Cook et le filtre indépendant ; ce n'est pas une
réexécution exacte du comportement DESeq2 de production lorsque le
remplacement est permis.

## Résultats réellement observés

Le tableau représente les **études avec au moins un faux appel BH**
(`q≤0,05`) parmi 32 études par scénario. Ces nombres ont
été **strictement identiques dans les quatre bras** du protocole :

| Scénario nul global | Études avec ≥1 faux appel, dans chacun des quatre bras |
|---|---:|
| NB classique, équilibré | **4/32 (12,50 %)** |
| NB classique, profondeur/dispersion hétérogènes | **3/32 (9,38 %)** |
| NB classique, groupes déséquilibrés/lot | **4/32 (12,50 %)** |
| Gamma–Poisson indépendant, équilibré | **2/32 (6,25 %)** |
| Gamma–Poisson indépendant, profondeur/dispersion élevées | **3/32 (9,38 %)** |
| Gamma–Poisson indépendant, groupes déséquilibrés/lot | **5/32 (15,63 %)** |

Au nul global, la FDR d'une étude est égale à sa FWER
(probabilité d'au moins une fausse découverte), **mais les IC
Clopper–Pearson individuels restent larges** avec 32 études :
aucun bras n'est certifié à 5 %.

Le filtrage indépendant n'a exclu **aucune p-value supplémentaire
finie** de l'ajustement BH parmi les 192 études :
sa mise en marche ou à l'arrêt n'a modifié ni le nombre de
q-values finies ni les faux appels au niveau de l'étude.
Le filtre Cook n'a écarté que **deux p-values supplémentaires**
au total, dans le scénario Gamma–Poisson à forte dispersion,
soit en moyenne **0,0625 feature par étude** sur 32 études.
Il n'a modifié aucun statut d'étude « au moins un faux appel ».

Les invariants suivants sont testés automatiquement : les
p-values d'un même modèle sont identiques quand seul le filtre
indépendant varie, activer Cook ne peut restaurer des p-values
précédemment exclues et chaque famille de q-values est
recalculée indépendamment par Benjamini–Hochberg.

## Interprétation et limites obligatoires

**Dans les mécanismes simulés testés**, ni le filtre de Cook
appliqué aux résultats ni le filtrage indépendant n'expliquent
l'augmentation antérieure des faux appels. L'axe scientifique
prioritaire reste donc la qualité des estimations de dispersions,
leur régularisation, et la calibration du test de Wald.

Cela ne signifie surtout **pas** que ces filtres seraient inutiles
ou qu'il faudrait les désactiver en recherche. Les 400 gènes
simulés sont suffisamment exprimés pour que le filtrage
indépendant soit essentiellement inactif, et les simulations
ne reproduisent ni les distributions complexes de valeurs
aberrantes ni les nombreux gènes à faible abondance de vraies
études RNA-seq.

L'expérience a désactivé le remplacement automatique des comptes :
son impact sur un jeu réel, avec des valeurs influentes, reste
à évaluer distinctement. Les méthodes de production DESeq2
de MultiOmics **ne sont pas modifiées**.

**Décision Research v1 :** qualification numérique et diagnostic
mécanistique documentés ; **validation scientifique
confirmatoire universelle non acquise**. L'[issue #44](https://github.com/rberrah/MultiOmics/issues/44)
reste ouverte.
