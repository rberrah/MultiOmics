# MultiOmics — deux limites confirmées par les preuves du 15 octobre 2026

Cette note documente les **échecs scientifiques observés**, sans les renommer en succès par le statut technique GitHub Actions. Aucun résultat ci-dessous n'est une validation générale du logiciel ou un résultat provenant de patients de l'utilisateur.

## 1. Fausse découverte lorsque des mesures sont manquantes de façon MNAR

Source réelle : [référence de calibration à vérité mixte](https://github.com/rberrah/rberrah.github.io/actions/runs/38711975291), [artefact JSON](https://github.com/rberrah/rberrah.github.io/actions/runs/38711975291/artifacts/11961067004).

Dans la simulation préspécifiée **MNAR opposite tails**, la moyenne du FDP (FDR empirique sur jeux indépendants contenant signal et nul) est **0,1908** avec N=24, IC affiché [0,1412 ; 0,2400], et **0,5855** avec N=80, IC [0,5522 ; 0,6164]. Le rappel des effets plantés est respectivement 0,120 et 0,697. Les cinq autres mécanismes simulés produisent environ 0,017–0,034 sur les mêmes tailles. Ces intervalles reflètent la simulation et ne doivent pas être interprétés comme une démonstration d'échantillonnage MAR/MNAR identifiable avec les observations.

**Implications** : absence de différence de pourcentage de valeurs manquantes entre groupes ≠ garantie MCAR. Le vrai mécanisme MNAR n'est pas identifiable à partir des seules données observées ; seuils et tests de missingness sont des *alertes*, pas une correction statistique. Toute conclusion confirmatoire requiert sensibilité à plusieurs hypothèses de non-observation, mesure des raisons de missingness, analyse avant/après imputation et revue indépendante. Un q≤0,05 ou une bonne puissance apparente ne certifie pas la validité.

État déjà implémenté : contrôles explicites dans \`scientific-assurance.js\`, statut de revue/blocage pour missingness, contrôle sur les maxima par échantillon même lorsque la médiane est nulle. Ne pas prétendre que l'on détecte automatiquement le MNAR.

## 2. Discrimination externe forte ≠ calibration de probabilité acceptable

Source réelle : [PMx figé entraîné MAINZ (n=200) → NKI (n=337), transfert inter-plateformes](https://github.com/rberrah/rberrah.github.io/actions/runs/38711975134), [recalibration TRANSBIG (n=198) → NKI (n=337)](https://github.com/rberrah/rberrah.github.io/actions/runs/38711975322).

Les mesures de classification ER de la voie **primary_rank** montrent sur NKI **AUC = 0,9855**, mais **Brier = 0,2793**, Brier skill = **−0,4104** comparé au prédicteur de référence utilisé dans cet audit ; le Brier empirique d'une autre référence descriptive est **0,197999**. Attention : les références rapportées dans les scripts ne sont pas nécessairement identiques ; ne pas recomposer artificiellement un nouveau score à partir d'un chiffre de l'autre protocole.

Une recalibration de l'intercept entraînée exclusivement sur TRANSBIG (\`−0,254831\`) est appliquée sans réentraînement du modèle sur NKI. Le Brier NKI **empire de 0,2793 à 0,3021**, différence **+0,022805** [bootstrap 95 % : +0,020651 ; +0,025202]. La discrimination AUC reste inchangée, car une transformation monotone de l'intercept ne modifie pas les rangs. **Le mot PASS du workflow signifie que le protocole a été exécuté et ses garde-fous respectés, pas que la recalibration était efficace.**

L'analyse TRANSBIG indépendante sans transfert de plateforme a montré une AUC de 0,8954 et un Brier de 0,1009 [workflow correspondant](https://github.com/rberrah/rberrah.github.io/actions/runs/38711975602). Cela ne dispense pas de quantifier la calibration sur la population de destination.

**Implications** : présenter AUC et score de probabilité (Brier et log-loss) côte à côte. Ne pas interpréter une AUC excellente comme des probabilités justes, une intervention utile ou un biomarqueur externe validé. Exiger une calibration externe (courbe, intercept/pente, Brier par rapport à une base préspécifiée et IC adaptés), transportabilité laboratoire/plateforme et évaluation d'utilité clinique avant usage décisionnel.

## 3. Correction d'affichage et frontière des claims

Pour toute prédiction **binaire** qui dispose des prédictions des folds externes de la *nested cross-validation*, le moteur détermine désormais \`brier = mean((prediction_probability − observed_binary)^2)\`, calculé au niveau des sujets effectivement tenus à l'écart, sans utiliser les valeurs des plis de test lors de l'apprentissage. L'UI FR/EN et le dossier Methods affichent AUC, Brier et log-loss quand disponibles, avec un avertissement sur les cohortes externes. Il ne s'agit **pas** d'une pente de calibration ; l'évaluation du modèle sur les mêmes sujets en nested CV reste une validation *interne*.

Dans les scénarios réels/indépendants, les valeurs numériques du présent document proviennent des scripts et logs R/JS identifiés ci-dessus ; elles ne sont **pas** injectées artificiellement dans les résultats de l'utilisateur. Les colonnes, nombres de patients et prévalences diffèrent entre benchmarks. Ne jamais fusionner ou moyenner des jeux indépendants sans leur protocole d'harmonisation.

**Statut** : qualification logicielle et evidence négative documentées ; **aucune certification confirmatoire universelle ou clinique**.
