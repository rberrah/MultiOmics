(() => {
  if (!window.location.pathname.includes('/multiomics/tool')) return;

  const terms = {
    fr: new Map([
      ['Effet biologique', 'Différence observée'],
      ['Significativité', 'Résultat après correction statistique'],
      ['Relations inter-omiques', 'Relations entre les omiques'],
      ['Axe latent principal', 'Tendance principale commune'],
      ['Covariables', 'Facteurs d’ajustement'],
      ['Covariables à ajuster', 'Facteurs à prendre en compte'],
      ['Association avec l’outcome', 'Lien avec le critère étudié'],
      ['Prédiction hors échantillon', 'Prédiction chez de nouveaux sujets'],
      ['Moteur R de référence', 'Analyse R de référence'],
      ['QC par couche', 'Contrôle qualité par type de données'],
      ['PCA QC', 'Vue globale des échantillons (ACP)'],
      ['screening', 'exploratoire'],
      ['Fold ratio', 'Amplitude (rapport)'],
      ['q BH', 'Résultat corrigé (q)'],
      ['Pattern', 'Profil'],
      ['r reference', 'Corrélation · groupe de référence'],
      ['r comparison', 'Corrélation · groupe comparé'],
      ['Δr', 'Différence de corrélation (Δr)'],
      ['Combined FDR', 'Résultat corrigé combiné'],
      ['FDR univers assay', 'Résultat corrigé avec les variables mesurées'],
      ['Nombre de couches concordantes', 'Nombre d’omics qui soutiennent la même voie'],
      ['Loadings', 'Variables qui contribuent le plus'],
      ['loadings', 'contributions des variables'],
      ['Outcome / critère', 'Critère étudié'],
      ['Type d’outcome principal', 'Type de critère étudié'],
      ['Temps omique utilisé pour l’outcome', 'Temps des mesures moléculaires utilisé pour étudier le critère'],
      ['Covariables / facteurs d’ajustement', 'Facteurs d’ajustement'],
      ['Batch connu ?', 'Plusieurs séries techniques sont-elles présentes ?'],
      ['Batch', 'Série technique'],
      ['Detected from uploaded features:', 'Type d’identifiant détecté dans le fichier :'],
      ['Column', 'Colonne'],
      ['Meaning', 'Signification'],
      ['Example', 'Exemple'],
      ['Rule', 'Règle'],
      ['Physical biological specimen', 'Prélèvement biologique'],
      ['Technical measurement / run', 'Mesure technique / analyse'],
      ['Measured layer', 'Type de données mesuré'],
      ['Experimental group', 'Groupe expérimental'],
      ['Visit / experimental time', 'Visite / temps expérimental'],
      ['Technical batch', 'Série technique'],
      ['Repeated technical assay', 'Répétition technique'],
      ['No file selected', 'Aucun fichier sélectionné'],
      ['Required fields mapped', 'Colonnes indispensables reconnues'],
      ['Mapping incomplete', 'Correspondance des colonnes incomplète'],
      ['None detected from repeated sample + omic pairs.', 'Aucune répétition technique détectée pour un même prélèvement et une même omique.'],
      ['not loaded', 'non chargée'],
      ['matched', 'appariés']
    ]),
    en: new Map([
      ['Biological effect', 'Observed difference'],
      ['Statistical evidence', 'Result after multiple-testing correction'],
      ['Cross-omics relationships', 'Relationships between omics'],
      ['Main latent axis', 'Main shared trend'],
      ['Covariates', 'Adjustment factors'],
      ['Outcome association', 'Association with the studied endpoint'],
      ['Out-of-sample prediction', 'Prediction in unseen subjects'],
      ['Reference R engine', 'Reference R analysis'],
      ['PCA QC', 'Global sample view (PCA)'],
      ['screening', 'exploratory'],
      ['Fold ratio', 'Magnitude (ratio)'],
      ['q BH', 'Adjusted result (q)'],
      ['Pattern', 'Profile'],
      ['r reference', 'Correlation · reference group'],
      ['r comparison', 'Correlation · comparison group'],
      ['Δr', 'Correlation difference (Δr)'],
      ['Combined FDR', 'Combined adjusted result'],
      ['Assay-universe FDR', 'Adjusted result using measured variables'],
      ['Covariates / adjustment factors', 'Adjustment factors']
    ])
  };

  const frenchFragments = new Map([
    ['Same value across all visits/omics for one participant, animal or culture.', 'Même valeur pour toutes les visites et toutes les omiques d’un même participant, animal ou modèle biologique.'],
    ['Same value only when assays come from the same specimen.', 'Même valeur uniquement lorsque les analyses proviennent du même prélèvement.'],
    ['Must match the corresponding matrix row/column identifier exactly.', 'Doit correspondre exactement à l’identifiant utilisé dans la matrice concernée.'],
    ['Canonical values: transcriptomics, proteomics, metabolomics.', 'Valeurs attendues : transcriptomics, proteomics ou metabolomics.'],
    ['Use one consistent vocabulary across subjects.', 'Utiliser les mêmes noms de groupes pour tous les sujets.'],
    ['Required for longitudinal designs.', 'Nécessaire pour les études longitudinales.'],
    ['Keep assay-specific batches even when different omics use different batches.', 'Conserver la série technique propre à chaque omique, même si elles diffèrent entre omiques.'],
    ['Distinct assay_id, same sample_id + omic.', 'assay_id différent, mais même sample_id et même omique.'],
    [' columns detected · delimiter: ', ' colonnes détectées · séparateur : '],
    [' assay rows', ' lignes de mesures'],
    ['Declared ', 'Déclaré : '],
    ['; metadata contains ', ' ; le tableau des échantillons contient '],
    [' data columns', ' colonnes de données'],
    [' expected assay IDs matched.', ' identifiants de mesure attendus retrouvés.'],
    ['Missing from matrix: ', 'Absents de la matrice : '],
    ['Not declared in metadata: ', 'Non déclarés dans le tableau des échantillons : '],
    ['Protocol: ', 'Plan d’étude : ']
  ]);

  function language() {
    return document.documentElement.lang === 'en' ? 'en' : 'fr';
  }

  function replaceText(root = document.body) {
    if (!root) return;
    const currentLanguage = language();
    const map = terms[currentLanguage];
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    for (const node of nodes) {
      const raw = node.nodeValue || '';
      const trimmed = raw.trim();
      const replacement = map.get(trimmed);
      if (replacement) {
        const leading = raw.match(/^\s*/)?.[0] || '';
        const trailing = raw.match(/\s*$/)?.[0] || '';
        node.nodeValue = leading + replacement + trailing;
        continue;
      }
      if (currentLanguage !== 'fr') continue;
      let updated = raw;
      for (const [source, target] of frenchFragments) updated = updated.replaceAll(source, target);
      if (updated !== raw) node.nodeValue = updated;
    }
  }

  let queued = false;
  const refresh = () => {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => {
      queued = false;
      replaceText();
    });
  };

  new MutationObserver(refresh).observe(document.documentElement, {
    childList: true,
    subtree: true,
    characterData: true
  });
  document.addEventListener('DOMContentLoaded', refresh, { once: true });
  refresh();
})();
