<script>
  import { language } from '$lib/stores/language';
  const terms = [
    ['Omique / Omics','Étude globale d’une classe de molécules mesurées (ARN, protéines, métabolites).','System-wide measurement of a class of molecules (RNA, proteins, metabolites).'],
    ['Échantillon biologique / Specimen','Prélèvement physique identifié par sample_id, qui peut être analysé plusieurs fois.','Physical biospecimen identified by sample_id, potentially assayed more than once.'],
    ['Assay','Analyse d’un échantillon associée à une colonne de matrice et à un identifiant unique.','Measurement of a specimen corresponding to one matrix column and a unique assay identifier.'],
    ['Répétition technique / Technical replicate','Nouvelle mesure du même prélèvement : ce n’est pas un sujet biologique supplémentaire.','Repeated measurement of the same specimen, not an additional biological subject.'],
    ['Effet différentiel / Differential effect','Différence de groupe ou de condition estimée par un modèle ; ce n’est pas une valeur individuelle.','Group- or condition-level model effect, not an individual measurement.'],
    ['FDR / q-value','Contrôle des découvertes multiples dépendant des hypothèses et de la définition de la famille de tests.','Multiple-testing adjustment whose interpretation depends on assumptions and the test family.'],
    ['QC / Contrôle qualité','Vérification de la qualité des mesures, réplicats, manquants et effets de lot.','Assessment of measurement quality, replicates, missingness and batch effects.'],
    ['MOFA','Modèle de facteurs latents multi-omiques ; une projection individuelle externe nécessite un modèle gelé.','Multi-omics latent factor model; external individual projection requires a frozen fitted model.'],
    ['DIABLO','Méthode supervisée multi-blocs ; l’évaluation prédictive nécessite une validation sans fuite.','Supervised multiblock method; predictive claims require leakage-safe validation.'],
    ['PK / Pharmacocinétique','Évolution de l’exposition à un médicament dans le temps.','Drug exposure and concentration changes over time.'],
    ['PD / Pharmacodynamie','Relation entre exposition, cible biologique et effet.','Relationship between drug exposure, biological target and effect.'],
    ['Provenance','Traçabilité des fichiers, échelles, versions et analyses exécutées.','Traceability of files, data scales, software versions and executed analyses.']
  ];
</script>
<svelte:head>
  <title>{$language === 'en' ? 'Glossary — MultiOmics' : 'Glossaire — MultiOmics'}</title>
  <meta name="description" content={$language === 'en' ? 'Omics and PK/PD terminology explained for researchers.' : 'Définitions pédagogiques en omiques et PK/PD.'}/>
</svelte:head>
<section class="content">
  <p class="eyebrow">MultiOmics · {$language === 'en' ? 'Resources' : 'Ressources'}</p>
  <h1>{$language === 'en' ? 'Omics glossary' : 'Glossaire des omiques'}</h1>
  <p>{$language === 'en' ? 'Definitions are educational; consult the methods and original references before reporting a scientific result.' : 'Ces définitions sont pédagogiques ; consulter les méthodes et les références originales avant de rapporter des résultats scientifiques.'}</p>
  <dl>{#each terms as term}
    <div class="entry"><dt>{term[0]}</dt><dd>{$language === 'en' ? term[2] : term[1]}</dd></div>
  {/each}</dl>
</section>
<style>
  .content{padding-block:clamp(2rem,5vw,4rem);max-width:1000px}
  .content>p{max-width:72ch;color:var(--text-secondary)}
  dl{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem;margin:2rem 0}
  .entry{border:1px solid var(--border-subtle);border-radius:var(--radius);padding:1.2rem;background:var(--bg-secondary)}
  dt{font-weight:700;color:var(--text-primary)}
  dd{margin:.5rem 0 0;color:var(--text-secondary)}
  @media(max-width:680px){dl{grid-template-columns:1fr}}
</style>
