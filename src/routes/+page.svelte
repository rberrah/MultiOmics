<script>
 import {base} from '$app/paths';
 import {language} from '$lib/stores/language';
 /** @param {string} fr @param {string} en */
 const t=(fr,en)=>$language==='en'?en:fr;
</script>
<svelte:head><title>MultiOmics — {t('Recherche reproductible','Reproducible research')}</title></svelte:head>
<section class="intro">
  <p class="tag">MultiOmics</p>
  <h1>{t('Explorer une ou plusieurs omiques','Explore single or multiple omics')}</h1>
  <p>{t('Transcriptomique, protéomique et métabolomique. Analyses exploratoires transparentes, résultats de référence R quand le moteur local est disponible, et rapports scientifiques reproductibles.', 'Transcriptomics, proteomics and metabolomics. Transparent exploratory analyses, native reference R results when a local engine is available, and reproducible research reports.')}</p>
  <div class="actions">
    <a class="primary" href={base+'/multiomics/tool/'}>{t('Ouvrir l’outil','Open analysis tool')}</a>
    <a class="secondary" href={base+'/multiomics/learn/'}>{t('Comprendre les méthodes','Learn the methods')}</a>
  </div>
  <p class="disclaimer">{t('Une CI réussie ne garantit ni une FDR universellement contrôlée ni une validité clinique.', 'Passing CI does not guarantee universally calibrated FDR or clinical validity.')}</p>
</section>
<style>
.intro{max-width:910px;padding:4rem 1.25rem;margin:auto}
.tag{font-weight:700;font-size:.85rem;letter-spacing:.1em;text-transform:uppercase;color:var(--accent-pk)}
h1{font-size:clamp(2rem,5vw,3.4rem);line-height:1.15}
.intro>p{max-width:740px}
.actions{display:flex;flex-wrap:wrap;gap:.75rem;margin:2rem 0}
.actions a{padding:.75rem 1rem;border-radius:5px;font-weight:650;text-decoration:none}
.primary{background:var(--accent-ai);color:#fff}.secondary{border:1px solid var(--border-strong);color:var(--text-primary)}
.disclaimer{font-size:.85rem;color:var(--text-secondary)}
</style>