<script>
  import {language} from '#lib/stores/language.js';
  const refs = [
    {name:'DESeq2 — Love, Huber & Anders (2014)',journal:'Genome Biology 15:550',url:'https://doi.org/10.1186/s13059-014-0550-8',scope:'RNA-seq differential expression'},
    {name:'limma — Ritchie et al. (2015)',journal:'Nucleic Acids Research 43:e47',url:'https://doi.org/10.1093/nar/gkv007',scope:'Linear modelling and empirical Bayes'},
    {name:'fgsea — Korotkevich et al. (2021)',journal:'bioRxiv preprint, 060012',url:'https://doi.org/10.1101/060012',scope:'Fast preranked gene set enrichment'},
    {name:'MOFA2 — Argelaguet et al. (2020)',journal:'Genome Biology 21:111',url:'https://doi.org/10.1186/s13059-020-02015-1',scope:'Multi-omics factor analysis'},
    {name:'DIABLO — Singh et al. (2019)',journal:'Bioinformatics 35:3055–3062',url:'https://doi.org/10.1093/bioinformatics/bty1054',scope:'Multi-block omics integration'},
    {name:'xcms — Smith et al. (2006)',journal:'Analytical Chemistry 78:779–787',url:'https://doi.org/10.1021/ac051437y',scope:'Mass-spectrometry feature detection'},
    {name:'Reactome — official documentation',journal:'Curated biological pathways',url:'https://reactome.org/',scope:'Pathway provenance and mappings'},
    {name:'MultiOmics — source and citation file',journal:'Versioned software and qualification workflows',url:'https://github.com/rberrah/MultiOmics/blob/main/CITATION.cff',scope:'Cite the actual software version'}
  ];
</script>
<svelte:head>
  <title>{$language === 'en' ? 'Scientific references — MultiOmics' : 'Références scientifiques — MultiOmics'}</title>
  <meta name="description" content={$language === 'en' ? 'Scientific references for the methods behind MultiOmics.' : 'Références scientifiques des méthodes utilisées par MultiOmics.'}/>
</svelte:head>
<section class="references">
  <p class="eyebrow">MultiOmics · {$language === 'en' ? 'Resources' : 'Ressources'}</p>
  <h1>{$language === 'en' ? 'Scientific references' : 'Références scientifiques'}</h1>
  <p>{$language === 'en' ? 'These are methodological references, not proof that a specific backend was executed or independently validated on your dataset. Consult the exported methods manifest for the computation actually used.' : 'Ces articles documentent les méthodes, sans prouver leur exécution ou leur validation indépendante sur vos données. Vérifiez le manifeste de méthodes exporté pour connaître les calculs réellement exécutés.'}</p>
  <ol>{#each refs as item}
    <li><strong>{item.name}</strong> — {item.journal}<p>{item.scope}</p>
      <a href={item.url} target="_blank" rel="noopener noreferrer">{item.url} ↗</a>
    </li>
  {/each}</ol>
</section>
<style>
  .references{padding-block:clamp(2rem,5vw,4rem);max-width:940px}
  .references>p{max-width:76ch;color:var(--text-secondary)}
  ol{padding-left:1.5rem;margin-top:2rem}
  li{padding:1rem .6rem;border-bottom:1px solid var(--border-subtle)}
  li p{margin:.3rem 0;color:var(--text-secondary)}
  li a{overflow-wrap:anywhere}
</style>
