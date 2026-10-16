(() => {
  if (!window.location.pathname.includes('/multiomics/tool')) return;

  const language = () => document.documentElement.lang === 'en' ? 'en' : 'fr';

  function renderFileGuide() {
    const firstInput = document.querySelector('input[type="file"][accept*=".csv"]');
    if (!(firstInput instanceof HTMLInputElement)) return;
    const panel = firstInput.closest('section.panel');
    if (!panel) return;

    let guide = document.querySelector('[data-testid="multiomics-file-guide"]');
    if (!guide) {
      guide = document.createElement('aside');
      guide.className = 'licence-file-guide';
      guide.dataset.testid = 'multiomics-file-guide';
      const head = panel.querySelector('.section-head');
      if (head) head.insertAdjacentElement('afterend', guide);
      else panel.insertAdjacentElement('afterbegin', guide);
    }

    const lang = language();
    if (guide.dataset.language === lang) return;

    guide.innerHTML = lang === 'en' ? `
      <details class="licence-file-extra">
        <summary>How should I prepare my files?</summary>
      <div class="licence-file-head"><span>Before choosing files</span><h3>What do I need to upload?</h3><b>1–3 omics matrices + a sample sheet (or guided creation)</b></div>
      <p>Start with one omic (RNA, proteins or metabolites). Add a second different omic for multi-omics; a third layer is optional. Do not duplicate the same layer to make it multi-omics.</p>
      <p>For exploration with one omic only, you can upload a single matrix and use the on-page helper to create the sample sheet <strong>after confirming that each column is an independent biological sample</strong>. This is not suitable for group or longitudinal inference: provide actual sample metadata for those questions.</p>
      <div class="licence-file-grid">
        <article><span>File 1</span><strong>Sample sheet</strong><p>One row describes one technical measurement: who was measured, which specimen, which omic, group and time point.</p><code>subject_id · sample_id · assay_id · omic</code></article>
        <article><span>Files 2–4 (choose only those measured)</span><strong>Omics matrices</strong><p>One matrix per omic. Rows are genes, proteins or metabolites; columns are measurements.</p><code>feature_id | RNA001 | RNA002 | …</code></article>
      </div>
      <div class="licence-id-chain"><strong>How the files connect</strong><div><span>P001<br><small>subject</small></span><i>→</i><span>P001_T0<br><small>specimen</small></span><i>→</i><span>RNA001<br><small>measurement</small></span><i>→</i><span>RNA matrix column<br><small>same assay_id</small></span></div></div>
      <details><summary>Single-omic example (RNA only)</summary><pre>P001  P001_T0  RNA001  transcriptomics
P002  P002_T0  RNA002  transcriptomics</pre><p>Upload the sample sheet and the RNA matrix only; each matrix column must match an assay_id. These are two different biological subjects, not two technical replicates.</p></details>
      <details><summary>Example with the same specimen measured by three omics</summary><pre>P001  P001_T0  RNA001   transcriptomics
P001  P001_T0  PROT001  proteomics
P001  P001_T0  MET001   metabolomics</pre><p>The subject and specimen IDs stay the same. The assay ID changes because these are three different measurements.</p></details>
      </details>
    ` : `
      <details class="licence-file-extra">
        <summary>Comment préparer mes fichiers ?</summary>
      <div class="licence-file-head"><span>Avant de choisir les fichiers</span><h3>Que faut-il charger ?</h3><b>1 à 3 matrices omiques + tableau des échantillons (ou assistant)</b></div>
      <p>Commencez avec une seule omique (ARN, protéines ou métabolites). Deux omiques différentes permettent une analyse multi-omique ; une troisième est facultative. Deux fichiers de la même omique ne constituent pas une analyse multi-omique.</p>
      <p>Pour explorer une seule omique, vous pouvez charger uniquement la matrice et utiliser l’assistant pour créer le tableau <strong>après avoir confirmé que chaque colonne correspond à un échantillon biologique indépendant</strong>. Pour comparer des groupes ou des temps, fournissez vos vraies informations d’échantillons : l’assistant ne doit pas les inventer.</p>
      <div class="licence-file-grid">
        <article><span>Fichier 1</span><strong>Tableau des échantillons</strong><p>Une ligne décrit une mesure technique : qui a été mesuré, quel prélèvement, quelle omique, quel groupe et quel temps.</p><code>subject_id · sample_id · assay_id · omic</code></article>
        <article><span>Fichiers 2–4 (uniquement les omiques mesurées)</span><strong>Matrices omiques</strong><p>Une matrice par omique. Les lignes sont les gènes, protéines ou métabolites ; les colonnes sont les mesures.</p><code>feature_id | RNA001 | RNA002 | …</code></article>
      </div>
      <div class="licence-id-chain"><strong>Comment les fichiers se relient</strong><div><span>P001<br><small>sujet</small></span><i>→</i><span>P001_T0<br><small>prélèvement</small></span><i>→</i><span>RNA001<br><small>mesure</small></span><i>→</i><span>colonne matrice RNA<br><small>même assay_id</small></span></div></div>
      <details><summary>Exemple mono-omique (ARN uniquement)</summary><pre>P001  P001_T0  RNA001  transcriptomics
P002  P002_T0  RNA002  transcriptomics</pre><p>Chargez le tableau des échantillons et la matrice ARN seulement ; chaque colonne de matrice doit correspondre à un assay_id. Ce sont deux sujets biologiques distincts, et non deux réplicats techniques.</p></details>
      <details><summary>Exemple : un même prélèvement mesuré par trois omiques</summary><pre>P001  P001_T0  RNA001   transcriptomics
P001  P001_T0  PROT001  proteomics
P001  P001_T0  MET001   metabolomics</pre><p>Le sujet et le prélèvement restent identiques. L’identifiant de mesure change car il s’agit de trois dosages différents.</p></details>
      </details>
    `;
    guide.dataset.language = lang;
  }

  function installStyle() {
    if (document.getElementById('multiomics-file-guide-style')) return;
    const style = document.createElement('style');
    style.id = 'multiomics-file-guide-style';
    style.textContent = `
      .licence-file-extra > summary{cursor:pointer;font-weight:700}
      .licence-file-guide{margin:14px 0 20px;padding:16px;border:1px solid var(--border-strong,#c9c9c9);border-radius:12px;background:var(--bg-secondary,#f7f7f7)}
      .licence-file-head{display:grid;grid-template-columns:1fr auto;gap:3px 14px;align-items:start}.licence-file-head>span{grid-column:1/-1;font-size:.74rem;text-transform:uppercase;letter-spacing:.08em;color:var(--text-secondary,#666)}.licence-file-head h3{margin:2px 0;font-size:1.08rem}.licence-file-head>b{border:1px solid var(--border-subtle,#ccc);border-radius:999px;padding:5px 9px;font-size:.8rem}.licence-file-guide details>p{max-width:88ch;margin:8px 0 13px;color:var(--text-secondary,#555);line-height:1.45}
      .licence-file-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}.licence-file-grid article{padding:12px;border:1px solid var(--border-subtle,#ddd);border-radius:9px;background:var(--bg-primary,#fff)}.licence-file-grid article>span{font-size:.72rem;text-transform:uppercase;letter-spacing:.06em;color:var(--text-secondary,#666)}.licence-file-grid strong,.licence-file-grid p,.licence-file-grid code{display:block}.licence-file-grid strong{margin-top:4px}.licence-file-grid p{margin:5px 0;font-size:.89rem;line-height:1.42;color:var(--text-secondary,#555)}.licence-file-grid code{margin-top:8px;font-size:.8rem;overflow-wrap:anywhere}
      .licence-id-chain{margin-top:10px;padding:12px;border:1px solid var(--border-subtle,#ddd);border-radius:9px;background:var(--bg-primary,#fff)}.licence-id-chain>strong{display:block;margin-bottom:9px}.licence-id-chain>div{display:flex;align-items:center;gap:8px;flex-wrap:wrap}.licence-id-chain span{padding:7px 9px;border:1px solid var(--border-subtle,#ddd);border-radius:8px;font-family:var(--font-mono,monospace);font-size:.82rem}.licence-id-chain small{font-family:inherit;color:var(--text-secondary,#666)}.licence-id-chain i{font-style:normal;opacity:.6}.licence-file-guide details{margin-top:11px}.licence-file-guide summary{cursor:pointer;font-weight:700}.licence-file-guide pre{overflow:auto;padding:10px;border:1px solid var(--border-subtle,#ddd);border-radius:8px;background:var(--bg-primary,#fff);font-size:.8rem}.licence-file-guide details p{font-size:.88rem;color:var(--text-secondary,#555)}
      @media(max-width:760px){.licence-file-head,.licence-file-grid{grid-template-columns:1fr}.licence-file-head>b{justify-self:start}.licence-id-chain>div{align-items:stretch;flex-direction:column}.licence-id-chain i{transform:rotate(90deg);align-self:flex-start;margin-left:20px}}
    `;
    document.head.appendChild(style);
  }

  let queued = false;
  function refresh() {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => {
      queued = false;
      installStyle();
      renderFileGuide();
    });
  }

  const observer = new MutationObserver(refresh);
  observer.observe(document.documentElement, { childList:true, subtree:true, attributes:true, attributeFilter:['lang'] });
  document.addEventListener('DOMContentLoaded', refresh, { once:true });
  refresh();
})();
