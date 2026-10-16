// @ts-nocheck
const LAYERS = ['transcriptomics', 'proteomics', 'metabolomics'];
const MULTIOMICS_ENGINE_VERSION = '1.3.0';
const CHEBI_SEARCH_URL = 'https://www.ebi.ac.uk/chebi/backend/api/public/es_search/';
const metaboliteResolutionCache = new Map();

const METABOLITE_ALIASES = new Map(Object.entries({
  'c14.0': 'myristic acid',
  'c16.0': 'palmitic acid',
  'c16.1n.7': 'palmitoleic acid',
  'c18.0': 'stearic acid',
  'c18.1n.9': 'oleic acid',
  'c18.1n.7': 'vaccenic acid',
  'c18.2n.6': 'linoleic acid',
  'c18.3n.6': 'gamma-linolenic acid',
  'c18.3n.3': 'alpha-linolenic acid',
  'c20.1n.9': 'gadoleic acid',
  'c20.2n.6': 'eicosadienoic acid',
  'c20.3n.6': 'dihomo-gamma-linolenic acid',
  'c20.4n.6': 'arachidonic acid',
  'c20.5n.3': 'eicosapentaenoic acid',
  'c22.4n.6': 'adrenic acid',
  'c22.6n.3': 'docosahexaenoic acid'
}));

function normaliseText(value) {
  return String(value ?? '').trim();
}

function lexicalKey(value) {
  return normaliseText(value)
    .toLowerCase()
    .replace(/[α]/g, 'alpha')
    .replace(/[β]/g, 'beta')
    .replace(/[γ]/g, 'gamma')
    .replace(/[^a-z0-9]+/g, ' ')
    .trim()
    .replace(/\s+/g, ' ');
}

function directIdentifierType(value) {
  const id = normaliseText(value);
  if (/^CHEBI:\d+$/i.test(id)) return 'chebi';
  if (/^HMDB\d+$/i.test(id)) return 'hmdb';
  if (/^C\d{5}$/i.test(id)) return 'kegg';
  if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/i.test(id)) return 'inchikey';
  return null;
}

function candidateChebiSource(result) {
  return result?._source || result?.source || result || null;
}

export async function resolveMetaboliteIdentifier(identifier, { fetchFn = fetch } = {}) {
  const original = normaliseText(identifier);
  if (!original) return { original, resolved: null, status: 'empty', method: 'none', query: '' };

  const direct = directIdentifierType(original);
  if (direct === 'chebi') {
    return { original, resolved: original.toUpperCase(), status: 'canonical', method: 'input', query: original };
  }

  if (direct === 'inchikey') {
    try {
      const response = await fetchFn(
        'https://www.ebi.ac.uk/unichem/rest/verbose_inchikey/' + encodeURIComponent(original.toUpperCase()),
        { headers: { Accept: 'application/json' } }
      );
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const json = await response.json();
      const rows = Array.isArray(json) ? json : [];
      const chebiIds = [...new Set(rows
        .filter((item) => String(item?.name || item?.name_label || '').toLowerCase() === 'chebi' || Number(item?.src_id) === 7)
        .flatMap((item) => Array.isArray(item?.src_compound_id) ? item.src_compound_id : [item?.src_compound_id])
        .map(normaliseText)
        .filter((id) => /^CHEBI:\d+$/i.test(id))
        .map((id) => id.toUpperCase()))];
      if (chebiIds.length === 1) {
        return {
          original,
          resolved: chebiIds[0],
          status: 'resolved_external_id',
          method: 'unichem_inchikey_to_chebi',
          query: original,
          crossReferences: rows.length
        };
      }
      return {
        original,
        resolved: null,
        status: chebiIds.length > 1 ? 'ambiguous' : 'unresolved',
        method: 'unichem_inchikey_to_chebi',
        query: original,
        candidates: chebiIds
      };
    } catch (error) {
      return {
        original,
        resolved: null,
        status: 'api_error',
        method: 'unichem_inchikey_to_chebi',
        query: original,
        error: error instanceof Error ? error.message : 'UniChem request failed'
      };
    }
  }

  const alias = METABOLITE_ALIASES.get(original.toLowerCase());
  const query = alias || original;
  const cacheKey = lexicalKey(query);
  if (metaboliteResolutionCache.has(cacheKey)) {
    const cached = metaboliteResolutionCache.get(cacheKey);
    return { ...cached, original, query };
  }

  try {
    const url = new URL(CHEBI_SEARCH_URL);
    url.searchParams.set('term', query);
    url.searchParams.set('page', '1');
    url.searchParams.set('size', '8');
    const response = await fetchFn(url.toString(), { headers: { Accept: 'application/json' } });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const json = await response.json();
    const results = Array.isArray(json?.results) ? json.results : [];
    const queryKey = lexicalKey(query);
    const sources = results.map(candidateChebiSource).filter(Boolean);
    const exact = direct === 'hmdb' || direct === 'kegg'
      ? sources[0]
      : sources.find((source) => {
          const names = [source.name, source.ascii_name].filter(Boolean).map(lexicalKey);
          return names.includes(queryKey);
        });

    if (!exact?.chebi_accession) {
      const unresolved = {
        resolved: null,
        status: 'unresolved',
        method: alias ? 'chebi_exact_after_alias' : 'chebi_exact',
        query,
        label: null
      };
      metaboliteResolutionCache.set(cacheKey, unresolved);
      return { original, ...unresolved };
    }

    const resolved = {
      resolved: String(exact.chebi_accession).toUpperCase(),
      status: direct === 'hmdb' || direct === 'kegg' ? 'resolved_external_id' : 'resolved',
      method: direct === 'hmdb' || direct === 'kegg'
        ? `chebi_xref_${direct}`
        : (alias ? 'chebi_exact_after_alias' : 'chebi_exact'),
      query,
      label: exact.ascii_name || exact.name || query,
      stars: exact.stars ?? null
    };
    metaboliteResolutionCache.set(cacheKey, resolved);
    return { original, ...resolved };
  } catch (error) {
    return {
      original,
      resolved: null,
      status: 'api_error',
      method: alias ? 'chebi_exact_after_alias' : 'chebi_exact',
      query,
      error: error instanceof Error ? error.message : 'ChEBI request failed'
    };
  }
}

export async function resolveMetaboliteIdentifiers(identifiers, { fetchFn = fetch, maxQueries = 30 } = {}) {
  const unique = [...new Set(identifiers.map(normaliseText).filter(Boolean))];
  const results = [];
  let networkQueries = 0;
  for (const identifier of unique) {
    const direct = directIdentifierType(identifier);
    const needsNetwork = direct !== 'chebi';
    const cached = metaboliteResolutionCache.has(lexicalKey(METABOLITE_ALIASES.get(identifier.toLowerCase()) || identifier));
    if (needsNetwork && !cached && networkQueries >= maxQueries) {
      results.push({ original: identifier, resolved: null, status: 'query_limit', method: 'none', query: identifier });
      continue;
    }
    if (needsNetwork && !cached) networkQueries += 1;
    results.push(await resolveMetaboliteIdentifier(identifier, { fetchFn }));
  }
  return {
    mappings: results,
    resolvedCount: results.filter((x) => x.resolved).length,
    unresolvedCount: results.filter((x) => !x.resolved).length,
    networkQueries
  };
}

const geneResolutionCache = new Map();
const proteinResolutionCache = new Map();

function organismLookup(organism) {
  const key = normaliseText(organism).toLowerCase();
  if (['human','homo_sapiens','homo sapiens'].includes(key)) return { ensembl: 'homo_sapiens', taxon: 9606 };
  if (['mouse','mus_musculus','mus musculus'].includes(key)) return { ensembl: 'mus_musculus', taxon: 10090 };
  if (['rat','rattus_norvegicus','rattus norvegicus'].includes(key)) return { ensembl: 'rattus_norvegicus', taxon: 10116 };
  return null;
}

export async function resolveGeneIdentifier(identifier, { identifierType = 'unknown', organism = 'human', fetchFn = fetch } = {}) {
  const original = normaliseText(identifier);
  if (!original) return { original, resolved: null, status: 'empty', method: 'none' };
  const stable = original.replace(/\.\d+$/, '');
  if (/^ENSG\d+$/i.test(stable) || /^ENSMUSG\d+$/i.test(stable) || /^ENSRNOG\d+$/i.test(stable)) {
    return { original, resolved: stable.toUpperCase(), status: 'canonical', method: 'input_ensembl' };
  }
  const org = organismLookup(organism);
  if (!org) return { original, resolved: null, status: 'unsupported_organism', method: 'none' };
  const cacheKey = org.ensembl + '|' + identifierType + '|' + original.toUpperCase();
  if (geneResolutionCache.has(cacheKey)) return { original, ...geneResolutionCache.get(cacheKey) };

  try {
    let url;
    if (identifierType === 'entrez' || /^\d+$/.test(original)) {
      url = 'https://rest.ensembl.org/xrefs/name/' + org.ensembl + '/' + encodeURIComponent(original) + '?content-type=application/json';
      const response = await fetchFn(url, { headers: { Accept: 'application/json' } });
      if (!response.ok) throw new Error('HTTP ' + response.status);
      const json = await response.json();
      const genes = (Array.isArray(json) ? json : []).filter((item) => String(item.type || '').toLowerCase() === 'gene' || /^ENS.*G\d+/i.test(String(item.id || '')));
      if (genes.length !== 1) {
        const result = { resolved: null, status: genes.length > 1 ? 'ambiguous' : 'unresolved', method: 'ensembl_xrefs_name', candidates: genes.slice(0,5).map((x) => x.id) };
        geneResolutionCache.set(cacheKey, result);
        return { original, ...result };
      }
      const result = { resolved: String(genes[0].id).replace(/\.\d+$/, ''), status: 'resolved', method: 'ensembl_xrefs_name', label: genes[0].display_id || original };
      geneResolutionCache.set(cacheKey, result);
      return { original, ...result };
    }

    url = 'https://rest.ensembl.org/lookup/symbol/' + org.ensembl + '/' + encodeURIComponent(original) + '?content-type=application/json';
    const response = await fetchFn(url, { headers: { Accept: 'application/json' } });
    if (response.status === 400 || response.status === 404) {
      const result = { resolved: null, status: 'unresolved', method: 'ensembl_lookup_symbol' };
      geneResolutionCache.set(cacheKey, result);
      return { original, ...result };
    }
    if (!response.ok) throw new Error('HTTP ' + response.status);
    const json = await response.json();
    const id = normaliseText(json?.id).replace(/\.\d+$/, '');
    const result = id
      ? { resolved: id, status: 'resolved', method: 'ensembl_lookup_symbol', label: json?.display_name || original, biotype: json?.biotype || null }
      : { resolved: null, status: 'unresolved', method: 'ensembl_lookup_symbol' };
    geneResolutionCache.set(cacheKey, result);
    return { original, ...result };
  } catch (error) {
    return { original, resolved: null, status: 'api_error', method: 'ensembl', error: error instanceof Error ? error.message : 'Ensembl request failed' };
  }
}

export async function resolveGeneIdentifiers(identifiers, options = {}) {
  const unique = [...new Set(identifiers.map(normaliseText).filter(Boolean))];
  const maxQueries = options.maxQueries ?? 30;
  const mappings = [];
  let networkQueries = 0;
  for (const identifier of unique) {
    const direct = /^ENS(?:G|MUSG|RNOG)\d+/i.test(identifier);
    if (!direct && networkQueries >= maxQueries) {
      mappings.push({ original: identifier, resolved: null, status: 'query_limit', method: 'none' });
      continue;
    }
    if (!direct) networkQueries += 1;
    mappings.push(await resolveGeneIdentifier(identifier, options));
  }
  return {
    mappings,
    resolvedCount: mappings.filter((item) => item.resolved).length,
    unresolvedCount: mappings.filter((item) => !item.resolved).length,
    networkQueries
  };
}

export async function resolveProteinIdentifier(identifier, { identifierType = 'unknown', organism = 'human', fetchFn = fetch } = {}) {
  const original = normaliseText(identifier);
  if (!original) return { original, resolved: null, status: 'empty', method: 'none' };
  const accession = original.replace(/-\d+$/, '');
  if (/^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])$/i.test(accession)) {
    return { original, resolved: accession.toUpperCase(), status: 'canonical', method: 'input_uniprot' };
  }
  const org = organismLookup(organism);
  if (!org) return { original, resolved: null, status: 'unsupported_organism', method: 'none' };
  const cacheKey = org.taxon + '|' + identifierType + '|' + original.toUpperCase();
  if (proteinResolutionCache.has(cacheKey)) return { original, ...proteinResolutionCache.get(cacheKey) };

  try {
    if (/^ENS.*P\d+/i.test(original) || identifierType === 'ensembl_protein') {
      const url = 'https://rest.ensembl.org/xrefs/id/' + encodeURIComponent(original.replace(/\.\d+$/, '')) + '?content-type=application/json';
      const response = await fetchFn(url, { headers: { Accept: 'application/json' } });
      if (!response.ok) throw new Error('HTTP ' + response.status);
      const json = await response.json();
      const candidates = (Array.isArray(json) ? json : [])
        .filter((item) => /uniprot|swiss/i.test(String(item.dbname || item.db_display_name || '')))
        .map((item) => normaliseText(item.primary_id || item.display_id))
        .filter((id) => /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])$/i.test(id));
      const unique = [...new Set(candidates)];
      const result = unique.length === 1
        ? { resolved: unique[0], status: 'resolved', method: 'ensembl_xrefs_uniprot' }
        : { resolved: null, status: unique.length > 1 ? 'ambiguous' : 'unresolved', method: 'ensembl_xrefs_uniprot', candidates: unique.slice(0,5) };
      proteinResolutionCache.set(cacheKey, result);
      return { original, ...result };
    }

    const query = '(gene_exact:' + original.replace(/[^A-Za-z0-9_.-]/g, '') + ')+AND+(organism_id:' + org.taxon + ')';
    const url = 'https://rest.uniprot.org/uniprotkb/search?query=' + encodeURIComponent(query) + '&fields=accession,gene_names,reviewed&format=json&size=5';
    const response = await fetchFn(url, { headers: { Accept: 'application/json' } });
    if (!response.ok) throw new Error('HTTP ' + response.status);
    const json = await response.json();
    const results = Array.isArray(json?.results) ? json.results : [];
    const reviewed = results.filter((item) => String(item.entryType || '').toLowerCase().includes('reviewed'));
    const pool = reviewed.length ? reviewed : results;
    const accessions = [...new Set(pool.map((item) => normaliseText(item.primaryAccession)).filter(Boolean))];
    const result = accessions.length === 1
      ? { resolved: accessions[0], status: 'resolved', method: 'uniprot_gene_search' }
      : { resolved: null, status: accessions.length > 1 ? 'ambiguous' : 'unresolved', method: 'uniprot_gene_search', candidates: accessions.slice(0,5) };
    proteinResolutionCache.set(cacheKey, result);
    return { original, ...result };
  } catch (error) {
    return { original, resolved: null, status: 'api_error', method: 'uniprot_or_ensembl', error: error instanceof Error ? error.message : 'Protein identifier request failed' };
  }
}

export async function resolveProteinIdentifiers(identifiers, options = {}) {
  const unique = [...new Set(identifiers.map(normaliseText).filter(Boolean))];
  const maxQueries = options.maxQueries ?? 30;
  const mappings = [];
  let networkQueries = 0;
  for (const identifier of unique) {
    const direct = /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])(?:-\d+)?$/i.test(identifier);
    if (!direct && networkQueries >= maxQueries) {
      mappings.push({ original: identifier, resolved: null, status: 'query_limit', method: 'none' });
      continue;
    }
    if (!direct) networkQueries += 1;
    mappings.push(await resolveProteinIdentifier(identifier, options));
  }
  return {
    mappings,
    resolvedCount: mappings.filter((item) => item.resolved).length,
    unresolvedCount: mappings.filter((item) => !item.resolved).length,
    networkQueries
  };
}

function canonicalOmic(value) {
  const key = normaliseText(value).toLowerCase().replace(/[^a-z0-9]+/g, '_');
  if (['transcriptomics','transcriptome','transcriptomique','rna','rnaseq','rna_seq','mrna','arn','gene_expression'].includes(key)) return 'transcriptomics';
  if (['proteomics','proteome','proteomique','protein','proteins','proteine','proteines','lfq'].includes(key)) return 'proteomics';
  if (['metabolomics','metabolome','metabolomique','metabolite','metabolites','met'].includes(key)) return 'metabolomics';
  return key;
}

function splitDelimitedLine(line, delimiter) {
  const out = [];
  let current = '';
  let quoted = false;
  for (let i = 0; i < line.length; i += 1) {
    const char = line[i];
    if (char === '"') {
      if (quoted && line[i + 1] === '"') {
        current += '"';
        i += 1;
      } else quoted = !quoted;
    } else if (char === delimiter && !quoted) {
      out.push(current.trim());
      current = '';
    } else current += char;
  }
  out.push(current.trim());
  return out;
}

function detectDelimiter(text) {
  const first = text.split(/\r?\n/).find((line) => line.trim()) ?? '';
  return [',','\t',';']
    .map((delimiter) => ({ delimiter, n: splitDelimitedLine(first, delimiter).length }))
    .sort((a,b) => b.n - a.n)[0]?.delimiter ?? ',';
}

export function parseDelimited(text) {
  const delimiter = detectDelimiter(text);
  const lines = text.split(/\r?\n/).filter((line) => line.trim());
  if (!lines.length) return { headers: [], rows: [], delimiter };
  const headers = splitDelimitedLine(lines[0], delimiter);
  const rows = lines.slice(1).map((line) => {
    const values = splitDelimitedLine(line, delimiter);
    return Object.fromEntries(headers.map((header, i) => [header, values[i] ?? '']));
  });
  return { headers, rows, delimiter };
}

function finiteNumber(value) {
  if (value === '' || value == null) return null;
  const x = Number(String(value).replace(',', '.'));
  return Number.isFinite(x) ? x : null;
}

function mean(values) {
  const x = values.filter(Number.isFinite);
  return x.length ? x.reduce((a,b) => a+b, 0) / x.length : NaN;
}

function variance(values) {
  const x = values.filter(Number.isFinite);
  if (x.length < 2) return NaN;
  const m = mean(x);
  return x.reduce((sum, value) => sum + (value-m)*(value-m), 0) / (x.length-1);
}

function numericTime(value) {
  const match = String(value ?? '').match(/-?\d+(?:\.\d+)?/);
  return match ? Number(match[0]) : NaN;
}

function linearSlope(points) {
  const valid = points.filter((point) => Number.isFinite(point.time) && Number.isFinite(point.value));
  if (valid.length < 2) return NaN;
  const mt = mean(valid.map((point) => point.time));
  const my = mean(valid.map((point) => point.value));
  let numerator = 0;
  let denominator = 0;
  for (const point of valid) {
    numerator += (point.time-mt)*(point.value-my);
    denominator += (point.time-mt)*(point.time-mt);
  }
  return denominator > 0 ? numerator/denominator : NaN;
}

function median(values) {
  const x = values.filter(Number.isFinite).slice().sort((a,b) => a-b);
  if (!x.length) return NaN;
  const m = Math.floor(x.length / 2);
  return x.length % 2 ? x[m] : (x[m-1] + x[m]) / 2;
}

function ranks(values) {
  const indexed = values.map((value, index) => ({ value, index })).sort((a,b) => a.value-b.value);
  const out = new Array(values.length);
  let i = 0;
  while (i < indexed.length) {
    let j = i + 1;
    while (j < indexed.length && indexed[j].value === indexed[i].value) j += 1;
    const rank = (i + 1 + j) / 2;
    for (let k = i; k < j; k += 1) out[indexed[k].index] = rank;
    i = j;
  }
  return out;
}

function pearson(x, y) {
  if (x.length !== y.length || x.length < 3) return NaN;
  const mx = mean(x);
  const my = mean(y);
  let num = 0;
  let dx = 0;
  let dy = 0;
  for (let i = 0; i < x.length; i += 1) {
    const a = x[i] - mx;
    const b = y[i] - my;
    num += a*b;
    dx += a*a;
    dy += b*b;
  }
  return dx > 0 && dy > 0 ? num / Math.sqrt(dx*dy) : NaN;
}

function spearman(x, y) {
  if (x.length !== y.length || x.length < 3) return NaN;
  return pearson(ranks(x), ranks(y));
}

function logGamma(z) {
  const coefficients = [
    676.5203681218851,-1259.1392167224028,771.32342877765313,
    -176.61502916214059,12.507343278686905,-0.13857109526572012,
    9.9843695780195716e-6,1.5056327351493116e-7
  ];
  if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI*z)) - logGamma(1-z);
  z -= 1;
  let x = 0.99999999999980993;
  for (let i=0;i<coefficients.length;i+=1) x += coefficients[i]/(z+i+1);
  const t = z + coefficients.length - 0.5;
  return 0.5*Math.log(2*Math.PI) + (z+0.5)*Math.log(t) - t + Math.log(x);
}

function betaContinuedFraction(a,b,x) {
  const maxIter = 200;
  const eps = 3e-12;
  const fpmin = 1e-300;
  const qab = a+b;
  const qap = a+1;
  const qam = a-1;
  let c = 1;
  let d = 1 - qab*x/qap;
  if (Math.abs(d) < fpmin) d = fpmin;
  d = 1/d;
  let h = d;
  for (let m=1;m<=maxIter;m+=1) {
    const m2=2*m;
    let aa = m*(b-m)*x/((qam+m2)*(a+m2));
    d = 1 + aa*d; if (Math.abs(d)<fpmin) d=fpmin;
    c = 1 + aa/c; if (Math.abs(c)<fpmin) c=fpmin;
    d=1/d; h*=d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1 + aa*d; if (Math.abs(d)<fpmin) d=fpmin;
    c = 1 + aa/c; if (Math.abs(c)<fpmin) c=fpmin;
    d=1/d;
    const del=d*c;
    h*=del;
    if (Math.abs(del-1)<eps) break;
  }
  return h;
}

function regularizedBeta(x,a,b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(logGamma(a+b)-logGamma(a)-logGamma(b)+a*Math.log(x)+b*Math.log(1-x));
  return x < (a+1)/(a+b+2)
    ? bt*betaContinuedFraction(a,b,x)/a
    : 1-bt*betaContinuedFraction(b,a,1-x)/b;
}

function studentTCdf(t,df) {
  if (!Number.isFinite(t) || !Number.isFinite(df) || df <= 0) {
    if (t === Infinity) return 1;
    if (t === -Infinity) return 0;
    return NaN;
  }
  const x=df/(df+t*t);
  const ib=regularizedBeta(x,df/2,0.5);
  return t >= 0 ? 1-0.5*ib : 0.5*ib;
}

// Inference using a Student-t p-value must use the same distribution for its
// confidence interval, particularly with small cohorts.
const t95Cache = new Map();
function tCritical95(df) {
  if (!(Number.isFinite(df) && df > 0)) return NaN;
  if (t95Cache.has(df)) return t95Cache.get(df);
  let lo = 0, hi = 20;
  for (let iter = 0; iter < 42; iter += 1) {
    const mid = (lo + hi) / 2;
    if (studentTCdf(mid, df) >= 0.975) hi = mid;
    else lo = mid;
  }
  const critical = (lo + hi) / 2;
  t95Cache.set(df, critical);
  return critical;
}

function fCdf(f,d1,d2) {
  if (f <= 0) return 0;
  if (f === Infinity) return 1;
  if (![f,d1,d2].every(Number.isFinite) || d1<=0 || d2<=0) return NaN;
  return regularizedBeta((d1*f)/(d1*f+d2), d1/2, d2/2);
}

function welchPValue(groupA,groupB) {
  const a=groupA.filter(Number.isFinite);
  const b=groupB.filter(Number.isFinite);
  if (a.length<2 || b.length<2) return null;
  const va=variance(a), vb=variance(b);
  const se2=va/a.length + vb/b.length;
  if (!(se2>0)) return mean(a)===mean(b) ? 1 : 0;
  const t=(mean(b)-mean(a))/Math.sqrt(se2);
  const aTerm=(va/a.length)**2/(a.length-1);
  const bTerm=(vb/b.length)**2/(b.length-1);
  const df=se2*se2/(aTerm+bTerm);
  const cdf=studentTCdf(Math.abs(t),df);
  return Number.isFinite(cdf) ? Math.max(0,Math.min(1,2*(1-cdf))) : null;
}

function pairedTPValue(differences) {
  const d=differences.filter(Number.isFinite);
  if (d.length<2) return null;
  const v=variance(d);
  if (!(v>0)) return mean(d)===0 ? 1 : 0;
  const t=mean(d)/(Math.sqrt(v/d.length));
  const cdf=studentTCdf(Math.abs(t),d.length-1);
  return Number.isFinite(cdf) ? Math.max(0,Math.min(1,2*(1-cdf))) : null;
}

function oneWayAnovaPValue(groups) {
  const clean=groups.map(g=>g.filter(Number.isFinite)).filter(g=>g.length);
  const k=clean.length;
  const n=clean.reduce((s,g)=>s+g.length,0);
  const f=oneWayF(clean);
  if (k<2 || n<=k || (!Number.isFinite(f) && f!==Infinity)) return null;
  const cdf=fCdf(f,k-1,n-k);
  return Number.isFinite(cdf) ? Math.max(0,Math.min(1,1-cdf)) : null;
}

function normalCdf(x) {
  const t = 1 / (1 + 0.2316419 * Math.abs(x));
  const d = 0.3989422804014327 * Math.exp(-x*x/2);
  const p = 1 - d*t*(0.319381530 + t*(-0.356563782 + t*(1.781477937 + t*(-1.821255978 + t*1.330274429))));
  return x >= 0 ? p : 1-p;
}

export function differentialCorrelationPair(referenceX, referenceY, comparisonX, comparisonY, method = 'spearman') {
  const corr = method === 'pearson' ? pearson : spearman;
  const rReference = corr(referenceX, referenceY);
  const rComparison = corr(comparisonX, comparisonY);
  if (!Number.isFinite(rReference) || !Number.isFinite(rComparison)) {
    return { rReference, rComparison, deltaR: NaN, z: NaN, pValue: null, method };
  }
  const clamp = (r) => Math.max(-0.999999, Math.min(0.999999, r));
  const zReference = Math.atanh(clamp(rReference));
  const zComparison = Math.atanh(clamp(rComparison));
  const se = referenceX.length > 3 && comparisonX.length > 3
    ? Math.sqrt(1/(referenceX.length-3) + 1/(comparisonX.length-3))
    : NaN;
  const z = Number.isFinite(se) && se > 0 ? (zComparison-zReference)/se : NaN;
  const pValue = Number.isFinite(z) ? Math.min(1, 2*(1-normalCdf(Math.abs(z)))) : null;
  return {
    rReference,
    rComparison,
    deltaR: rComparison-rReference,
    z,
    pValue,
    method
  };
}

function hashString(value) {
  let h = 2166136261;
  for (let i = 0; i < value.length; i += 1) {
    h ^= value.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return h >>> 0;
}

function rngFromSeed(seed) {
  let state = seed >>> 0;
  return () => {
    state += 0x6D2B79F5;
    let t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function shuffled(array, rng) {
  const out = array.slice();
  for (let i = out.length - 1; i > 0; i -= 1) {
    const j = Math.floor(rng() * (i + 1));
    [out[i], out[j]] = [out[j], out[i]];
  }
  return out;
}

function chooseCount(n, k) {
  k = Math.min(k, n-k);
  let value = 1;
  for (let i = 1; i <= k; i += 1) value = value * (n-k+i) / i;
  return Math.round(value);
}

function combinations(n, k, limit = 5000) {
  const total = chooseCount(n,k);
  if (total > limit) return null;
  const out = [];
  const current = [];
  function walk(start, left) {
    if (left === 0) { out.push(current.slice()); return; }
    for (let i = start; i <= n-left; i += 1) {
      current.push(i);
      walk(i+1, left-1);
      current.pop();
    }
  }
  walk(0,k);
  return out;
}

function permutationPValue(groupA, groupB, seedKey) {
  const a = groupA.filter(Number.isFinite);
  const b = groupB.filter(Number.isFinite);
  if (a.length < 2 || b.length < 2) return null;
  const observed = Math.abs(mean(b) - mean(a));
  const all = [...a, ...b];
  const nA = a.length;
  const exact = combinations(all.length, nA, 4096);
  let extreme = 0;
  let total = 0;
  if (exact) {
    for (const indices of exact) {
      const set = new Set(indices);
      const pa = all.filter((_,i) => set.has(i));
      const pb = all.filter((_,i) => !set.has(i));
      if (Math.abs(mean(pb) - mean(pa)) >= observed - 1e-12) extreme += 1;
      total += 1;
    }
    return total ? extreme / total : null;
  }
  const rng = rngFromSeed(hashString(seedKey));
  total = 4096;
  for (let iter = 0; iter < total; iter += 1) {
    const perm = shuffled(all, rng);
    const pa = perm.slice(0,nA);
    const pb = perm.slice(nA);
    if (Math.abs(mean(pb) - mean(pa)) >= observed - 1e-12) extreme += 1;
  }
  return (extreme + 1) / (total + 1);
}

function oneWayF(groups) {
  const clean = groups.map((group) => group.filter(Number.isFinite)).filter((group) => group.length);
  const k = clean.length;
  const n = clean.reduce((sum, group) => sum + group.length, 0);
  if (k < 2 || n <= k) return NaN;
  const grand = mean(clean.flat());
  let ssBetween = 0;
  let ssWithin = 0;
  for (const group of clean) {
    const m = mean(group);
    ssBetween += group.length * (m-grand) * (m-grand);
    for (const value of group) ssWithin += (value-m) * (value-m);
  }
  if (ssWithin <= 0) return ssBetween > 0 ? Infinity : 0;
  return (ssBetween/(k-1)) / (ssWithin/(n-k));
}

function multiGroupPermutationPValue(groups, seedKey) {
  const clean = groups.map((group) => group.filter(Number.isFinite));
  if (clean.length < 3 || clean.some((group) => group.length < 2)) return null;
  const observed = oneWayF(clean);
  if (!Number.isFinite(observed) && observed !== Infinity) return null;
  const sizes = clean.map((group) => group.length);
  const values = clean.flat();
  const rng = rngFromSeed(hashString(seedKey));
  const total = 4096;
  let extreme = 0;
  for (let iter = 0; iter < total; iter += 1) {
    const perm = shuffled(values, rng);
    const permGroups = [];
    let cursor = 0;
    for (const size of sizes) {
      permGroups.push(perm.slice(cursor, cursor + size));
      cursor += size;
    }
    const stat = oneWayF(permGroups);
    if (stat >= observed - 1e-12) extreme += 1;
  }
  return (extreme + 1)/(total + 1);
}

function independentMultiGroupValues(aggregated, feature) {
  const source = aggregated.values.get(feature);
  const timepoints = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.timepoint));
  const targetTime = timepoints.length ? timepoints[timepoints.length-1] : '';
  const conditions = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.condition));
  const perConditionSubjects = new Map(conditions.map((condition) => [condition, new Map()]));

  for (const [sampleId, row] of aggregated.sampleMeta.entries()) {
    const value = source.get(sampleId);
    if (!Number.isFinite(value)) continue;
    if (targetTime && row.timepoint !== targetTime) continue;
    if (!perConditionSubjects.has(row.condition)) continue;
    const subjects = perConditionSubjects.get(row.condition);
    if (!subjects.has(row.subjectId)) subjects.set(row.subjectId, []);
    subjects.get(row.subjectId).push(value);
  }

  const groups = conditions.map((condition) =>
    [...perConditionSubjects.get(condition).values()].map((values) => mean(values))
  );
  return { conditions, groups, timepoints, targetTime };
}

function pairedPermutationPValue(differences, seedKey) {
  const d = differences.filter(Number.isFinite);
  if (d.length < 2) return null;
  const observed = Math.abs(mean(d));
  if (d.length <= 12) {
    const total = 2 ** d.length;
    let extreme = 0;
    for (let mask = 0; mask < total; mask += 1) {
      const signed = d.map((value, i) => (mask & (1 << i)) ? value : -value);
      if (Math.abs(mean(signed)) >= observed - 1e-12) extreme += 1;
    }
    return extreme / total;
  }
  const rng = rngFromSeed(hashString(seedKey));
  const total = 4096;
  let extreme = 0;
  for (let iter = 0; iter < total; iter += 1) {
    const signed = d.map((value) => rng() < 0.5 ? value : -value);
    if (Math.abs(mean(signed)) >= observed - 1e-12) extreme += 1;
  }
  return (extreme + 1)/(total + 1);
}

function bhAdjust(rows) {
  const valid = rows.map((row,index) => ({ index, p: row.pValue })).filter((x) => Number.isFinite(x.p)).sort((a,b) => a.p-b.p);
  const m = valid.length;
  let previous = 1;
  for (let i = m-1; i >= 0; i -= 1) {
    const q = Math.min(previous, valid[i].p * m / (i+1), 1);
    rows[valid[i].index].qValue = q;
    previous = q;
  }
  return rows;
}

function naturalOrder(values) {
  return [...new Set(values.filter(Boolean))].sort((a,b) => {
    const na = Number(String(a).match(/-?\d+(?:\.\d+)?/)?.[0]);
    const nb = Number(String(b).match(/-?\d+(?:\.\d+)?/)?.[0]);
    if (Number.isFinite(na) && Number.isFinite(nb) && na !== nb) return na-nb;
    return String(a).localeCompare(String(b), undefined, {numeric:true});
  });
}

function canonicalMetadata(metadataRows, columnMapping, covariateColumns = []) {
  return metadataRows.map((row) => {
    const get = (key) => {
      const column = columnMapping[key];
      return column ? normaliseText(row[column]) : '';
    };
    const covariates = Object.fromEntries(
      covariateColumns
        .filter((column) => column && Object.prototype.hasOwnProperty.call(row, column))
        .map((column) => [column, normaliseText(row[column])])
    );
    return {
      subjectId: get('subject_id'),
      sampleId: get('sample_id'),
      assayId: get('assay_id'),
      omic: canonicalOmic(get('omic')),
      condition: get('condition'),
      timepoint: get('timepoint'),
      batch: get('batch'),
      technicalReplicate: get('technical_replicate'),
      outcome: get('outcome'),
      survivalTime: get('survival_time'),
      survivalEvent: get('survival_event'),
      sampleType: get('sample_type'),
      injectionOrder: get('injection_order'),
      covariates
    };
  }).filter((row) => row.subjectId && row.sampleId && row.assayId && LAYERS.includes(row.omic));
}

function assertSampleSheetIdentifiers(metadata) {
  const assayKeys = new Set();
  const samples = new Map();
  for (const row of metadata) {
    // Assay IDs are unique within each omic. Different omics may legitimately
    // use the same machine identifier, but never for two rows of one matrix.
    const assayKey = JSON.stringify([row.omic, row.assayId]);
    if (assayKeys.has(assayKey)) {
      throw new Error('Duplicate assay_id ' + row.assayId + ' within ' + row.omic
        + '. Every measurement column must map to exactly one sample-sheet row.');
    }
    assayKeys.add(assayKey);
    const old = samples.get(row.sampleId);
    if (old && (old.subjectId !== row.subjectId || old.timepoint !== row.timepoint)) {
      throw new Error('Conflicting sample_id ' + row.sampleId
        + ': the same biological sample cannot belong to different subject_id or timepoint values.');
    }
    samples.set(row.sampleId, { subjectId: row.subjectId, timepoint: row.timepoint });
  }
}

function matrixFromText(text, expectedAssays) {
  const parsed = parseDelimited(text);
  if (parsed.headers.length < 2) throw new Error('Matrix requires a feature column and at least one assay.');
  const first = parsed.headers[0];
  const headerAssays = parsed.headers.slice(1);
  if (headerAssays.some((id) => !normaliseText(id)) ||
      new Set(headerAssays.map(normaliseText)).size !== headerAssays.length) {
    throw new Error('Matrix contains missing or duplicate column identifiers. Sample or feature IDs must be unique.');
  }
  const rowIds = parsed.rows.map((row) => normaliseText(row[first])).filter(Boolean);
  const expected = new Set(expectedAssays);
  const columnMatches = headerAssays.filter((id) => expected.has(id)).length;
  const rowMatches = rowIds.filter((id) => expected.has(id)).length;

  if (rowMatches > columnMatches) {
    const featureHeaders = parsed.headers.slice(1);
    const assayRows = parsed.rows.filter((row) => expected.has(normaliseText(row[first])));
    const assayRowIds = assayRows.map((row) => normaliseText(row[first]));
    if (new Set(assayRowIds).size !== assayRowIds.length) {
      throw new Error('Transposed matrix contains duplicate assay_id rows. Each assay must appear exactly once.');
    }
    const values = new Map();
    for (const feature of featureHeaders) {
      const featureValues = new Map();
      for (const row of assayRows) {
        const assay = normaliseText(row[first]);
        featureValues.set(assay, finiteNumber(row[feature]));
      }
      values.set(feature, featureValues);
    }
    return { features: featureHeaders, assays: assayRows.map((row) => normaliseText(row[first])), values, transposed: true };
  }

  const values = new Map();
  for (const row of parsed.rows) {
    const feature = normaliseText(row[first]);
    if (!feature) continue;
    if (values.has(feature)) {
      throw new Error('Matrix contains duplicate feature_id ' + feature
        + '. Duplicate biological identifiers must be resolved explicitly before import.');
    }
    const featureValues = new Map();
    for (const assay of headerAssays) featureValues.set(assay, finiteNumber(row[assay]));
    values.set(feature, featureValues);
  }
  return { features: [...values.keys()], assays: headerAssays, values, transposed: false };
}

function mad(values) {
  const x = values.filter(Number.isFinite);
  if (!x.length) return NaN;
  const m = median(x);
  return median(x.map((value) => Math.abs(value - m)));
}

function robustOutlierFlags(values, threshold = 4) {
  const finite = values.filter(Number.isFinite);
  const center = median(finite);
  const spread = mad(finite);
  if (!Number.isFinite(center) || !Number.isFinite(spread) || spread === 0) return values.map(() => false);
  const scale = 1.4826 * spread;
  return values.map((value) => Number.isFinite(value) && Math.abs(value - center) / scale > threshold);
}

function assayGroupsForQc(metadataRows, assays) {
  const assaySet = new Set(assays);
  const rows = metadataRows.filter((row) => assaySet.has(row.assayId));
  const byCondition = new Map();
  for (const row of rows) {
    const key = row.condition || '__all__';
    if (!byCondition.has(key)) byCondition.set(key, []);
    byCondition.get(key).push(row.assayId);
  }
  if (!byCondition.size) byCondition.set('__all__', assays.slice());
  return byCondition;
}

function canonicalSampleType(value) {
  const key = lexicalKey(value).replace(/ /g, '_');
  if (!key) return 'biological';
  if (['blank','solvent_blank','process_blank','extraction_blank','method_blank'].includes(key)) return 'blank';
  if (['qc','pooled_qc','quality_control','quality_control_pool','pooled','pool_qc'].includes(key)) return 'qc';
  return 'biological';
}

function simpleQuantile(values, probability) {
  const x = values.filter(Number.isFinite).slice().sort((a,b) => a-b);
  if (!x.length) return NaN;
  if (x.length === 1) return x[0];
  const pos = Math.max(0, Math.min(1, probability)) * (x.length - 1);
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  if (lo === hi) return x[lo];
  return x[lo] + (pos - lo) * (x[hi] - x[lo]);
}

function localLinearQcTrend(xs, ys, target, span = 0.6) {
  const pairs = xs.map((x,i) => ({x,y:ys[i],d:Math.abs(x-target)}))
    .filter((item) => Number.isFinite(item.x) && Number.isFinite(item.y))
    .sort((a,b) => a.d-b.d);
  if (pairs.length < 3) return NaN;
  const k = Math.min(pairs.length, Math.max(4, Math.ceil(pairs.length * span)));
  const near = pairs.slice(0,k);
  const maxDistance = Math.max(...near.map((item) => item.d), 1e-12);
  const weights = near.map((item) => {
    const u = Math.min(1, item.d / maxDistance);
    return Math.pow(1 - Math.pow(u,3),3);
  });
  const sw = weights.reduce((a,b)=>a+b,0);
  if (!(sw > 0)) return mean(near.map((item)=>item.y));
  const xbar = near.reduce((sum,item,i)=>sum+weights[i]*item.x,0)/sw;
  const ybar = near.reduce((sum,item,i)=>sum+weights[i]*item.y,0)/sw;
  let num = 0;
  let den = 0;
  for (let i=0;i<near.length;i+=1) {
    num += weights[i]*(near[i].x-xbar)*(near[i].y-ybar);
    den += weights[i]*(near[i].x-xbar)*(near[i].x-xbar);
  }
  const slope = den > 1e-12 ? num/den : 0;
  return ybar + slope*(target-xbar);
}

function applyMetabolomicsMsQc(matrix, metadataRows, config = {}) {
  const metaByAssay = new Map(metadataRows.map((row) => [
    row.assayId,
    {
      ...row,
      canonicalSampleType: canonicalSampleType(row.sampleType),
      numericInjectionOrder: finiteNumber(row.injectionOrder)
    }
  ]));
  const blankAssays = matrix.assays.filter((assay) => metaByAssay.get(assay)?.canonicalSampleType === 'blank');
  const qcAssays = matrix.assays.filter((assay) => metaByAssay.get(assay)?.canonicalSampleType === 'qc');
  const biologicalAssays = matrix.assays.filter((assay) => {
    const type = metaByAssay.get(assay)?.canonicalSampleType || 'biological';
    return type === 'biological';
  });
  if (!biologicalAssays.length) {
    throw new Error('Metabolomics QC metadata contain no biological assays after excluding blank/QC injections.');
  }

  const blankFold = Number.isFinite(Number(config.blankFold)) ? Math.max(1, Number(config.blankFold)) : 5;
  const rsdThreshold = Number.isFinite(Number(config.qcRsdThreshold)) ? Math.max(0.01, Number(config.qcRsdThreshold)) : 0.30;
  const requestedBlankMode = String(config.blankMode ?? '').toLowerCase();
  const blankMode = ['flag','remove','off','no','yes'].includes(requestedBlankMode)
    ? (requestedBlankMode === 'no' || requestedBlankMode === 'off'
        ? 'off'
        : requestedBlankMode === 'yes' ? 'remove' : requestedBlankMode)
    : (config.blankFilter === false ? 'off' : config.blankFilter === true ? 'remove' : 'flag');
  const blankFilterEnabled = blankMode !== 'off';
  const driftEnabled = config.driftCorrection !== false;
  const qcRsdEnabled = config.qcRsdFilter !== false;
  const mnarStrategy = config.mnarStrategy || 'none';

  let features = matrix.features.slice();
  const values = new Map(matrix.features.map((feature) => [
    feature,
    new Map(matrix.assays.map((assay) => [assay, matrix.values.get(feature)?.get(assay)]))
  ]));

  const blankStats = [];
  if (blankFilterEnabled && blankAssays.length >= 2) {
    const keep = [];
    for (const feature of features) {
      const map = values.get(feature);
      const blanks = blankAssays.map((assay)=>map.get(assay)).filter(Number.isFinite);
      const biological = biologicalAssays.map((assay)=>map.get(assay)).filter(Number.isFinite);
      const blankMedian = median(blanks);
      const bioMedian = median(biological);
      const contaminated = Number.isFinite(blankMedian) && blankMedian > 0 &&
        Number.isFinite(bioMedian) && bioMedian < blankFold * blankMedian;
      blankStats.push({ feature, blankMedian, biologicalMedian: bioMedian, contaminated });
      if (!contaminated || blankMode === 'flag') keep.push(feature);
    }
    features = keep;
  }

  let driftCorrected = 0;
  let medianAbsoluteDriftLog = null;
  const driftMagnitudes = [];
  const qcOrders = qcAssays.map((assay)=>metaByAssay.get(assay)?.numericInjectionOrder);
  const orderedQcCount = qcOrders.filter(Number.isFinite).length;
  const allOrdered = matrix.assays.filter((assay)=>Number.isFinite(metaByAssay.get(assay)?.numericInjectionOrder));
  const qcOrderSupport = qcOrders.filter(Number.isFinite);
  const qcMinOrder = qcOrderSupport.length ? Math.min(...qcOrderSupport) : Infinity;
  const qcMaxOrder = qcOrderSupport.length ? Math.max(...qcOrderSupport) : -Infinity;
  const outsideQcSupport = biologicalAssays.filter((assay) => {
    const order = metaByAssay.get(assay)?.numericInjectionOrder;
    return Number.isFinite(order) && (order < qcMinOrder || order > qcMaxOrder);
  }).length;

  if (driftEnabled && qcAssays.length >= 5 && orderedQcCount >= 5 && allOrdered.length >= 5) {
    for (const feature of features) {
      const map = values.get(feature);
      const positives = matrix.assays.map((assay)=>map.get(assay)).filter((value)=>Number.isFinite(value) && value > 0);
      if (positives.length < 5) continue;
      const pseudocount = Math.max(Math.min(...positives) * 0.5, 1e-12);
      const qx = [];
      const qy = [];
      for (const assay of qcAssays) {
        const order = metaByAssay.get(assay)?.numericInjectionOrder;
        const value = map.get(assay);
        if (Number.isFinite(order) && Number.isFinite(value) && value >= 0) {
          qx.push(order);
          qy.push(Math.log(value + pseudocount));
        }
      }
      if (qx.length < 5) continue;
      const qcPredictions = qx.map((order)=>localLinearQcTrend(qx,qy,order));
      const reference = median(qcPredictions);
      if (!Number.isFinite(reference)) continue;

      let featureCorrected = false;
      for (const assay of allOrdered) {
        const raw = map.get(assay);
        const order = metaByAssay.get(assay)?.numericInjectionOrder;
        if (!Number.isFinite(raw) || raw < 0 || !Number.isFinite(order)) continue;
        // A correction must interpolate between observed QC positions; do
        // not extrapolate outside the supporting pooled-QC injection sequence.
        if (order < Math.min(...qx) || order > Math.max(...qx)) continue;
        const pred = localLinearQcTrend(qx,qy,order);
        if (!Number.isFinite(pred)) continue;
        const corrected = Math.max(0, Math.exp(Math.log(raw+pseudocount) - pred + reference) - pseudocount);
        map.set(assay, corrected);
        driftMagnitudes.push(Math.abs(pred-reference));
        featureCorrected = true;
      }
      if (featureCorrected) driftCorrected += 1;
    }
    medianAbsoluteDriftLog = median(driftMagnitudes);
  }

  const qcRsdStats = [];
  if (qcRsdEnabled && qcAssays.length >= 3) {
    const keep = [];
    for (const feature of features) {
      const map = values.get(feature);
      const qcValues = qcAssays.map((assay)=>map.get(assay)).filter((value)=>Number.isFinite(value) && value >= 0);
      let rsd = NaN;
      if (qcValues.length >= 3) {
        const m = mean(qcValues);
        rsd = m > 0 ? Math.sqrt(Math.max(0,variance(qcValues)))/m : NaN;
      }
      const unstable = Number.isFinite(rsd) && rsd > rsdThreshold;
      qcRsdStats.push({ feature, rsd, unstable });
      if (!unstable) keep.push(feature);
    }
    features = keep;
  }

  let mnarImputed = 0;
  const mnarByFeature = [];
  if (mnarStrategy === 'left_censored') {
    for (const feature of features) {
      const map = values.get(feature);
      const observed = biologicalAssays.map((assay)=>map.get(assay)).filter((value)=>Number.isFinite(value) && value > 0);
      if (observed.length < 3) continue;
      const logs = observed.map(Math.log);
      const robustSpread = 1.4826 * mad(logs);
      const fallbackSpread = Math.sqrt(Math.max(0,variance(logs)));
      const spread = Math.max(1e-6, Number.isFinite(robustSpread) && robustSpread > 0 ? robustSpread : (fallbackSpread || 1));
      const minPositive = Math.min(...observed);
      const q10 = simpleQuantile(logs,0.10);
      const low = Math.max(
        minPositive * 0.01,
        Math.min(minPositive * 0.8, Math.exp(Math.min(Math.log(minPositive)-0.2*spread, q10-1.28*spread)))
      );
      let count = 0;
      for (const assay of biologicalAssays) {
        const value = map.get(assay);
        if (!Number.isFinite(value)) {
          map.set(assay, low);
          count += 1;
          mnarImputed += 1;
        }
      }
      if (count) mnarByFeature.push({ feature, count, imputedValue: low });
    }
  }

  const warnings = [];
  if (blankFilterEnabled && blankAssays.length < 2) warnings.push('Blank assessment requested but fewer than 2 blank injections were annotated.');
  if (blankMode === 'remove' && blankAssays.length >= 2) warnings.push('Blank-associated features were excluded using the declared biological/blank ratio. Review flagged features and retain a flag-only sensitivity analysis for confirmatory work.');
  if (driftEnabled && orderedQcCount < 5) warnings.push('QC drift correction requested but fewer than 5 pooled-QC injections with numeric injection order were available.');
  if (driftEnabled && outsideQcSupport > 0 && orderedQcCount >= 5) warnings.push(outsideQcSupport + ' biological injections fall outside pooled-QC injection-order support; no drift extrapolation was applied to them.');
  if (qcRsdEnabled && qcAssays.length < 3) warnings.push('QC RSD filtering requested but fewer than 3 pooled-QC injections were annotated.');
  if (mnarStrategy === 'left_censored') warnings.push('Left-censored MNAR imputation was applied only to missing biological metabolomics values; perform a no-imputation sensitivity analysis for confirmatory work.');

  return {
    matrix: {
      ...matrix,
      features,
      assays: biologicalAssays,
      values: new Map(features.map((feature)=>[
        feature,
        new Map(biologicalAssays.map((assay)=>[assay,values.get(feature).get(assay)]))
      ]))
    },
    qc: {
      applied: blankAssays.length > 0 || qcAssays.length > 0 || mnarStrategy !== 'none',
      blankAssays: blankAssays.length,
      qcAssays: qcAssays.length,
      biologicalAssays: biologicalAssays.length,
      blankMode,
      blankFold,
      blankFlaggedFeatures: blankStats.filter((item)=>item.contaminated).length,
      blankRemovedFeatures: blankMode === 'remove' ? blankStats.filter((item)=>item.contaminated).length : 0,
      blankFilteredFeatures: blankMode === 'remove' ? blankStats.filter((item)=>item.contaminated).length : 0,
      qcRsdThreshold: rsdThreshold,
      qcRsdFilteredFeatures: qcRsdStats.filter((item)=>item.unstable).length,
      driftCorrection: {
        requested: driftEnabled,
        applied: driftCorrected > 0,
        correctedFeatures: driftCorrected,
        orderedQcCount,
        outsideQcSupport,
        interpolationOnly: true,
        medianAbsoluteLogCorrection: medianAbsoluteDriftLog,
        method: 'pooled-QC local linear log-intensity correction using injection order'
      },
      mnar: {
        strategy: mnarStrategy,
        imputedValues: mnarImputed,
        affectedFeatures: mnarByFeature.length,
        method: mnarStrategy === 'left_censored'
          ? 'deterministic low-tail left-censored imputation bounded below each feature minimum'
          : 'no imputation'
      },
      warnings,
      topBlankContaminants: blankStats.filter((item)=>item.contaminated)
        .sort((a,b)=>(b.blankMedian||0)-(a.blankMedian||0)).slice(0,20),
      topUnstableQcFeatures: qcRsdStats.filter((item)=>item.unstable)
        .sort((a,b)=>(b.rsd||0)-(a.rsd||0)).slice(0,20)
    }
  };
}

function prepareMatrixQc(matrix, layer, valueType, metadataRows = []) {
  const assays = matrix.assays;
  const featureCountBefore = matrix.features.length;
  const sampleMetrics = assays.map((assay) => {
    const raw = matrix.features.map((feature) => matrix.values.get(feature)?.get(assay));
    const observed = raw.filter(Number.isFinite);
    const positives = observed.filter((value) => value > 0);
    return {
      assayId: assay,
      observedFeatures: observed.length,
      detectedFeatures: positives.length,
      missingFraction: featureCountBefore ? 1 - observed.length / featureCountBefore : null,
      zeroFraction: observed.length ? observed.filter((value) => value === 0).length / observed.length : null,
      totalSignal: observed.reduce((sum, value) => sum + Math.max(0, value), 0),
      medianSignal: median(observed)
    };
  });

  const logTotals = sampleMetrics.map((metric) => Math.log1p(metric.totalSignal));
  const detected = sampleMetrics.map((metric) => metric.detectedFeatures);
  const totalOutliers = robustOutlierFlags(logTotals, 4);
  const detectionOutliers = robustOutlierFlags(detected, 4);
  sampleMetrics.forEach((metric, index) => {
    metric.outlier = totalOutliers[index] || detectionOutliers[index] || (metric.missingFraction != null && metric.missingFraction > 0.40);
    metric.outlierReasons = [
      totalOutliers[index] ? 'total signal' : '',
      detectionOutliers[index] ? 'detected feature count' : '',
      metric.missingFraction != null && metric.missingFraction > 0.40 ? 'missingness >40%' : ''
    ].filter(Boolean);
  });

  const groups = assayGroupsForQc(metadataRows, assays);
  const totals = new Map(sampleMetrics.map((metric) => [metric.assayId, metric.totalSignal || 1]));
  const isCounts = layer === 'transcriptomics' && valueType === 'raw_counts';
  const isSpectral = layer === 'proteomics' && valueType === 'spectral_count';

  const keep = [];
  const removalReasons = { lowAbundance: 0, excessiveMissingness: 0, constantOrEmpty: 0 };
  for (const feature of matrix.features) {
    const values = matrix.values.get(feature);
    const all = assays.map((assay) => values.get(assay));
    const finite = all.filter(Number.isFinite);
    if (!finite.length || new Set(finite).size <= 1) {
      removalReasons.constantOrEmpty += 1;
      continue;
    }

    if (isCounts || isSpectral) {
      let detectedInGroup = false;
      for (const groupAssays of groups.values()) {
        const required = Math.max(2, Math.ceil(groupAssays.length * 0.20));
        const detectedCount = groupAssays.filter((assay) => {
          const raw = values.get(assay);
          if (!Number.isFinite(raw)) return false;
          const cpm = raw / (totals.get(assay) || 1) * 1e6;
          return cpm >= 1;
        }).length;
        if (detectedCount >= Math.min(required, groupAssays.length)) {
          detectedInGroup = true;
          break;
        }
      }
      if (!detectedInGroup) {
        removalReasons.lowAbundance += 1;
        continue;
      }
    } else {
      let adequatelyObserved = false;
      for (const groupAssays of groups.values()) {
        const required = Math.max(2, Math.ceil(groupAssays.length * 0.50));
        const observed = groupAssays.filter((assay) => Number.isFinite(values.get(assay))).length;
        if (observed >= Math.min(required, groupAssays.length)) {
          adequatelyObserved = true;
          break;
        }
      }
      if (!adequatelyObserved) {
        removalReasons.excessiveMissingness += 1;
        continue;
      }
    }
    keep.push(feature);
  }

  const filteredValues = new Map(keep.map((feature) => [feature, matrix.values.get(feature)]));
  const differentialMissingness = [];
  if (groups.size >= 2) {
    for (const feature of keep) {
      const values = filteredValues.get(feature);
      const byGroup = {};
      for (const [group, groupAssays] of groups.entries()) {
        const n = groupAssays.length;
        const missing = groupAssays.filter((assay) => !Number.isFinite(values?.get(assay))).length;
        byGroup[group] = { n, missing, fraction: n ? missing / n : null };
      }
      const fractions = Object.values(byGroup).map((item) => item.fraction).filter(Number.isFinite);
      if (fractions.length < 2) continue;
      const imbalance = Math.max(...fractions) - Math.min(...fractions);
      if (imbalance >= 0.30) {
        differentialMissingness.push({ feature, imbalance, byGroup });
      }
    }
    differentialMissingness.sort((a,b) => b.imbalance - a.imbalance);
  }

  const replicateCorrelations = [];
  const bySample = new Map();
  for (const row of metadataRows) {
    if (!assays.includes(row.assayId)) continue;
    if (!bySample.has(row.sampleId)) bySample.set(row.sampleId, []);
    bySample.get(row.sampleId).push(row.assayId);
  }
  for (const [sampleId, replicateAssays] of bySample.entries()) {
    if (replicateAssays.length < 2) continue;
    for (let i = 0; i < replicateAssays.length; i += 1) {
      for (let j = i + 1; j < replicateAssays.length; j += 1) {
        const a = [];
        const b = [];
        for (const feature of keep) {
          const va = filteredValues.get(feature)?.get(replicateAssays[i]);
          const vb = filteredValues.get(feature)?.get(replicateAssays[j]);
          if (Number.isFinite(va) && Number.isFinite(vb)) {
            a.push(va);
            b.push(vb);
          }
        }
        const correlation = a.length >= 5 ? spearman(a,b) : NaN;
        replicateCorrelations.push({
          sampleId,
          assayA: replicateAssays[i],
          assayB: replicateAssays[j],
          nFeatures: a.length,
          correlation,
          warning: Number.isFinite(correlation) && correlation < 0.80
        });
      }
    }
  }

  const missingFractions = sampleMetrics.map((metric) => metric.missingFraction).filter(Number.isFinite);
  const warnings = [];
  const outlierSamples = sampleMetrics.filter((metric) => metric.outlier);
  if (outlierSamples.length) warnings.push(outlierSamples.length + ' assay(s) flagged by robust signal/detection/missingness QC.');
  const lowReplicatePairs = replicateCorrelations.filter((item) => item.warning);
  if (lowReplicatePairs.length) warnings.push(lowReplicatePairs.length + ' technical replicate pair(s) have Spearman r < 0.80.');
  if (differentialMissingness.length) warnings.push(differentialMissingness.length + ' retained feature(s) have ≥30 percentage-point missingness imbalance between biological conditions; inspect missingness before interpreting group effects.');
  if (keep.length < Math.max(1, featureCountBefore * 0.20)) warnings.push('More than 80% of features were removed by the conservative detection filter.');

  return {
    matrix: { ...matrix, features: keep, values: filteredValues },
    qc: {
      layer,
      valueType,
      featuresBefore: featureCountBefore,
      featuresAfter: keep.length,
      removedFeatures: featureCountBefore - keep.length,
      removalReasons,
      assays: assays.length,
      medianMissingFraction: median(missingFractions),
      sampleMetrics,
      outlierSamples,
      replicateCorrelations,
      differentialMissingness: {
        threshold: 0.30,
        flaggedFeatures: differentialMissingness.length,
        topFeatures: differentialMissingness.slice(0, 25),
        interpretation: 'QC flag only: differential observation/missingness can create or mask apparent biological differences and should be investigated before confirmatory interpretation.'
      },
      warnings,
      filterPolicy: isCounts || isSpectral
        ? 'retain non-constant features with CPM ≥1 in at least max(2, 20% of assays) within at least one biological condition'
        : 'retain non-constant features observed in at least max(2, 50% of assays) within at least one biological condition'
    }
  };
}

function minimumPositive(matrix) {
  let min = Infinity;
  for (const values of matrix.values.values()) {
    for (const value of values.values()) if (Number.isFinite(value) && value > 0 && value < min) min = value;
  }
  return Number.isFinite(min) ? min : 1;
}

function preprocessMatrix(matrix, layer, valueType) {
  const assays = matrix.assays;
  const processed = new Map();
  const steps = [];
  const finiteValues = [];
  for (const values of matrix.values.values()) {
    for (const value of values.values()) if (Number.isFinite(value)) finiteValues.push(value);
  }
  const minPositive = minimumPositive(matrix);
  const pseudo = Math.max(minPositive / 2, 1e-12);
  const negativeValues = finiteValues.filter((value) => value < 0).length;
  const zeroValues = finiteValues.filter((value) => value === 0).length;
  const nonIntegerValues = finiteValues.filter((value) => Math.abs(value - Math.round(value)) > 1e-8).length;
  const hasNegative = negativeValues > 0;

  const isCounts = layer === 'transcriptomics' && valueType === 'raw_counts';
  const isSpectral = layer === 'proteomics' && valueType === 'spectral_count';
  const explicitlyLog = ['log_expression','log_intensity','log_abundance'].includes(valueType);
  const asSupplied = ['normalized','unknown'].includes(valueType) || (hasNegative && !isCounts && !isSpectral && !explicitlyLog);
  const logPositive = ['tpm','lfq_intensity','peak_area','concentration'].includes(valueType);
  const medianCenter = (layer === 'proteomics' && valueType === 'lfq_intensity') || (layer === 'metabolomics' && valueType === 'peak_area');
  let scale = 'as_supplied';
  const auditWarnings = [];

  if ((isCounts || isSpectral) && negativeValues > 0) {
    throw new Error(`${layer}: ${valueType} was declared but ${negativeValues} negative value(s) were observed. Count-like inputs must be non-negative; correct the declared data type or the input matrix before analysis.`);
  }
  if ((isCounts || isSpectral) && nonIntegerValues > 0) {
    auditWarnings.push(`${valueType} was declared but ${nonIntegerValues} non-integer value(s) were observed. The browser can screen these values, but reference count-model assumptions and upstream quantification should be reviewed.`);
  }
  if (logPositive && negativeValues > 0) {
    auditWarnings.push(`${valueType} was declared as a positive abundance scale but ${negativeValues} negative value(s) were observed. The browser keeps the supplied values instead of forcing a log transform; verify whether the matrix was already transformed or centred.`);
  }
  if (valueType === 'unknown') {
    auditWarnings.push('The measurement scale is declared unknown. Values are kept as supplied; log2 fold-ratio interpretation is disabled unless the scale is clarified.');
  }

  const totals = new Map();
  if (isCounts || isSpectral) {
    scale = 'log2';
    for (const assay of assays) {
      let total = 0;
      for (const values of matrix.values.values()) {
        const v = values.get(assay);
        if (Number.isFinite(v) && v > 0) total += v;
      }
      totals.set(assay, total || 1);
    }
    steps.push(isCounts ? 'library-size normalisation to CPM + log2(CPM + 0.5)' : 'library-size normalisation + log2 transform');
  } else if (explicitlyLog) {
    scale = 'log2';
    steps.push('values treated as already log-transformed');
  } else if (asSupplied) {
    scale = hasNegative ? 'transformed_unknown' : 'as_supplied';
    steps.push(hasNegative
      ? 'negative values detected: values kept as supplied; effect is a difference on the supplied scale'
      : 'normalised/unknown values kept as supplied; effect is a difference on the supplied scale');
  } else if (logPositive) {
    scale = 'log2';
    steps.push(`log2(value + ${pseudo.toPrecision(3)})`);
  } else {
    scale = 'as_supplied';
    steps.push('values kept as supplied');
  }

  for (const [feature, values] of matrix.values.entries()) {
    const next = new Map();
    for (const assay of assays) {
      const raw = values.get(assay);
      if (!Number.isFinite(raw)) { next.set(assay, null); continue; }
      if (isCounts || isSpectral) {
        const cpm = raw / totals.get(assay) * 1e6;
        next.set(assay, Math.log2(cpm + 0.5));
      } else if (explicitlyLog || asSupplied || !logPositive) next.set(assay, raw);
      else next.set(assay, Math.log2(Math.max(0, raw) + pseudo));
    }
    processed.set(feature, next);
  }

  if (medianCenter) {
    const medians = new Map();
    for (const assay of assays) {
      const vals = [...processed.values()].map((values) => values.get(assay)).filter(Number.isFinite);
      medians.set(assay, median(vals));
    }
    const target = median([...medians.values()]);
    for (const values of processed.values()) {
      for (const assay of assays) {
        const v = values.get(assay);
        if (Number.isFinite(v)) values.set(assay, v - medians.get(assay) + target);
      }
    }
    steps.push('sample-wise median centering on the log scale');
  }

  const inputMin = finiteValues.length ? Math.min(...finiteValues) : null;
  const inputMax = finiteValues.length ? Math.max(...finiteValues) : null;
  const audit = {
    layer,
    declaredValueType: valueType,
    observedFiniteValues: finiteValues.length,
    negativeValues,
    zeroValues,
    nonIntegerValues,
    inputRange: { min: inputMin, max: inputMax },
    minimumPositive: Number.isFinite(minPositive) ? minPositive : null,
    pseudoCount: (logPositive && !asSupplied) ? pseudo : null,
    outputScale: scale,
    medianCentered: medianCenter,
    status: auditWarnings.length ? 'review_required' : 'compatible',
    warnings: auditWarnings,
    policy: 'Transformations follow the explicitly declared measurement type. The browser does not infer or silently apply ComBat, quantile normalisation, VST or other study-specific preprocessing.'
  };

  return { ...matrix, values: processed, steps, scale, audit };
}

/**
 * Prediction inputs cannot inherit cohort-wide feature filtering, centering,
 * pseudocount estimation or feature selection performed before nested CV.
 * Only explicitly declared per-observation/per-assay transformations are safe.
 * All feature scoring, nuisance fitting, missingness and scaling take place
 * on each training fold in selectPredictionFeatures / predictionMatrix.
 */
export function prepareFoldIsolatedPredictionMatrix(matrix, layer, valueType) {
  const asSupplied = ['log_expression','log_intensity','log_abundance'].includes(valueType);
  const rawRna = layer === 'transcriptomics' && valueType === 'raw_counts';
  // Unlike the exploratory QC path, these intensity-scale predictive recipes
  // do not estimate a cohort-wide pseudocount, abundance filter or median.
  // log2(x + 1) fixes the offset before outcome data are seen. For LFQ and
  // MS peak areas, sample median centering uses ONLY the current assay.
  // This deliberately does not pretend to reproduce pooled-QC drift correction.
  const sampleLocalIntensity =
    (layer === 'proteomics' && valueType === 'lfq_intensity') ||
    (layer === 'metabolomics' && valueType === 'peak_area');
  if (!asSupplied && !rawRna && !sampleLocalIntensity) return null;
  if (sampleLocalIntensity) {
    for (const feature of matrix.features) {
      for (const assay of matrix.assays) {
        const value = matrix.values.get(feature)?.get(assay);
        if (Number.isFinite(value) && value < 0) {
          // Already transformed or background-subtracted data were mislabeled.
          // Do not silently treat them as positive raw intensity.
          return null;
        }
      }
    }
  }
  const assayMedians = new Map();
  if (sampleLocalIntensity) {
    for (const assay of matrix.assays) {
      const logged = matrix.features.map((feature) => matrix.values.get(feature)?.get(assay))
        .filter(Number.isFinite).map((value) => Math.log2(value + 1));
      assayMedians.set(assay, logged.length ? median(logged) : 0);
    }
  }
  const totals = new Map();
  if (rawRna) {
    for (const assay of matrix.assays) {
      let total = 0;
      for (const values of matrix.values.values()) {
        const value = values.get(assay);
        if (Number.isFinite(value)) {
          if (value < 0 || !Number.isInteger(value)) {
            return null; // screening of fractional counts is not a validated predictive input
          }
          total += value;
        }
      }
      totals.set(assay, total);
    }
  }
  const values = new Map();
  for (const feature of matrix.features) {
    const row = new Map();
    for (const assay of matrix.assays) {
      const value = matrix.values.get(feature)?.get(assay);
      row.set(assay, !Number.isFinite(value) ? null
        : rawRna
          ? Math.log2((value / Math.max(1, totals.get(assay))) * 1e6 + 0.5)
          : sampleLocalIntensity
            ? Math.log2(value + 1) - assayMedians.get(assay)
            : value);
    }
    values.set(feature, row);
  }
  return {
    ...matrix,
    values,
    // No feature removal/selection and no learned cohort-level parameters.
    predictionPreparation: rawRna
      ? 'per-assay library-size CPM and fixed log2(CPM + 0.5), no global feature filtering'
      : sampleLocalIntensity
        ? 'fixed log2(raw intensity + 1), current-assay-only median centering, no cohort-wide feature filtering or QC drift fit'
        : 'as-supplied values, no global feature filtering'
  };
}

function aggregateTechnicalReplicates(matrix, metadata, layer) {
  const rows = metadata.filter((row) => row.omic === layer);
  const bySample = new Map();
  for (const row of rows) {
    if (!bySample.has(row.sampleId)) bySample.set(row.sampleId, []);
    bySample.get(row.sampleId).push(row.assayId);
  }
  const values = new Map();
  for (const feature of matrix.features) {
    const source = matrix.values.get(feature);
    const next = new Map();
    for (const [sampleId, assays] of bySample.entries()) {
      const vals = assays.map((assay) => source.get(assay)).filter(Number.isFinite);
      next.set(sampleId, vals.length ? mean(vals) : null);
    }
    values.set(feature, next);
  }
  const sampleMeta = new Map();
  for (const row of rows) {
    if (!sampleMeta.has(row.sampleId)) sampleMeta.set(row.sampleId, row);
  }
  const replicateGroups = [...bySample.entries()].filter(([,assays]) => assays.length > 1).map(([sampleId,assays]) => ({sampleId, assays}));
  return { features: matrix.features, values, sampleMeta, replicateGroups, steps: matrix.steps, scale: matrix.scale };
}


function solveLinearSystem(matrix, vector, ridge = 1e-8) {
  const n = matrix.length;
  if (!n || vector.length !== n) return null;
  const a = matrix.map((row, i) => row.map((value, j) => value + (i === j ? ridge : 0)).concat([vector[i]]));
  for (let col = 0; col < n; col += 1) {
    let pivot = col;
    for (let row = col + 1; row < n; row += 1) {
      if (Math.abs(a[row][col]) > Math.abs(a[pivot][col])) pivot = row;
    }
    if (Math.abs(a[pivot][col]) < 1e-12) return null;
    [a[col], a[pivot]] = [a[pivot], a[col]];
    const scale = a[col][col];
    for (let j = col; j <= n; j += 1) a[col][j] /= scale;
    for (let row = 0; row < n; row += 1) {
      if (row === col) continue;
      const factor = a[row][col];
      if (factor === 0) continue;
      for (let j = col; j <= n; j += 1) a[row][j] -= factor * a[col][j];
    }
  }
  return a.map((row) => row[n]);
}

function invertMatrix(matrix, ridge = 1e-8) {
  const n = matrix.length;
  if (!n) return null;
  const inverse = Array.from({ length: n }, () => Array(n).fill(0));
  for (let col = 0; col < n; col += 1) {
    const unit = Array(n).fill(0);
    unit[col] = 1;
    const solution = solveLinearSystem(matrix, unit, ridge);
    if (!solution) return null;
    for (let row = 0; row < n; row += 1) inverse[row][col] = solution[row];
  }
  return inverse;
}

function crossProductMatrix(design, weights = null) {
  const p = design[0]?.length || 0;
  const out = Array.from({ length: p }, () => Array(p).fill(0));
  for (let i = 0; i < design.length; i += 1) {
    const w = weights ? weights[i] : 1;
    for (let a = 0; a < p; a += 1) {
      for (let b = 0; b < p; b += 1) out[a][b] += w * design[i][a] * design[i][b];
    }
  }
  return out;
}

function crossProductVector(design, response, weights = null) {
  const p = design[0]?.length || 0;
  const out = Array(p).fill(0);
  for (let i = 0; i < design.length; i += 1) {
    const w = weights ? weights[i] : 1;
    for (let a = 0; a < p; a += 1) out[a] += w * design[i][a] * response[i];
  }
  return out;
}

function multiplyMatrixVector(matrix, vector) {
  return matrix.map((row) => row.reduce((sum, value, index) => sum + value * vector[index], 0));
}

function buildCovariateEncoder(rows, covariateColumns = [], { includeBatch = true } = {}) {
  const specs = [];
  if (includeBatch) {
    const levels = naturalOrder(rows.map((row) => row.batch || '__MISSING__'));
    if (levels.length > 1) specs.push({ name: 'batch', type: 'categorical', levels, reference: levels[0] });
  }

  for (const column of covariateColumns) {
    const values = rows.map((row) => row.covariates?.[column] ?? '').filter((value) => value !== '');
    if (!values.length) continue;
    const numericValues = values.map(finiteNumber);
    const numeric = numericValues.every(Number.isFinite) && new Set(numericValues).size > 1;
    if (numeric) {
      const center = mean(numericValues);
      const sd = Math.sqrt(variance(numericValues));
      specs.push({
        name: column,
        type: 'numeric',
        center,
        scale: Number.isFinite(sd) && sd > 0 ? sd : 1,
        median: median(numericValues)
      });
    } else {
      const levels = naturalOrder(rows.map((row) => row.covariates?.[column] || '__MISSING__'));
      if (levels.length > 1) specs.push({ name: column, type: 'categorical', levels, reference: levels[0] });
    }
  }

  const columnNames = ['intercept'];
  for (const spec of specs) {
    if (spec.type === 'numeric') {
      columnNames.push(spec.name);
      const missingName = spec.name + ':missing';
      const hasMissing = rows.some((row) => !Number.isFinite(finiteNumber(row.covariates?.[spec.name] ?? '')));
      if (hasMissing) columnNames.push(missingName);
    } else {
      for (const level of spec.levels.slice(1)) columnNames.push(spec.name + '=' + level);
    }
  }

  const encode = (row) => {
    const vector = [1];
    for (const spec of specs) {
      if (spec.type === 'numeric') {
        const parsed = finiteNumber(row.covariates?.[spec.name] ?? '');
        const missing = !Number.isFinite(parsed);
        const value = missing ? spec.median : parsed;
        vector.push((value - spec.center) / spec.scale);
        if (columnNames.includes(spec.name + ':missing')) vector.push(missing ? 1 : 0);
      } else {
        const value = spec.name === 'batch'
          ? (row.batch || '__MISSING__')
          : (row.covariates?.[spec.name] || '__MISSING__');
        for (const level of spec.levels.slice(1)) vector.push(value === level ? 1 : 0);
      }
    }
    return vector;
  };

  return { specs, columnNames, encode };
}

function residualizeAggregated(aggregated, covariateColumns = []) {
  const samples = [...aggregated.sampleMeta.entries()];
  const rows = samples.map(([, row]) => row);
  const encoder = buildCovariateEncoder(rows, covariateColumns, { includeBatch: true });
  if (encoder.columnNames.length === 1) {
    return {
      aggregated,
      adjustment: {
        applied: false,
        method: 'none',
        columns: [],
        covariates: covariateColumns,
        note: 'No varying batch or selected covariate required adjustment.'
      }
    };
  }

  const designAll = rows.map(encoder.encode);
  const values = new Map();
  let failedFeatures = 0;
  for (const feature of aggregated.features) {
    const source = aggregated.values.get(feature);
    const indexes = [];
    const response = [];
    for (let i = 0; i < samples.length; i += 1) {
      const sampleId = samples[i][0];
      const value = source.get(sampleId);
      if (!Number.isFinite(value)) continue;
      indexes.push(i);
      response.push(value);
    }
    const design = indexes.map((index) => designAll[index]);
    if (response.length <= encoder.columnNames.length + 1) {
      values.set(feature, new Map(source));
      failedFeatures += 1;
      continue;
    }
    const xtx = crossProductMatrix(design);
    const xty = crossProductVector(design, response);
    const beta = solveLinearSystem(xtx, xty);
    if (!beta) {
      values.set(feature, new Map(source));
      failedFeatures += 1;
      continue;
    }
    const fitted = multiplyMatrixVector(design, beta);
    const center = mean(response);
    const next = new Map(source);
    for (let k = 0; k < indexes.length; k += 1) {
      const sampleId = samples[indexes[k]][0];
      next.set(sampleId, response[k] - fitted[k] + center);
    }
    values.set(feature, next);
  }

  return {
    aggregated: { ...aggregated, values },
    adjustment: {
      applied: true,
      method: 'feature-wise OLS residualisation before biological contrast',
      columns: encoder.columnNames.slice(1),
      covariates: covariateColumns,
      failedFeatures,
      note: 'Batch and selected covariates are removed feature-wise before group, longitudinal or exploratory summaries. Complete confounding is blocked upstream.'
    }
  };
}

function subjectEndpointRows(aggregated, feature, requestedTimepoint = '') {
  const source = aggregated.values.get(feature);
  if (!source) return [];
  const timepoints = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.timepoint));
  const targetTime = requestedTimepoint || (timepoints.length ? timepoints[timepoints.length - 1] : '');
  const perSubject = new Map();
  for (const [sampleId, row] of aggregated.sampleMeta.entries()) {
    const value = source.get(sampleId);
    if (!Number.isFinite(value)) continue;
    if (targetTime && row.timepoint !== targetTime) continue;
    if (!perSubject.has(row.subjectId)) perSubject.set(row.subjectId, []);
    perSubject.get(row.subjectId).push({ value, row });
  }
  return [...perSubject.entries()].map(([subjectId, entries]) => ({
    subjectId,
    value: mean(entries.map((entry) => entry.value)),
    row: entries[0].row
  }));
}

function featureOutcomeRows(aggregated, feature, outcomeType, targetTimepoint = '') {
  return subjectEndpointRows(aggregated, feature, targetTimepoint).map(({ subjectId, value, row }) => ({
    subjectId,
    feature: value,
    outcome: row.outcome,
    survivalTime: finiteNumber(row.survivalTime),
    survivalEvent: row.survivalEvent,
    row
  })).filter((entry) => {
    if (outcomeType === 'survival') {
      const event = finiteNumber(entry.survivalEvent);
      return Number.isFinite(entry.feature) && Number.isFinite(entry.survivalTime) && Number.isFinite(event);
    }
    return Number.isFinite(entry.feature) && entry.outcome !== '';
  });
}

function subjectCovariateDesign(entries, covariateColumns = [], { includeIntercept = true, includeBatch = false } = {}) {
  const rows = entries.map((entry) => entry.row);
  const encoder = buildCovariateEncoder(rows, covariateColumns, { includeBatch });
  const columnNames = encoder.columnNames.slice(includeIntercept ? 0 : 1);
  const design = rows.map((row) => {
    const vector = encoder.encode(row);
    return includeIntercept ? vector : vector.slice(1);
  });
  return { design, columnNames, encoder };
}

// Inference must not use the tiny numerical ridge added by solveLinearSystem
// to disguise exact aliases (e.g. a technical batch perfectly coding group).
// This rank check uses scale-aware Gaussian elimination on the *raw* design.
function hasFullColumnRank(design) {
  const n = design.length;
  const p = design[0]?.length || 0;
  if (!n || !p || n < p || design.some((row) => row.length !== p || row.some((x) => !Number.isFinite(x)))) return false;
  const a = design.map((row) => row.slice());
  let rank = 0;
  for (let col = 0; col < p; col += 1) {
    let pivot = -1;
    let best = 0;
    for (let i = rank; i < n; i += 1) {
      const magnitude = Math.abs(a[i][col]);
      if (magnitude > best) { best = magnitude; pivot = i; }
    }
    const columnMax = Math.max(...design.map((row) => Math.abs(row[col])));
    if (pivot < 0 || best <= 1e-10 * Math.max(columnMax, 1e-12)) continue;
    [a[rank], a[pivot]] = [a[pivot], a[rank]];
    const scale = a[rank][col];
    for (let j = col; j < p; j += 1) a[rank][j] /= scale;
    for (let i = rank + 1; i < n; i += 1) {
      const factor = a[i][col];
      for (let j = col; j < p; j += 1) a[i][j] -= factor * a[rank][j];
    }
    rank += 1;
  }
  return rank === p;
}

function ordinaryLeastSquares(design, response, coefficientIndex) {
  if (!design.length || design.length <= (design[0]?.length || 0)) return null;
  if (!hasFullColumnRank(design)) return null;
  const xtx = crossProductMatrix(design);
  const xty = crossProductVector(design, response);
  const beta = solveLinearSystem(xtx, xty);
  const inverse = invertMatrix(xtx);
  if (!beta || !inverse) return null;
  const fitted = multiplyMatrixVector(design, beta);
  const residuals = response.map((value, i) => value - fitted[i]);
  const df = response.length - design[0].length;
  const rss = residuals.reduce((sum, value) => sum + value * value, 0);
  const sigma2 = df > 0 ? rss / df : NaN;
  const varianceBeta = Number.isFinite(sigma2) ? sigma2 * inverse[coefficientIndex][coefficientIndex] : NaN;
  const se = varianceBeta > 0 ? Math.sqrt(varianceBeta) : NaN;
  const t = Number.isFinite(se) && se > 0 ? beta[coefficientIndex] / se : NaN;
  const cdf = Number.isFinite(t) ? studentTCdf(Math.abs(t), df) : NaN;
  const pValue = Number.isFinite(cdf) ? Math.max(0, Math.min(1, 2 * (1 - cdf))) : null;
  return { beta, se, statistic: t, pValue, df, fitted, residuals };
}

function ordinaryLeastSquaresRobust(design, response, coefficientIndex) {
  const base = ordinaryLeastSquares(design, response, coefficientIndex);
  if (!base) return null;
  const xtx = crossProductMatrix(design);
  const bread = invertMatrix(xtx);
  if (!bread) return base;
  const p = design[0].length;
  const meat = Array.from({ length: p }, () => Array(p).fill(0));
  for (let i = 0; i < design.length; i += 1) {
    const x = design[i];
    const bx = multiplyMatrixVector(bread, x);
    const leverage = Math.max(0, Math.min(0.999999, x.reduce((sum, value, j) => sum + value * bx[j], 0)));
    const adjustedResidual = base.residuals[i] / Math.max(1e-6, 1 - leverage);
    const weight = adjustedResidual * adjustedResidual;
    for (let a = 0; a < p; a += 1) {
      for (let b = 0; b < p; b += 1) meat[a][b] += weight * x[a] * x[b];
    }
  }
  const temp = Array.from({ length: p }, () => Array(p).fill(0));
  const sandwich = Array.from({ length: p }, () => Array(p).fill(0));
  for (let i = 0; i < p; i += 1) {
    for (let j = 0; j < p; j += 1) {
      for (let k = 0; k < p; k += 1) temp[i][j] += bread[i][k] * meat[k][j];
    }
  }
  for (let i = 0; i < p; i += 1) {
    for (let j = 0; j < p; j += 1) {
      for (let k = 0; k < p; k += 1) sandwich[i][j] += temp[i][k] * bread[k][j];
    }
  }
  const varianceBeta = sandwich[coefficientIndex][coefficientIndex];
  const se = varianceBeta > 0 ? Math.sqrt(varianceBeta) : NaN;
  const t = Number.isFinite(se) && se > 0 ? base.beta[coefficientIndex] / se : NaN;
  const cdf = Number.isFinite(t) ? studentTCdf(Math.abs(t), base.df) : NaN;
  const pValue = Number.isFinite(cdf) ? Math.max(0, Math.min(1, 2 * (1 - cdf))) : null;
  return { ...base, se, statistic: t, pValue, robust: 'HC3' };
}

function fitGeneralizedLinear(design, response, family, coefficientIndex) {
  const n = response.length;
  const p = design[0]?.length || 0;
  if (!n || n <= p || !hasFullColumnRank(design)) return null;
  let beta = Array(p).fill(0);
  if (family === 'poisson') {
    const start = Math.log(Math.max(mean(response), 1e-6));
    beta[0] = Number.isFinite(start) ? start : 0;
  }
  let information = null;
  for (let iter = 0; iter < 40; iter += 1) {
    const eta = multiplyMatrixVector(design, beta);
    const weights = [];
    const z = [];
    for (let i = 0; i < n; i += 1) {
      if (family === 'binomial') {
        const e = Math.max(-30, Math.min(30, eta[i]));
        const m = 1 / (1 + Math.exp(-e));
        const w = Math.max(m * (1 - m), 1e-6);
        weights.push(w);
        z.push(eta[i] + (response[i] - m) / w);
      } else {
        const m = Math.max(1e-8, Math.min(1e8, Math.exp(Math.max(-20, Math.min(20, eta[i])))));
        const w = Math.max(m, 1e-6);
        weights.push(w);
        z.push(eta[i] + (response[i] - m) / w);
      }
    }
    const xtwx = crossProductMatrix(design, weights);
    const xtwz = crossProductVector(design, z, weights);
    const next = solveLinearSystem(xtwx, xtwz, 1e-6);
    if (!next) return null;
    const delta = Math.max(...next.map((value, index) => Math.abs(value - beta[index])));
    beta = next;
    information = xtwx;
    if (delta < 1e-7) break;
  }
  const inverse = information ? invertMatrix(information, 1e-6) : null;
  if (!inverse) return null;
  const se2 = inverse[coefficientIndex][coefficientIndex];
  const se = se2 > 0 ? Math.sqrt(se2) : NaN;
  const zStat = Number.isFinite(se) && se > 0 ? beta[coefficientIndex] / se : NaN;
  const pValue = Number.isFinite(zStat) ? Math.min(1, 2 * (1 - normalCdf(Math.abs(zStat)))) : null;
  return { beta, se, statistic: zStat, pValue };
}


/**
 * Breslow Cox partial-likelihood score and observed information.
 * All events at the SAME time share one risk-set denominator:
 * U(t) = sum_event(x) - d(t)*E_risk[x].
 * The former browser implementation iterated each event with its own
 * denominator, accidentally applying the d(t)-fold term while revisiting
 * identical risk sets. It was algebraically Breslow-like but lacked an
 * explicit, auditable grouped-ties contract.
 *
 * Risk weights are shifted by max linear predictor to avoid exp overflow.
 */
function coxBreslowScoreInformation(design, times, events, beta) {
  const n=times.length,p=beta.length;
  const deathTimes=[...new Set(times.filter((t,i)=>events[i]===1))].sort((a,b)=>a-b);
  const score=Array(p).fill(0);
  const information=Array.from({length:p},()=>Array(p).fill(0));
  for(const time of deathTimes) {
    const deaths=[];
    const risk=[];
    for(let i=0;i<n;i+=1) {
      if(times[i]>=time)risk.push(i);
      if(events[i]===1&&times[i]===time)deaths.push(i);
    }
    if(!deaths.length||!risk.length)return null;
    const eta=risk.map(i=>design[i].reduce((v,x,j)=>v+x*beta[j],0));
    const pivot=Math.max(...eta);
    const weights=eta.map(v=>Math.exp(v-pivot));
    const denom=weights.reduce((x,y)=>x+y,0);
    if(!(denom>0)||!Number.isFinite(denom))return null;
    const mu=Array(p).fill(0);
    const second=Array.from({length:p},()=>Array(p).fill(0));
    for(let k=0;k<risk.length;k+=1) {
      const x=design[risk[k]],w=weights[k]/denom;
      for(let j=0;j<p;j+=1) {
        mu[j]+=w*x[j];
        for(let h=0;h<p;h+=1)second[j][h]+=w*x[j]*x[h];
      }
    }
    for(const idx of deaths) {
      for(let j=0;j<p;j+=1)score[j]+=design[idx][j];
    }
    for(let j=0;j<p;j+=1) {
      score[j]-=deaths.length*mu[j];
      for(let h=0;h<p;h+=1) {
        information[j][h]+=deaths.length*(second[j][h]-mu[j]*mu[h]);
      }
    }
  }
  if(score.some(x=>!Number.isFinite(x))||
      information.some(row=>row.some(x=>!Number.isFinite(x))))return null;
  return {score,information,deathTimeGroups:deathTimes.length};
}

function validCoxInput(design,times,events) {
  const n=times.length,p=design[0]?.length||0;
  return n>p&&p>0&&events.length===n&&
    events.reduce((sum,x)=>sum+x,0)>=3&&
    times.every(x=>Number.isFinite(x)&&x>0)&&
    events.every(x=>x===0||x===1)&&
    design.length===n&&design.every(row=>row.length===p&&row.every(Number.isFinite));
}

export function coxRegression(design,times,events,coefficientIndex=0) {
  if(!validCoxInput(design,times,events)||
    !hasFullColumnRank(design)||
    !Number.isInteger(coefficientIndex)||coefficientIndex<0||
    coefficientIndex>=design[0].length)return null;
  const p=design[0].length;
  let beta=Array(p).fill(0);
  let converged=false;
  for(let iter=0;iter<70;iter+=1) {
    const parts=coxBreslowScoreInformation(design,times,events,beta);
    if(!parts)return null;
    const step=solveLinearSystem(parts.information,parts.score,1e-8);
    if(!step||step.some(x=>!Number.isFinite(x)))return null;
    const maxStep=Math.max(...step.map(Math.abs));
    const fraction=Math.min(1,3/Math.max(1e-12,maxStep));
    beta=beta.map((v,j)=>v+fraction*step[j]);
    if(maxStep*fraction<1e-8){converged=true;break;}
    if(beta.some(x=>!Number.isFinite(x)||Math.abs(x)>100))return null;
  }
  if(!converged)return null; // no inferential Wald statistics from a failed fit
  const parts=coxBreslowScoreInformation(design,times,events,beta);
  if(!parts)return null;
  const inv=invertMatrix(parts.information,1e-10);
  if(!inv)return null;
  const v=inv[coefficientIndex][coefficientIndex];
  const se=v>0?Math.sqrt(v):NaN;
  const statistic=Number.isFinite(se)&&se>0?beta[coefficientIndex]/se:NaN;
  const pValue=Number.isFinite(statistic)
    ?Math.max(0,Math.min(1,2*(1-normalCdf(Math.abs(statistic))))):null;
  return Number.isFinite(statistic)
    ?{beta,se,statistic,pValue,ties:'breslow',converged:true,
       deathTimeGroups:parts.deathTimeGroups}
    :null;
}

function coxRidgeRegression(design,times,events,lambda=1) {
  if(!validCoxInput(design,times,events)||!(Number.isFinite(lambda)&&lambda>=0))return null;
  const p=design[0].length;
  let beta=Array(p).fill(0);
  for(let iter=0;iter<70;iter+=1) {
    const parts=coxBreslowScoreInformation(design,times,events,beta);
    if(!parts)return null;
    const information=parts.information.map((row,i)=>row.map((value,j)=>
      value+(i===j?lambda:0)));
    const score=parts.score.map((value,j)=>value-lambda*beta[j]);
    const step=solveLinearSystem(information,score,1e-8);
    if(!step||step.some(x=>!Number.isFinite(x)))return null;
    const maxStep=Math.max(...step.map(Math.abs));
    const fraction=Math.min(1,3/Math.max(1e-12,maxStep));
    beta=beta.map((v,j)=>v+fraction*step[j]);
    if(maxStep*fraction<1e-8)return {beta,lambda,ties:'breslow',converged:true};
    if(beta.some(x=>!Number.isFinite(x)||Math.abs(x)>100))return null;
  }
  return null; // no model if ridge fit did not converge
}

function harrellCIndex(times, events, risks) {
  let comparable = 0;
  let concordant = 0;
  for (let i = 0; i < times.length; i += 1) {
    for (let j = i + 1; j < times.length; j += 1) {
      let early = -1;
      let late = -1;
      if (times[i] < times[j] && events[i] === 1) {
        early = i; late = j;
      } else if (times[j] < times[i] && events[j] === 1) {
        early = j; late = i;
      } else {
        continue;
      }
      comparable += 1;
      if (risks[early] > risks[late]) concordant += 1;
      else if (risks[early] === risks[late]) concordant += 0.5;
    }
  }
  return comparable ? concordant / comparable : NaN;
}

function fitRandomInterceptGls(design, response, subjectIds, coefficientIndex) {
  const n = response.length;
  const p = design[0]?.length || 0;
  const subjects = [...new Set(subjectIds)];
  if (!n || n <= p || subjects.length < 3 || !hasFullColumnRank(design)) return null;

  let beta = solveLinearSystem(crossProductMatrix(design), crossProductVector(design, response));
  if (!beta) return null;
  let sigmaWithin = 1;
  let sigmaBetween = 0;

  for (let iter = 0; iter < 6; iter += 1) {
    const fitted = multiplyMatrixVector(design, beta);
    const residuals = response.map((value, i) => value - fitted[i]);
    const groups = new Map();
    for (let i = 0; i < n; i += 1) {
      if (!groups.has(subjectIds[i])) groups.set(subjectIds[i], []);
      groups.get(subjectIds[i]).push(residuals[i]);
    }

    let withinSs = 0;
    let withinDf = 0;
    const means = [];
    const invSizes = [];
    for (const values of groups.values()) {
      const m = mean(values);
      means.push(m);
      invSizes.push(1 / values.length);
      for (const value of values) withinSs += (value - m) * (value - m);
      withinDf += Math.max(0, values.length - 1);
    }
    sigmaWithin = withinDf > 0 ? Math.max(withinSs / withinDf, 1e-8) : Math.max(variance(residuals), 1e-8);
    const meanVariance = variance(means);
    sigmaBetween = Number.isFinite(meanVariance)
      ? Math.max(0, meanVariance - sigmaWithin * mean(invSizes))
      : 0;

    const xtvx = Array.from({ length: p }, () => Array(p).fill(0));
    const xtvy = Array(p).fill(0);
    for (const subject of subjects) {
      const idx = subjectIds.map((id, i) => id === subject ? i : -1).filter((i) => i >= 0);
      const m = idx.length;
      const a = 1 / sigmaWithin;
      const b = sigmaBetween > 0
        ? sigmaBetween / (sigmaWithin * (sigmaWithin + m * sigmaBetween))
        : 0;
      for (const ii of idx) {
        for (let col = 0; col < p; col += 1) {
          let vy = a * response[ii];
          let sumY = 0;
          for (const jj of idx) sumY += response[jj];
          vy -= b * sumY;
          xtvy[col] += design[ii][col] * vy;
          for (let col2 = 0; col2 < p; col2 += 1) {
            let vx = a * design[ii][col2];
            let sumX = 0;
            for (const jj of idx) sumX += design[jj][col2];
            vx -= b * sumX;
            xtvx[col][col2] += design[ii][col] * vx;
          }
        }
      }
    }
    const next = solveLinearSystem(xtvx, xtvy);
    if (!next) return null;
    const delta = Math.max(...next.map((value, i) => Math.abs(value - beta[i])));
    beta = next;
    if (delta < 1e-8) break;
  }

  const xtvx = Array.from({ length: p }, () => Array(p).fill(0));
  for (const subject of subjects) {
    const idx = subjectIds.map((id, i) => id === subject ? i : -1).filter((i) => i >= 0);
    const m = idx.length;
    const a = 1 / sigmaWithin;
    const b = sigmaBetween > 0
      ? sigmaBetween / (sigmaWithin * (sigmaWithin + m * sigmaBetween))
      : 0;
    for (const ii of idx) {
      for (let col = 0; col < p; col += 1) {
        for (let col2 = 0; col2 < p; col2 += 1) {
          let vx = a * design[ii][col2];
          let sumX = 0;
          for (const jj of idx) sumX += design[jj][col2];
          vx -= b * sumX;
          xtvx[col][col2] += design[ii][col] * vx;
        }
      }
    }
  }
  const covariance = invertMatrix(xtvx);
  if (!covariance) return null;
  const varianceBeta = covariance[coefficientIndex][coefficientIndex];
  const se = varianceBeta > 0 ? Math.sqrt(varianceBeta) : NaN;
  const statistic = Number.isFinite(se) && se > 0 ? beta[coefficientIndex] / se : NaN;
  const df = Math.max(1, subjects.length - 2);
  const cdf = Number.isFinite(statistic) ? studentTCdf(Math.abs(statistic), df) : NaN;
  const pValue = Number.isFinite(cdf) ? Math.max(0, Math.min(1, 2 * (1 - cdf))) : null;
  return {
    beta,
    se,
    statistic,
    pValue,
    df,
    sigmaWithin,
    sigmaBetween,
    covariance,
    intraclassCorrelation: sigmaBetween / (sigmaBetween + sigmaWithin)
  };
}

function analyseLongitudinalMixedLayer(aggregated, layer, { covariateColumns = [] } = {}) {
  const rows = [];
  const sampleRows = [...aggregated.sampleMeta.entries()].map(([sampleId, row]) => ({ sampleId, ...row }));
  const conditions = naturalOrder(sampleRows.map((row) => row.condition));
  const times = naturalOrder(sampleRows.map((row) => row.timepoint));
  const numericTimes = times.map(numericTime);
  if (conditions.length < 2) return { error: 'Longitudinal random-intercept model requires at least two conditions.', rows: [], selected: [], groupSizes: [] };
  if (numericTimes.some((value) => !Number.isFinite(value))) return { error: 'Longitudinal random-intercept model requires numeric or numeric-labelled time points.', rows: [], selected: [], groupSizes: [] };

  const reference = conditions[0];
  const comparisons = conditions.slice(1);
  const minTime = Math.min(...numericTimes);
  const maxTime = Math.max(...numericTimes);
  const timeRange = Math.max(1e-12, maxTime - minTime);
  const multiGroup = conditions.length > 2;

  for (const feature of aggregated.features) {
    const source = aggregated.values.get(feature);
    const entries = sampleRows
      .map((row) => ({ row, value: source.get(row.sampleId), time: numericTime(row.timepoint) }))
      .filter((entry) => Number.isFinite(entry.value) && Number.isFinite(entry.time));
    const subjectCounts = new Map();
    for (const entry of entries) subjectCounts.set(entry.row.subjectId, (subjectCounts.get(entry.row.subjectId) || 0) + 1);
    const usable = entries.filter((entry) => subjectCounts.get(entry.row.subjectId) >= 2);
    const subjectIds = usable.map((entry) => entry.row.subjectId);
    const subjectN = new Set(subjectIds).size;
    if (subjectN < Math.max(4, conditions.length + 1)) continue;
    if (conditions.some((condition) => !usable.some((entry) => entry.row.condition === condition))) continue;

    const nuisance = buildCovariateEncoder(usable.map((entry) => entry.row), covariateColumns, { includeBatch: true });
    const conditionColumns = comparisons.length;
    const timeIndex = 1 + conditionColumns;
    const interactionStart = timeIndex + 1;
    const design = usable.map((entry) => {
      const time = entry.time - minTime;
      const conditionDummies = comparisons.map((condition) => entry.row.condition === condition ? 1 : 0);
      const interactions = conditionDummies.map((dummy) => dummy * time);
      return [1, ...conditionDummies, time, ...interactions, ...nuisance.encode(entry.row).slice(1)];
    });
    const response = usable.map((entry) => entry.value);
    const fit = fitRandomInterceptGls(design, response, subjectIds, interactionStart);
    if (!fit) continue;

    if (!multiGroup) {
      const interaction = fit.beta[interactionStart];
      const effect = interaction * timeRange;
      const effectSe = fit.se * timeRange;
      rows.push({
        feature,
        effect,
        foldRatio: aggregated.scale === 'log2' ? Math.pow(2, effect) : null,
        effectScale: aggregated.scale,
        // Browser GLS random-intercept Wald inference exceeded the nominal
        // null rejection rate in predeclared simulations. Retain only the
        // descriptive estimated effect until lmerTest confirms uncertainty.
        pValue: null,
        qValue: null,
        statistic: null,
        standardError: null,
        ciLow: null,
        ciHigh: null,
        inferentialStatus: 'uncalibrated_reference_R_required',
        intraclassCorrelation: fit.intraclassCorrelation,
        sigmaWithin: fit.sigmaWithin,
        sigmaBetween: fit.sigmaBetween,
        nSubjects: subjectN,
        nObservations: response.length,
        nReference: new Set(usable.filter((entry) => entry.row.condition === reference).map((entry) => entry.row.subjectId)).size,
        nComparison: new Set(usable.filter((entry) => entry.row.condition === comparisons[0]).map((entry) => entry.row.subjectId)).size,
        model: 'EXPLORATORY random-intercept GLS effect (browser p/q suppressed): feature ~ condition * time + batch + selected covariates + (1|subject)'
      });
      continue;
    }

    const interactionIndices = comparisons.map((_, index) => interactionStart + index);
    const interactionBeta = interactionIndices.map((index) => fit.beta[index]);
    const interactionCovariance = interactionIndices.map((rowIndex) =>
      interactionIndices.map((columnIndex) => fit.covariance[rowIndex][columnIndex])
    );
    const inverseInteractionCovariance = invertMatrix(interactionCovariance, 1e-10);
    if (!inverseInteractionCovariance) continue;
    const weighted = multiplyMatrixVector(inverseInteractionCovariance, interactionBeta);
    const wald = Math.max(0, interactionBeta.reduce((sum, value, index) => sum + value * weighted[index], 0));
    const df1 = interactionIndices.length;
    const df2 = Math.max(1, subjectN - conditions.length);
    const fStatistic = wald / Math.max(1, df1);
    const pValue = Math.max(0, Math.min(1, 1 - fCdf(fStatistic, df1, df2)));

    const referenceSlope = fit.beta[timeIndex];
    const slopes = [
      referenceSlope,
      ...interactionBeta.map((interaction) => referenceSlope + interaction)
    ];
    const endToEndChanges = slopes.map((slope) => slope * timeRange);
    const effect = Math.max(...endToEndChanges) - Math.min(...endToEndChanges);

    rows.push({
      feature,
      effect,
      foldRatio: aggregated.scale === 'log2' ? Math.pow(2, effect) : null,
      effectScale: aggregated.scale,
      pValue: null,
      qValue: null,
      statistic: null,
      testDf: { numerator: df1, denominator: df2 },
      inferentialStatus: 'uncalibrated_reference_R_required',
      groupSlopes: Object.fromEntries(conditions.map((condition, index) => [condition, slopes[index]])),
      groupEndToEndChanges: Object.fromEntries(conditions.map((condition, index) => [condition, endToEndChanges[index]])),
      intraclassCorrelation: fit.intraclassCorrelation,
      sigmaWithin: fit.sigmaWithin,
      sigmaBetween: fit.sigmaBetween,
      nSubjects: subjectN,
      nObservations: response.length,
      groupSizes: Object.fromEntries(conditions.map((condition) => [
        condition,
        new Set(usable.filter((entry) => entry.row.condition === condition).map((entry) => entry.row.subjectId)).size
      ])),
      nReference: null,
      nComparison: null,
      model: 'EXPLORATORY random-intercept GLS omnibus effect (browser p/q suppressed): feature ~ condition * time + batch + selected covariates + (1|subject)'
    });
  }

  // Deliberately no BH correction for unsupported browser Wald p-values.
  // The reference R backend must supply inferential statistics separately.
  bhAdjust(rows);
  rows.sort((a,b) => {
    const aq = Number.isFinite(a.qValue) ? a.qValue : 1;
    const bq = Number.isFinite(b.qValue) ? b.qValue : 1;
    if (aq !== bq) return aq - bq;
    return Math.abs(b.effect) - Math.abs(a.effect);
  });
  const significant = rows.filter((row) =>
    Number.isFinite(row.qValue) && row.qValue <= 0.10 &&
    (aggregated.scale !== 'log2' || Math.abs(row.effect) >= Math.log2(1.2))
  );
  const selected = (significant.length >= 10 ? significant : rows).slice(0, significant.length >= 10 ? 50 : 25);
  const groupSizes = conditions.map((condition) =>
    new Set(sampleRows.filter((row) => row.condition === condition).map((row) => row.subjectId)).size
  );

  return {
    rows,
    selected,
    selectionRule: significant.length >= 10
      ? 'q ≤ 0.10 with longitudinal condition × time effect, capped at 50 features'
      : 'top ranked longitudinal interaction features retained for exploratory pathway mapping',
    effectScale: aggregated.scale,
    contrast: multiGroup
      ? 'omnibus condition × time interaction across ' + conditions.length + ' groups over ' + minTime + '→' + maxTime
      : 'condition × time interaction: ' + comparisons[0] + ' vs ' + reference + ' over ' + minTime + '→' + maxTime,
    mode: multiGroup ? 'random-intercept-longitudinal-omnibus' : 'random-intercept-longitudinal-model',
    inferenceMethod: 'Exploratory iterative random-intercept GLS effect ONLY; browser Wald/t p-values and q-values suppressed after null-calibration concern. Use lmerTest in reference R backend for confirmatory longitudinal inference.',
    inferenceStatus: 'reference_required_for_p_q_and_ci',
    referenceMethodRequired: 'lmerTest',
    groupSizes,
    steps: aggregated.steps
  };
}
function analyseIndependentAdjustedLayer(aggregated, layer, { covariateColumns = [] } = {}) {
  const rows = [];
  const allConditions = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.condition));
  if (allConditions.length < 2) {
    return { error: 'At least two conditions are required.', rows: [], selected: [], groupSizes: [] };
  }

  for (const feature of aggregated.features) {
    const entries = subjectEndpointRows(aggregated, feature);
    const conditions = naturalOrder(entries.map((entry) => entry.row.condition));
    if (conditions.length !== allConditions.length) continue;
    const nuisance = subjectCovariateDesign(entries, covariateColumns, { includeIntercept: true, includeBatch: true });
    const response = entries.map((entry) => entry.value);

    if (conditions.length === 2) {
      const reference = conditions[0];
      const comparison = conditions[1];
      const design = entries.map((entry, i) => [
        1,
        entry.row.condition === comparison ? 1 : 0,
        ...nuisance.design[i].slice(1)
      ]);
      const fit = ordinaryLeastSquaresRobust(design, response, 1);
      if (!fit) continue;
      const effect = fit.beta[1];
      rows.push({
        feature,
        effect,
        foldRatio: aggregated.scale === 'log2' ? Math.pow(2, effect) : null,
        effectScale: aggregated.scale,
        pValue: fit.pValue,
        qValue: null,
        statistic: fit.statistic,
        standardError: fit.se,
        ciLow: Number.isFinite(fit.se) ? effect - tCritical95(fit.df) * fit.se : null,
        ciHigh: Number.isFinite(fit.se) ? effect + tCritical95(fit.df) * fit.se : null,
        nReference: entries.filter((entry) => entry.row.condition === reference).length,
        nComparison: entries.filter((entry) => entry.row.condition === comparison).length,
        model: 'feature ~ condition + batch + selected covariates (OLS HC3)'
      });
    } else {
      const classLevels = conditions.slice(1);
      const reducedDesign = nuisance.design;
      const fullDesign = entries.map((entry, i) => [
        1,
        ...classLevels.map((group) => entry.row.condition === group ? 1 : 0),
        ...nuisance.design[i].slice(1)
      ]);
      const reduced = ordinaryLeastSquares(reducedDesign, response, 0);
      const full = ordinaryLeastSquares(fullDesign, response, 0);
      if (!reduced || !full) continue;
      const rssReduced = reduced.residuals.reduce((sum, value) => sum + value * value, 0);
      const rssFull = full.residuals.reduce((sum, value) => sum + value * value, 0);
      const df1 = classLevels.length;
      const df2 = response.length - fullDesign[0].length;
      if (!(df1 > 0) || !(df2 > 0)) continue;
      const fStatistic = (rssFull / df2) > 0 ? Math.max(0, (rssReduced - rssFull) / df1) / (rssFull / df2) : Infinity;
      const pValue = Number.isFinite(fStatistic) || fStatistic === Infinity
        ? Math.max(0, Math.min(1, 1 - fCdf(fStatistic, df1, df2)))
        : null;
      const adjustedMeans = [full.beta[0], ...classLevels.map((_, index) => full.beta[0] + full.beta[index + 1])];
      rows.push({
        feature,
        effect: Math.max(...adjustedMeans) - Math.min(...adjustedMeans),
        foldRatio: aggregated.scale === 'log2' ? Math.pow(2, Math.max(...adjustedMeans) - Math.min(...adjustedMeans)) : null,
        effectScale: aggregated.scale,
        pValue,
        qValue: null,
        statistic: fStatistic,
        groupMeans: Object.fromEntries(conditions.map((condition, index) => [condition, adjustedMeans[index]])),
        groupSizes: Object.fromEntries(conditions.map((condition) => [condition, entries.filter((entry) => entry.row.condition === condition).length])),
        nReference: null,
        nComparison: null,
        model: 'feature ~ condition + batch + selected covariates (ANCOVA partial F-test)'
      });
    }
  }

  bhAdjust(rows);
  rows.sort((a,b) => {
    const aq = Number.isFinite(a.qValue) ? a.qValue : 1;
    const bq = Number.isFinite(b.qValue) ? b.qValue : 1;
    if (aq !== bq) return aq - bq;
    return Math.abs(b.effect) - Math.abs(a.effect);
  });
  const significant = rows.filter((row) =>
    Number.isFinite(row.qValue) &&
    row.qValue <= 0.10 &&
    (aggregated.scale !== 'log2' || Math.abs(row.effect) >= Math.log2(1.2))
  );
  const selected = (significant.length >= 10 ? significant : rows).slice(0, significant.length >= 10 ? 50 : 25);
  return {
    rows,
    selected,
    selectionRule: significant.length >= 10
      ? (aggregated.scale === 'log2' ? 'q ≤ 0.10 and |fold change| ≥ 1.2, capped at 50 features' : 'q ≤ 0.10, capped at 50 features')
      : 'top ranked adjusted features retained for exploratory biological mapping because fewer than 10 features passed the FDR/effect threshold',
    effectScale: aggregated.scale,
    contrast: allConditions.length === 2
      ? allConditions[1] + ' vs ' + allConditions[0] + ' adjusted for batch/covariates'
      : 'omnibus adjusted condition effect across ' + allConditions.length + ' groups',
    mode: allConditions.length === 2 ? 'adjusted-two-group-model' : 'adjusted-multi-group-model',
    inferenceMethod: allConditions.length === 2 ? 'OLS with HC3 robust standard errors' : 'ANCOVA partial F-test',
    groupSizes: allConditions.map((condition) =>
      new Set([...aggregated.sampleMeta.values()].filter((row) => row.condition === condition).map((row) => row.subjectId)).size
    ),
    steps: aggregated.steps
  };
}

function analyseOutcomeLayer(aggregated, layer, { outcomeType = 'continuous', covariateColumns = [], targetTimepoint = '' } = {}) {
  const rows = [];
  for (const feature of aggregated.features) {
    const entries = featureOutcomeRows(aggregated, feature, outcomeType, targetTimepoint);
    if (entries.length < 4) continue;
    if (outcomeType === 'multiclass') {
      const groups = naturalOrder(entries.map((entry) => entry.outcome));
      if (groups.length < 2) continue;
      const nuisance = subjectCovariateDesign(entries, covariateColumns, { includeIntercept: true, includeBatch: true });
      const classLevels = groups.slice(1);
      const reducedDesign = nuisance.design;
      const fullDesign = entries.map((entry, i) => [
        1,
        ...classLevels.map((group) => entry.outcome === group ? 1 : 0),
        ...nuisance.design[i].slice(1)
      ]);
      const response = entries.map((entry) => entry.feature);
      const reduced = ordinaryLeastSquares(reducedDesign, response, 0);
      const full = ordinaryLeastSquares(fullDesign, response, 0);
      if (!reduced || !full) continue;
      const rssReduced = reduced.residuals.reduce((sum, value) => sum + value * value, 0);
      const rssFull = full.residuals.reduce((sum, value) => sum + value * value, 0);
      const df1 = classLevels.length;
      const df2 = response.length - fullDesign[0].length;
      if (!(df1 > 0) || !(df2 > 0)) continue;
      const numerator = Math.max(0, (rssReduced - rssFull) / df1);
      const denominator = rssFull / df2;
      const fStatistic = denominator > 0 ? numerator / denominator : Infinity;
      const pValue = Number.isFinite(fStatistic) || fStatistic === Infinity
        ? Math.max(0, Math.min(1, 1 - fCdf(fStatistic, df1, df2)))
        : null;
      const adjustedMeans = [full.beta[0], ...classLevels.map((_, index) => full.beta[0] + full.beta[index + 1])];
      rows.push({
        feature,
        effect: Math.max(...adjustedMeans) - Math.min(...adjustedMeans),
        effectScale: 'adjusted between-class difference',
        pValue,
        qValue: null,
        statistic: fStatistic,
        groupMeans: Object.fromEntries(groups.map((group, index) => [group, adjustedMeans[index]])),
        n: entries.length,
        model: 'ANCOVA partial F-test (class + batch + covariates)'
      });
      continue;
    }

    const cov = subjectCovariateDesign(entries, covariateColumns, { includeIntercept: true, includeBatch: true });
    const design = entries.map((entry, i) => [1, entry.feature, ...cov.design[i].slice(1)]);
    let fit = null;
    let effectScale = 'regression coefficient';
    let model = '';

    if (outcomeType === 'binary') {
      const levels = naturalOrder(entries.map((entry) => entry.outcome));
      if (levels.length !== 2) continue;
      const response = entries.map((entry) => entry.outcome === levels[1] ? 1 : 0);
      fit = fitGeneralizedLinear(design, response, 'binomial', 1);
      effectScale = 'log odds ratio';
      model = 'logistic regression (' + levels[1] + ' vs ' + levels[0] + ')';
    } else if (outcomeType === 'count') {
      const response = entries.map((entry) => finiteNumber(entry.outcome));
      if (response.some((value) => !Number.isFinite(value) || value < 0)) continue;
      fit = fitGeneralizedLinear(design, response, 'poisson', 1);
      effectScale = 'log rate ratio';
      model = 'Poisson regression';
    } else if (outcomeType === 'survival') {
      const responseDesign = entries.map((entry, i) => [entry.feature, ...cov.design[i].slice(1)]);
      const times = entries.map((entry) => entry.survivalTime);
      const events = entries.map((entry) => finiteNumber(entry.survivalEvent) > 0 ? 1 : 0);
      fit = coxRegression(responseDesign, times, events, 0);
      effectScale = 'log hazard ratio';
      model = 'Cox proportional hazards (Breslow ties; exploratory coefficient)';
      if (fit) fit.beta = [fit.beta[0]];
    } else {
      const response = entries.map((entry) => finiteNumber(entry.outcome));
      if (response.some((value) => !Number.isFinite(value))) continue;
      fit = ordinaryLeastSquares(design, response, 1);
      effectScale = 'outcome units per feature unit';
      model = 'linear regression';
    }

    if (!fit) continue;
    const coefficient = fit.beta?.[1] ?? fit.beta?.[0];
    // OLS uses t(df), while GLM/Cox use asymptotic normal Wald intervals.
    const critical = outcomeType === 'continuous' ? tCritical95(fit.df) : 1.96;
    const ciLow = Number.isFinite(fit.se) ? coefficient - critical * fit.se : null;
    const ciHigh = Number.isFinite(fit.se) ? coefficient + critical * fit.se : null;
    const multiplicative = ['binary','count','survival'].includes(outcomeType);
    rows.push({
      feature,
      effect: coefficient,
      effectScale,
      // Browser Cox cannot assess proportional hazards, censoring mechanism,
      // time-dependent exposure or informative censoring. Even after
      // verifying Breslow coefficients against R, its Wald p/q and CI are
      // intentionally unavailable until a method-specific R diagnostic.
      pValue: outcomeType === 'survival' ? null : fit.pValue,
      qValue: null,
      statistic: outcomeType === 'survival' ? null : fit.statistic,
      standardError: outcomeType === 'survival' ? null : fit.se,
      ciLow: outcomeType === 'survival' ? null : ciLow,
      ciHigh: outcomeType === 'survival' ? null : ciHigh,
      inferentialStatus: outcomeType === 'survival'
        ? 'exploratory_cox_reference_R_and_PH_diagnostics_required' : 'browser_screening',
      tieMethod: outcomeType === 'survival' ? 'breslow' : null,
      n: entries.length,
      model,
      exponentiatedEffect: multiplicative ? Math.exp(coefficient) : null,
      exponentiatedCiLow: multiplicative && outcomeType !== 'survival' && Number.isFinite(ciLow) ? Math.exp(ciLow) : null,
      exponentiatedCiHigh: multiplicative && outcomeType !== 'survival' && Number.isFinite(ciHigh) ? Math.exp(ciHigh) : null
    });
  }

  bhAdjust(rows);
  rows.sort((a,b) => {
    const aq = Number.isFinite(a.qValue) ? a.qValue : 1;
    const bq = Number.isFinite(b.qValue) ? b.qValue : 1;
    if (aq !== bq) return aq - bq;
    return Math.abs(b.effect) - Math.abs(a.effect);
  });
  const significant = rows.filter((row) => Number.isFinite(row.qValue) && row.qValue <= 0.10);
  const selected = (significant.length >= 10 ? significant : rows).slice(0, significant.length >= 10 ? 50 : 25);
  return {
    rows,
    selected,
    selectionRule: outcomeType === 'survival'
      ? 'Exploratory Cox effect ranking only: browser p/q-values and confidence intervals are withheld until R survival and proportional-hazards checks.'
      : significant.length >= 10
        ? 'q ≤ 0.10, capped at 50 features'
        : 'only ' + significant.length + ' FDR-qualified feature(s); top ' + selected.length + ' ranked features retained for exploratory biological mapping',
    effectScale: rows[0]?.effectScale || 'outcome association',
    contrast: 'association with ' + outcomeType + ' outcome',
    mode: 'outcome-' + outcomeType,
    inferenceMethod: outcomeType === 'survival'
      ? 'Exploratory Cox Breslow coefficient; inferential p/q/CI suppressed, independent R survival PH diagnostics required'
      : rows[0]?.model || 'outcome regression',
    inferenceStatus: outcomeType === 'survival'
      ? 'reference_R_survival_and_ph_diagnostics_required'
      : 'browser_screening',
    groupSizes: [],
    steps: aggregated.steps
  };
}

function endpointSubjectValueMap(aggregated, feature) {
  return new Map(subjectEndpointRows(aggregated, feature).map((entry) => [entry.subjectId, entry.value]));
}

function powerEigen(matrix, seed = 1, orthogonalTo = []) {
  const n = matrix.length;
  if (!n) return null;
  let vector = Array.from({ length: n }, (_, i) => Math.sin((i + 1) * (seed + 0.37)));
  const normalise = (v) => {
    for (const basis of orthogonalTo) {
      const dot = v.reduce((sum, value, i) => sum + value * basis[i], 0);
      for (let i = 0; i < v.length; i += 1) v[i] -= dot * basis[i];
    }
    const norm = Math.sqrt(v.reduce((sum, value) => sum + value * value, 0));
    return norm > 0 ? v.map((value) => value / norm) : null;
  };
  vector = normalise(vector);
  if (!vector) return null;
  for (let iter = 0; iter < 200; iter += 1) {
    let next = multiplyMatrixVector(matrix, vector);
    next = normalise(next);
    if (!next) return null;
    const delta = Math.sqrt(next.reduce((sum, value, i) => sum + (value - vector[i]) ** 2, 0));
    vector = next;
    if (delta < 1e-9) break;
  }
  const mv = multiplyMatrixVector(matrix, vector);
  const eigenvalue = vector.reduce((sum, value, i) => sum + value * mv[i], 0);
  return { vector, eigenvalue };
}

function layerQcPca(aggregated, maxFeatures = 100) {
  const samples = [...aggregated.sampleMeta.entries()];
  if (samples.length < 3) return { scores: [], explained: [], method: 'not enough samples' };
  const features = topVariableFeatures(aggregated, maxFeatures);
  const matrix = samples.map(() => []);
  for (const feature of features) {
    const source = aggregated.values.get(feature);
    const observed = samples.map(([sampleId]) => source.get(sampleId)).filter(Number.isFinite);
    if (observed.length < 2) continue;
    const center = mean(observed);
    const sd = Math.sqrt(variance(observed));
    if (!(sd > 0)) continue;
    const fill = median(observed);
    for (let i = 0; i < samples.length; i += 1) {
      const raw = source.get(samples[i][0]);
      matrix[i].push(((Number.isFinite(raw) ? raw : fill) - center) / sd);
    }
  }
  if (!matrix[0]?.length) return { scores: [], explained: [], method: 'no variable features' };

  const gram = Array.from({ length: samples.length }, () => Array(samples.length).fill(0));
  let totalVariance = 0;
  for (let i = 0; i < samples.length; i += 1) {
    for (let j = i; j < samples.length; j += 1) {
      let value = 0;
      for (let k = 0; k < matrix[i].length; k += 1) value += matrix[i][k] * matrix[j][k];
      gram[i][j] = value;
      gram[j][i] = value;
    }
    totalVariance += gram[i][i];
  }

  const components = [];
  const basis = [];
  for (let component = 0; component < Math.min(2, samples.length - 1); component += 1) {
    const eig = powerEigen(gram, 101 + component, basis);
    if (!eig || !(eig.eigenvalue > 1e-10)) continue;
    basis.push(eig.vector);
    components.push({
      eigenvalue: eig.eigenvalue,
      explainedFraction: totalVariance > 0 ? eig.eigenvalue / totalVariance : null,
      scores: eig.vector.map((value) => value * Math.sqrt(eig.eigenvalue))
    });
  }
  const scores = samples.map(([sampleId, meta], i) => ({
    sampleId,
    subjectId: meta.subjectId,
    condition: meta.condition,
    batch: meta.batch,
    pc1: components[0]?.scores[i] ?? 0,
    pc2: components[1]?.scores[i] ?? 0
  }));
  return {
    scores,
    explained: components.map((component) => component.explainedFraction),
    featuresUsed: matrix[0].length,
    method: 'QC PCA on processed features; feature-wise median imputation is used only for this visualization'
  };
}

function analyseExploratoryIntegration(aggregatedByLayer, loadedLayers, maxFeaturesPerLayer = 50, allowPartialBlocks = false) {
  const subjectSets = Object.fromEntries(loadedLayers.map((layer) => [
    layer,
    new Set([...aggregatedByLayer[layer].sampleMeta.values()].map((row) => row.subjectId))
  ]));
  const unionSubjects = [...new Set(loadedLayers.flatMap((layer) => [...subjectSets[layer]]))].sort();
  const completeSubjects = unionSubjects.filter((subject) => loadedLayers.every((layer) => subjectSets[layer].has(subject)));
  const subjects = allowPartialBlocks ? unionSubjects : completeSubjects;
  if (subjects.length < 3) {
    return {
      error: loadedLayers.length === 1
        ? 'Single-omics exploratory PCA requires at least three independent subjects.'
        : allowPartialBlocks
          ? 'Exploratory multi-omics integration requires at least three subjects represented in at least one loaded layer.'
          : 'Exploratory multi-omics integration requires at least three subjects shared across all loaded layers.',
      subjects: subjects.length,
      components: [],
      layers: {}
    };
  }

  const subjectCoverage = subjects.map((subject) => ({
    subjectId: subject,
    observedLayers: loadedLayers.filter((layer) => subjectSets[layer].has(subject)),
    missingLayers: loadedLayers.filter((layer) => !subjectSets[layer].has(subject))
  }));
  const blocks = {};
  const concatenated = subjects.map(() => []);
  const featureMap = [];
  let totalVariance = 0;
  let imputedCells = 0;

  for (const layer of loadedLayers) {
    const features = topVariableFeatures(aggregatedByLayer[layer], maxFeaturesPerLayer);
    const standardizedColumns = [];
    const keptFeatures = [];
    for (const feature of features) {
      const map = endpointSubjectValueMap(aggregatedByLayer[layer], feature);
      const rawValues = subjects.map((subject) => map.get(subject));
      const observed = rawValues.filter(Number.isFinite);
      if (observed.length < 2) continue;
      if (!allowPartialBlocks && rawValues.some((value) => !Number.isFinite(value))) continue;
      const center = mean(observed);
      const sd = Math.sqrt(variance(observed));
      if (!(sd > 0)) continue;
      keptFeatures.push(feature);
      standardizedColumns.push(rawValues.map((value) => {
        if (!Number.isFinite(value)) {
          imputedCells += 1;
          return 0;
        }
        return (value - center) / sd;
      }));
    }
    const blockScale = keptFeatures.length ? 1 / Math.sqrt(keptFeatures.length) : 1;
    for (let fIndex = 0; fIndex < keptFeatures.length; fIndex += 1) {
      const values = standardizedColumns[fIndex].map((value) => value * blockScale);
      featureMap.push({ layer, feature: keptFeatures[fIndex], values });
      for (let i = 0; i < subjects.length; i += 1) concatenated[i].push(values[i]);
      totalVariance += values.reduce((sum, value) => sum + value * value, 0);
    }
    blocks[layer] = {
      featureCount: keptFeatures.length,
      features: keptFeatures,
      subjectsObserved: subjects.filter((subject) => subjectSets[layer].has(subject)).length
    };
  }

  if (!concatenated[0]?.length) {
    return { error: 'No variable features were available for exploratory integration.', subjects: subjects.length, components: [], layers: blocks };
  }

  const gram = Array.from({ length: subjects.length }, () => Array(subjects.length).fill(0));
  for (let i = 0; i < subjects.length; i += 1) {
    for (let j = i; j < subjects.length; j += 1) {
      let sum = 0;
      for (let k = 0; k < concatenated[i].length; k += 1) sum += concatenated[i][k] * concatenated[j][k];
      gram[i][j] = sum;
      gram[j][i] = sum;
    }
  }

  const components = [];
  const basis = [];
  for (let component = 0; component < Math.min(2, subjects.length - 1); component += 1) {
    const eig = powerEigen(gram, component + 1, basis);
    if (!eig || !(eig.eigenvalue > 1e-10)) continue;
    basis.push(eig.vector);
    const scores = eig.vector.map((value) => value * Math.sqrt(eig.eigenvalue));
    const loadings = featureMap.map(({ layer, feature, values }) => {
      const loading = values.reduce((sum, value, i) => sum + value * eig.vector[i], 0) / Math.sqrt(eig.eigenvalue);
      return { layer, feature, loading };
    });
    components.push({
      component: component + 1,
      explainedFraction: totalVariance > 0 ? eig.eigenvalue / totalVariance : null,
      scores: subjects.map((subject, i) => ({ subjectId: subject, score: scores[i] })),
      topLoadings: loadings.slice().sort((a,b) => Math.abs(b.loading) - Math.abs(a.loading)).slice(0,30),
      layerTopLoadings: Object.fromEntries(loadedLayers.map((layer) => [
        layer,
        loadings.filter((item) => item.layer === layer)
          .sort((a,b) => Math.abs(b.loading) - Math.abs(a.loading))
          .slice(0,25)
      ])),
      layerContribution: Object.fromEntries(loadedLayers.map((layer) => {
        const sumSquares = loadings.filter((item) => item.layer === layer).reduce((sum, item) => sum + item.loading * item.loading, 0);
        return [layer, sumSquares];
      }))
    });
  }

  return {
    method: loadedLayers.length === 1
      ? 'single-omics exploratory PCA: top-variable features standardized within their layer (1/sqrt(p) scaling is common and does not affect within-layer PCA); descriptive only, no multi-omics integration'
      : allowPartialBlocks
        ? 'balanced multi-block PCA with partial-block support: features are standardized within observed values; missing feature/block entries are set to the standardized mean (0) only for latent exploration, never for differential tests'
        : 'balanced multi-block PCA: top-variable features are standardized within layer, each layer is scaled by 1/sqrt(p), then deterministic PCA is computed on shared subjects',
    missingDataPolicy: allowPartialBlocks
      ? 'partial blocks allowed; standardized-mean fill only inside unsupervised latent integration'
      : 'complete subjects across all loaded layers',
    subjects: subjects.length,
    completeSubjects: completeSubjects.length,
    subjectIds: subjects,
    subjectCoverage,
    imputedCells,
    analysisMode: loadedLayers.length === 1 ? 'single_omic' : 'multiomics',
    layers: blocks,
    components
  };
}

function analyseSupervisedMultiblock(aggregatedByLayer, layers, loadedLayers, metadata, protocol) {
  const supported = protocol.objective === 'groups'
    ? true
    : protocol.objective === 'outcome' && ['binary','continuous','count'].includes(protocol.outcomeType);
  if (!supported) return null;

  const targetBySubject = new Map();
  if (protocol.objective === 'groups') {
    const conditions = naturalOrder(metadata.map((row) => row.condition));
    if (conditions.length !== 2) return null;
    for (const row of metadata) targetBySubject.set(row.subjectId, row.condition === conditions[1] ? 1 : 0);
  } else if (protocol.outcomeType === 'binary') {
    const levels = naturalOrder(metadata.map((row) => row.outcome));
    if (levels.length !== 2) return null;
    for (const row of metadata) if (row.outcome) targetBySubject.set(row.subjectId, row.outcome === levels[1] ? 1 : 0);
  } else {
    for (const row of metadata) {
      const value = finiteNumber(row.outcome);
      if (Number.isFinite(value)) targetBySubject.set(row.subjectId, value);
    }
  }

  const subjectSets = Object.fromEntries(loadedLayers.map((layer) => [
    layer,
    new Set([...aggregatedByLayer[layer].sampleMeta.values()].map((row) => row.subjectId))
  ]));
  const allSubjects = [...targetBySubject.keys()]
    .filter((subject) => loadedLayers.some((layer) => subjectSets[layer].has(subject)))
    .sort();
  if (allSubjects.length < 6) return null;

  const yRaw = allSubjects.map((subject) => targetBySubject.get(subject));
  const yMean = mean(yRaw);
  const ySd = Math.sqrt(variance(yRaw));
  if (!(ySd > 0)) return null;
  const y = yRaw.map((value) => (value - yMean) / ySd);

  const columns = [];
  const matrix = allSubjects.map(() => []);
  for (const layer of loadedLayers) {
    const candidates = integrationCandidates(aggregatedByLayer[layer], layers[layer], 30);
    const blockColumns = [];
    for (const feature of candidates) {
      const targetTimepoint = protocol.objective === 'outcome' ? (protocol.outcomeTimepoint || '') : '';
      const entries = subjectEndpointRows(aggregatedByLayer[layer], feature, targetTimepoint);
      const map = new Map(entries.map((entry) => [entry.subjectId, entry.value]));
      const observed = allSubjects.map((subject) => map.get(subject)).filter(Number.isFinite);
      if (observed.length < 4) continue;
      const center = mean(observed);
      const sd = Math.sqrt(variance(observed));
      if (!(sd > 0)) continue;
      blockColumns.push({
        layer,
        feature,
        values: allSubjects.map((subject) => {
          const value = map.get(subject);
          return Number.isFinite(value) ? (value - center) / sd : 0;
        })
      });
    }
    const scale = blockColumns.length ? 1 / Math.sqrt(blockColumns.length) : 1;
    for (const column of blockColumns) {
      const values = column.values.map((value) => value * scale);
      columns.push({ layer: column.layer, feature: column.feature, values });
      for (let i = 0; i < allSubjects.length; i += 1) matrix[i].push(values[i]);
    }
  }
  if (!columns.length) return null;

  let weights = columns.map((column) =>
    column.values.reduce((sum, value, i) => sum + value * y[i], 0) / Math.max(1, allSubjects.length - 1)
  );
  const norm = Math.sqrt(weights.reduce((sum, value) => sum + value * value, 0));
  if (!(norm > 0)) return null;
  weights = weights.map((value) => value / norm);
  const scores = matrix.map((row) => row.reduce((sum, value, j) => sum + value * weights[j], 0));
  const scoreTargetCorrelation = pearson(scores, yRaw);

  const weighted = columns.map((column, index) => ({
    layer: column.layer,
    feature: column.feature,
    weight: weights[index],
    absoluteWeight: Math.abs(weights[index])
  })).sort((a,b) => b.absoluteWeight - a.absoluteWeight);

  return {
    method: 'balanced multiblock PLS1-style supervised component: standardized candidate features, 1/sqrt(p) block balancing, covariance weights with the declared target',
    target: protocol.objective === 'groups' ? 'condition' : protocol.outcomeType + ' outcome',
    subjects: allSubjects.length,
    scoreTargetCorrelation,
    scoreTargetR2: Number.isFinite(scoreTargetCorrelation) ? scoreTargetCorrelation * scoreTargetCorrelation : null,
    scores: allSubjects.map((subject, i) => ({
      subjectId: subject,
      score: scores[i],
      target: yRaw[i]
    })),
    topWeights: weighted.slice(0,40),
    layerContribution: Object.fromEntries(loadedLayers.map((layer) => [
      layer,
      weighted.filter((item) => item.layer === layer).reduce((sum, item) => sum + item.weight * item.weight, 0)
    ])),
    caveat: 'This is a descriptive supervised latent component, not a cross-validated predictive model and not the mixOmics DIABLO implementation.'
  };
}

function ridgeLinearFit(design, response, lambda) {
  const xtx = crossProductMatrix(design);
  const xty = crossProductVector(design, response);
  for (let j = 1; j < xtx.length; j += 1) xtx[j][j] += lambda;
  const beta = solveLinearSystem(xtx, xty, 1e-10);
  return beta ? { beta } : null;
}

function ridgeGlmFit(design, response, family, lambda) {
  const n = response.length;
  const p = design[0]?.length || 0;
  if (!n || !p) return null;
  let beta = Array(p).fill(0);
  if (family === 'poisson') beta[0] = Math.log(Math.max(mean(response), 1e-6));
  for (let iter = 0; iter < 50; iter += 1) {
    const eta = multiplyMatrixVector(design, beta);
    const weights = [];
    const z = [];
    for (let i = 0; i < n; i += 1) {
      if (family === 'binomial') {
        const e = Math.max(-30, Math.min(30, eta[i]));
        const mu = 1 / (1 + Math.exp(-e));
        const w = Math.max(mu * (1 - mu), 1e-6);
        weights.push(w);
        z.push(eta[i] + (response[i] - mu) / w);
      } else {
        const mu = Math.max(1e-8, Math.min(1e8, Math.exp(Math.max(-20, Math.min(20, eta[i])))));
        weights.push(Math.max(mu, 1e-6));
        z.push(eta[i] + (response[i] - mu) / Math.max(mu, 1e-6));
      }
    }
    const xtwx = crossProductMatrix(design, weights);
    const xtwz = crossProductVector(design, z, weights);
    for (let j = 1; j < p; j += 1) xtwx[j][j] += lambda;
    const next = solveLinearSystem(xtwx, xtwz, 1e-8);
    if (!next) return null;
    const delta = Math.max(...next.map((value, j) => Math.abs(value - beta[j])));
    beta = next;
    if (delta < 1e-7) break;
  }
  return { beta };
}

function aucScore(labels, scores) {
  const pairs = labels.map((label, i) => ({ label, score: scores[i] }))
    .filter((item) => Number.isFinite(item.score) && (item.label === 0 || item.label === 1))
    .sort((a,b) => a.score - b.score);
  const n1 = pairs.filter((item) => item.label === 1).length;
  const n0 = pairs.length - n1;
  if (!n1 || !n0) return NaN;
  let rankSum = 0;
  let i = 0;
  while (i < pairs.length) {
    let j = i + 1;
    while (j < pairs.length && pairs[j].score === pairs[i].score) j += 1;
    const avgRank = (i + 1 + j) / 2;
    for (let k = i; k < j; k += 1) if (pairs[k].label === 1) rankSum += avgRank;
    i = j;
  }
  return (rankSum - n1 * (n1 + 1) / 2) / (n1 * n0);
}

function deterministicFolds(subjects, targetBySubject, k, categorical = false) {
  const folds = Array.from({ length: Math.max(2, Math.min(k, subjects.length)) }, () => []);
  if (categorical) {
    const groups = new Map();
    for (const subject of subjects) {
      const key = String(targetBySubject.get(subject));
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(subject);
    }
    for (const group of groups.values()) {
      group.sort();
      group.forEach((subject, index) => folds[index % folds.length].push(subject));
    }
  } else {
    subjects.slice().sort().forEach((subject, index) => folds[index % folds.length].push(subject));
  }
  return folds.filter((fold) => fold.length);
}

function predictionTarget(metadata, protocol) {
  const targetBySubject = new Map();
  if (protocol.outcomeType === 'binary' || protocol.outcomeType === 'multiclass') {
    const levels = naturalOrder(metadata.map((row) => row.outcome));
    for (const row of metadata) if (row.outcome) targetBySubject.set(row.subjectId, row.outcome);
    return { targetBySubject, levels, categorical: true };
  }
  if (['continuous','count'].includes(protocol.outcomeType)) {
    for (const row of metadata) {
      const value = finiteNumber(row.outcome);
      if (Number.isFinite(value)) targetBySubject.set(row.subjectId, value);
    }
    return { targetBySubject, levels: [], categorical: false };
  }
  if (protocol.outcomeType === 'survival') {
    for (const row of metadata) {
      const time = finiteNumber(row.survivalTime);
      const event = finiteNumber(row.survivalEvent);
      if (Number.isFinite(time) && time > 0 && Number.isFinite(event)) {
        targetBySubject.set(row.subjectId, { time, event: event > 0 ? 1 : 0 });
      }
    }
    return { targetBySubject, levels: [], categorical: true, survival: true };
  }
  return { targetBySubject, levels: [], categorical: false };
}

function foldAdjustedFeatureEntries(aggregated, feature, trainSubjects, subjects, protocol) {
  const entries = subjectEndpointRows(aggregated, feature, protocol.outcomeTimepoint || '');
  const bySubject = new Map(entries.map((entry) => [entry.subjectId, entry]));
  const trainEntries = trainSubjects.map((subject) => bySubject.get(subject)).filter(Boolean);
  const requestedEntries = subjects.map((subject) => bySubject.get(subject)).filter(Boolean);
  if (!trainEntries.length) return [];

  const encoder = buildCovariateEncoder(
    trainEntries.map((entry) => entry.row),
    protocol.covariateColumns || [],
    { includeBatch: true }
  );

  let beta = null;
  if (encoder.columnNames.length > 1 && trainEntries.length > encoder.columnNames.length + 1) {
    const trainDesign = trainEntries.map((entry) => encoder.encode(entry.row));
    const response = trainEntries.map((entry) => entry.value);
    beta = solveLinearSystem(crossProductMatrix(trainDesign), crossProductVector(trainDesign, response));
  }

  const trainResiduals = trainEntries.map((entry) => {
    if (!beta) return entry.value;
    const fitted = encoder.encode(entry.row).reduce((sum, value, j) => sum + value * beta[j], 0);
    return entry.value - fitted;
  });
  const residualCenter = mean(trainResiduals);

  return requestedEntries.map((entry) => {
    if (!beta) return entry;
    const fitted = encoder.encode(entry.row).reduce((sum, value, j) => sum + value * beta[j], 0);
    return { ...entry, value: entry.value - fitted + residualCenter };
  });
}

function selectPredictionFeatures(aggregatedByLayer, loadedLayers, trainSubjects, target, protocol, maxPerLayer = 12) {
  const trainSet = new Set(trainSubjects);
  const selected = [];
  for (const layer of loadedLayers) {
    const scored = [];
    for (const feature of aggregatedByLayer[layer].features) {
      const entries = foldAdjustedFeatureEntries(
        aggregatedByLayer[layer],
        feature,
        trainSubjects,
        trainSubjects,
        protocol
      );
      const usable = entries.filter((entry) => trainSet.has(entry.subjectId) && target.targetBySubject.has(entry.subjectId));
      if (usable.length < Math.max(5, Math.ceil(trainSubjects.length * 0.50))) continue;
      const x = usable.map((entry) => entry.value);
      let score = 0;
      if (protocol.outcomeType === 'binary') {
        const levels = target.levels;
        const a = usable.filter((entry) => target.targetBySubject.get(entry.subjectId) === levels[0]).map((entry) => entry.value);
        const b = usable.filter((entry) => target.targetBySubject.get(entry.subjectId) === levels[1]).map((entry) => entry.value);
        if (a.length < 2 || b.length < 2) continue;
        const pooled = Math.sqrt(Math.max(1e-12, ((a.length - 1) * variance(a) + (b.length - 1) * variance(b)) / Math.max(1, a.length + b.length - 2)));
        score = Math.abs(mean(b) - mean(a)) / pooled;
      } else if (protocol.outcomeType === 'multiclass') {
        const groups = target.levels.map((level) =>
          usable.filter((entry) => target.targetBySubject.get(entry.subjectId) === level).map((entry) => entry.value)
        );
        if (groups.some((group) => group.length < 2)) continue;
        score = oneWayF(groups);
      } else if (protocol.outcomeType === 'survival') {
        const truth = usable.map((entry) => target.targetBySubject.get(entry.subjectId));
        const times = truth.map((value) => value.time);
        const events = truth.map((value) => value.event);
        if (events.reduce((sum, value) => sum + value, 0) < 3) continue;
        const center = mean(x);
        const sd = Math.sqrt(variance(x));
        if (!(sd > 0)) continue;
        const design = x.map((value) => [(value - center) / sd]);
        const fit = coxRegression(design, times, events, 0);
        score = fit && Number.isFinite(fit.statistic) ? Math.abs(fit.statistic) : NaN;
      } else {
        const y = usable.map((entry) => target.targetBySubject.get(entry.subjectId));
        score = Math.abs(pearson(x,y));
      }
      if (Number.isFinite(score)) scored.push({ layer, feature, score });
    }
    scored.sort((a,b) => b.score - a.score);
    selected.push(...scored.slice(0, maxPerLayer));
  }
  return selected;
}

function predictionMatrix(selected, aggregatedByLayer, subjects, trainSubjects, protocol) {
  const columns = [];
  for (const item of selected) {
    const requestedEntries = foldAdjustedFeatureEntries(
      aggregatedByLayer[item.layer],
      item.feature,
      trainSubjects,
      subjects,
      protocol
    );
    const trainEntries = foldAdjustedFeatureEntries(
      aggregatedByLayer[item.layer],
      item.feature,
      trainSubjects,
      trainSubjects,
      protocol
    );
    const map = new Map(requestedEntries.map((entry) => [entry.subjectId, entry.value]));
    const trainMap = new Map(trainEntries.map((entry) => [entry.subjectId, entry.value]));
    const trainValues = trainSubjects.map((subject) => trainMap.get(subject)).filter(Number.isFinite);
    if (trainValues.length < 2) continue;
    const center = mean(trainValues);
    const sd = Math.sqrt(variance(trainValues));
    if (!(sd > 0)) continue;
    columns.push({
      ...item,
      center,
      sd,
      values: subjects.map((subject) => {
        const value = map.get(subject);
        return Number.isFinite(value) ? (value - center) / sd : 0;
      })
    });
  }
  const design = subjects.map((_, i) => [1, ...columns.map((column) => column.values[i])]);
  return { design, columns };
}

function predictFromBeta(design, beta, family) {
  const eta = multiplyMatrixVector(design, beta);
  if (family === 'binomial') return eta.map((value) => 1 / (1 + Math.exp(-Math.max(-30, Math.min(30, value)))));
  if (family === 'poisson') return eta.map((value) => Math.exp(Math.max(-20, Math.min(20, value))));
  return eta;
}

function validationLoss(protocol, truth, predictions, levels = []) {
  if (protocol.outcomeType === 'binary') {
    const y = truth.map((value) => value === levels[1] ? 1 : 0);
    const eps = 1e-8;
    return -mean(y.map((value, i) => value * Math.log(Math.max(eps, predictions[i])) + (1 - value) * Math.log(Math.max(eps, 1 - predictions[i]))));
  }
  if (protocol.outcomeType === 'multiclass') {
    return 1 - mean(truth.map((value, i) => predictions[i] === value ? 1 : 0));
  }
  if (protocol.outcomeType === 'survival') {
    const times = truth.map((value) => value.time);
    const events = truth.map((value) => value.event);
    const cIndex = harrellCIndex(times, events, predictions);
    return Number.isFinite(cIndex) ? 1 - cIndex : NaN;
  }
  return Math.sqrt(mean(truth.map((value, i) => (Number(value) - predictions[i]) ** 2)));
}

function fitPredictiveModel(design, truth, protocol, levels, lambda) {
  if (protocol.outcomeType === 'binary') {
    const response = truth.map((value) => value === levels[1] ? 1 : 0);
    return ridgeGlmFit(design, response, 'binomial', lambda);
  }
  if (protocol.outcomeType === 'count') return ridgeGlmFit(design, truth.map(Number), 'poisson', lambda);
  if (protocol.outcomeType === 'continuous') return ridgeLinearFit(design, truth.map(Number), lambda);
  if (protocol.outcomeType === 'survival') {
    const noIntercept = design.map((row) => row.slice(1));
    const times = truth.map((value) => value.time);
    const events = truth.map((value) => value.event);
    const fit = coxRidgeRegression(noIntercept, times, events, lambda);
    return fit ? { ...fit, survival: true } : null;
  }
  if (protocol.outcomeType === 'multiclass') {
    const models = [];
    for (const level of levels) {
      const response = truth.map((value) => value === level ? 1 : 0);
      const fit = ridgeGlmFit(design, response, 'binomial', lambda);
      if (!fit) return null;
      models.push({ level, beta: fit.beta });
    }
    return { models };
  }
  return null;
}

function predictModel(model, design, protocol) {
  if (!model) return [];
  if (protocol.outcomeType === 'multiclass') {
    const probabilities = model.models.map((entry) => ({
      level: entry.level,
      values: predictFromBeta(design, entry.beta, 'binomial')
    }));
    return design.map((_, i) => probabilities.slice().sort((a,b) => b.values[i] - a.values[i])[0].level);
  }
  if (protocol.outcomeType === 'survival') {
    return multiplyMatrixVector(design.map((row) => row.slice(1)), model.beta);
  }
  const family = protocol.outcomeType === 'binary' ? 'binomial' : protocol.outcomeType === 'count' ? 'poisson' : 'gaussian';
  return predictFromBeta(design, model.beta, family);
}

/**
 * Frozen binary prediction contract for a source-defined train/test partition.
 * Uses the SAME fold-feature ranking idea and ridgeGlmFit binomial solver as
 * analysePredictiveOutcome. It is intentionally *not* a fitted cross-omics
 * model: caller provides a single log-scale numerical feature matrix.
 *
 * Data: {features: string[], rows: {id: string, outcome?: string,
 * values: (number|null)[]}[]}. New/test subjects must never influence
 * feature selection, centering, scaling, lambda tuning or fitted coefficients.
 * No result here is automatically certified for clinical decisions.
 */
export function fitFrozenBinaryPredictor({
  features, rows, positiveLabel, negativeLabel, maxFeatures = 12
}) {
  if(!Array.isArray(features)||features.length<2||
     features.some(f=>typeof f!=='string'||!f.trim())||
     new Set(features).size!==features.length)
    throw new Error('Frozen predictor requires unique, nonempty feature identifiers.');
  if(!Array.isArray(rows)||rows.length<20)
    throw new Error('Frozen predictor requires at least 20 training subjects.');
  const ids=rows.map(r=>r.id);
  if(ids.some(x=>typeof x!=='string'||!x.trim())||new Set(ids).size!==ids.length)
    throw new Error('Frozen predictor training subject IDs must be unique and nonempty.');
  if(!positiveLabel||!negativeLabel||positiveLabel===negativeLabel)
    throw new Error('A fixed positive and negative outcome label is required.');
  if(rows.some(r=>r.outcome!==positiveLabel&&r.outcome!==negativeLabel))
    throw new Error('Training rows contain unexpected outcome classes.');
  if(rows.some(r=>!Array.isArray(r.values)||r.values.length!==features.length||
       r.values.some(v=>v!==null&&!Number.isFinite(v))))
    throw new Error('Training matrix must contain numeric values or explicit null, with exact feature order.');
  const classes=[negativeLabel,positiveLabel];
  const counts=classes.map(x=>rows.filter(r=>r.outcome===x).length);
  if(counts.some(n=>n<8))throw new Error('At least eight training subjects per class are required.');
  if(!Number.isInteger(maxFeatures)||maxFeatures<1||maxFeatures>48)
    throw new Error('Invalid frozen-predictor maximum feature count.');
  const byId=new Map(rows.map(r=>[r.id,r]));
  const targetBySubject=new Map(rows.map(r=>[r.id,r.outcome]));
  const orderedIds=ids.slice().sort();
  const folds=deterministicFolds(orderedIds,targetBySubject,Math.min(4,rows.length),true);
  const lambdas=[0.01,0.1,1,10];
  const label=(id)=>byId.get(id).outcome===positiveLabel?1:0;

  // Recompute selection and every preprocessing statistic on each
  // inner-training subset; never on the full cohort before tuning.
  const prepare=(trainIds)=>{
    const jScore=[];
    for(let j=0;j<features.length;j++){
      const observed=trainIds.map(id=>({value:byId.get(id).values[j],y:label(id)}))
        .filter(x=>Number.isFinite(x.value));
      if(observed.length<Math.max(5,Math.ceil(trainIds.length*.5)))continue;
      const a=observed.filter(x=>x.y===0).map(x=>x.value);
      const b=observed.filter(x=>x.y===1).map(x=>x.value);
      if(a.length<2||b.length<2)continue;
      const sd=Math.sqrt(Math.max(1e-12,
        ((a.length-1)*variance(a)+(b.length-1)*variance(b))/
        Math.max(1,a.length+b.length-2)));
      const score=Math.abs(mean(b)-mean(a))/sd;
      if(Number.isFinite(score))jScore.push({j,score,feature:features[j]});
    }
    jScore.sort((a,b)=>b.score-a.score||a.feature.localeCompare(b.feature));
    const chosen=[];
    for(const item of jScore.slice(0,maxFeatures)){
      const values=trainIds.map(id=>byId.get(id).values[item.j])
        .filter(Number.isFinite);
      const center=mean(values),scale=Math.sqrt(variance(values));
      if(!Number.isFinite(center)||!Number.isFinite(scale)||scale<1e-10)continue;
      chosen.push({index:item.j,feature:item.feature,center,scale});
    }
    if(!chosen.length)return null;
    const design=(sampleIds)=>sampleIds.map(id=>[
      1,...chosen.map(f=>{
        const value=byId.get(id).values[f.index];
        return Number.isFinite(value)?(value-f.center)/f.scale:0;
      })
    ]);
    return {chosen,design};
  };
  let bestLambda=null,bestLoss=Infinity;
  for(const lambda of lambdas){
    const losses=[];
    for(const validationIds of folds){
      const testSet=new Set(validationIds);
      const trainingIds=orderedIds.filter(id=>!testSet.has(id));
      if(!trainingIds.some(id=>label(id)===0)||!trainingIds.some(id=>label(id)===1))
        throw new Error('Internal tuning fold has only one biological outcome class.');
      const state=prepare(trainingIds);
      if(!state)throw new Error('No stable predictors remain within training-only tuning fold.');
      const truth=trainingIds.map(label);
      const fit=ridgeGlmFit(state.design(trainingIds),truth,'binomial',lambda);
      if(!fit||fit.beta.some(x=>!Number.isFinite(x)))
        throw new Error('Training-only ridge fit is nonestimable.');
      const predictions=predictFromBeta(state.design(validationIds),fit.beta,'binomial');
      if(predictions.some(x=>!Number.isFinite(x)))throw new Error('Invalid validation-fold probabilities.');
      losses.push(-mean(validationIds.map((id,i)=>{
        const y=label(id),p=Math.max(1e-8,Math.min(1-1e-8,predictions[i]));
        return y*Math.log(p)+(1-y)*Math.log(1-p);
      })));
    }
    const loss=mean(losses);
    if(!Number.isFinite(loss))throw new Error('Internal validation loss is nonfinite.');
    if(loss<bestLoss){bestLoss=loss;bestLambda=lambda;}
  }
  const state=prepare(orderedIds);
  if(!state||bestLambda===null)throw new Error('Cannot freeze selected features.');
  const fit=ridgeGlmFit(state.design(orderedIds),orderedIds.map(label),'binomial',bestLambda);
  if(!fit||fit.beta.length!==state.chosen.length+1||
     fit.beta.some(x=>!Number.isFinite(x)))
    throw new Error('Final PMx binomial ridge fit is nonestimable.');
  return {
    schema:'pmx_frozen_binary_ridge_v1',
    engine:'PMx deterministic-core ridgeGlmFit',
    featureRanking:'training-only absolute pooled-standardized mean difference',
    penaltyTuning:'4-fold stratified training-only cross-validation, log loss',
    positiveLabel,negativeLabel,
    fullTrainingSampleCount:rows.length,
    trainingClassCounts:{positive:counts[1],negative:counts[0]},
    trainingSubjectIds:orderedIds,
    requiredFeatures:state.chosen.map(({feature,center,scale})=>({feature,center,scale})),
    coefficients:fit.beta.slice(),
    lambda:bestLambda,innerLogLoss:bestLoss,
    inferentialScope:'Exploratory research only; training and test preprocessing provenance must be audited.'
  };
}

export function scoreFrozenBinaryPredictor(model,{features,rows}) {
  if(model?.schema!=='pmx_frozen_binary_ridge_v1'||!Array.isArray(model.requiredFeatures)||
     !Array.isArray(model.coefficients)||model.coefficients.length!==model.requiredFeatures.length+1)
    throw new Error('Incompatible frozen PMx binary model format.');
  if(!Array.isArray(features)||new Set(features).size!==features.length)
    throw new Error('Invalid evaluation feature schema.');
  if(!Array.isArray(rows)||!rows.length)
    throw new Error('No independent subjects to score.');
  const byFeature=new Map(features.map((f,i)=>[f,i]));
  const indices=model.requiredFeatures.map(f=>{
    if(!byFeature.has(f.feature))throw new Error('Missing frozen predictor feature: '+f.feature);
    return byFeature.get(f.feature);
  });
  const trainIds=new Set(model.trainingSubjectIds||[]);
  const scoredIds=new Set();
  const matrix=rows.map(row=>{
    if(typeof row.id!=='string'||!row.id.trim()||scoredIds.has(row.id))
      throw new Error('External subject identifiers must be nonempty and unique.');
    if(trainIds.has(row.id))throw new Error('Training subject leaked into the heldout evaluation.');
    scoredIds.add(row.id);
    if(!Array.isArray(row.values)||row.values.length!==features.length)
      throw new Error('External patient feature panel does not match its declared schema.');
    return [1,...model.requiredFeatures.map((f,k)=>{
      const v=row.values[indices[k]];
      if(v!==null&&!Number.isFinite(v))
        throw new Error('Invalid external measurement for '+f.feature);
      return v===null?0:(v-f.center)/f.scale;
    })];
  });
  const p=predictFromBeta(matrix,model.coefficients,'binomial');
  if(p.some(v=>!Number.isFinite(v)||v<0||v>1))
    throw new Error('Frozen PMx predictor produced an invalid probability.');
  return rows.map((r,i)=>({subject_id:r.id,prediction:p[i]}));
}

function analysePredictiveOutcome(aggregatedByLayer, loadedLayers, metadata, protocol) {
  if (protocol.objective !== 'outcome' || !['binary','continuous','count','multiclass','survival'].includes(protocol.outcomeType)) {
    return protocol.objective === 'outcome'
      ? { status: 'not_available', reason: 'Cross-validated prediction is unavailable for this outcome type.' }
      : null;
  }
  const target = predictionTarget(metadata, protocol);
  const subjects = [...target.targetBySubject.keys()].filter((subject) =>
    loadedLayers.some((layer) => [...aggregatedByLayer[layer].sampleMeta.values()].some((row) => row.subjectId === subject))
  ).sort();
  if (subjects.length < 12) return { status: 'not_available', reason: 'At least 12 subjects with outcome data are required for nested cross-validation.' };
  // Stratified five-fold validation needs each outcome class represented in
  // each fold. Otherwise AUC and inner-loop model selection are unstable.
  if (['binary','multiclass'].includes(protocol.outcomeType)) {
    const counts = target.levels.map((level) => subjects.filter((id) => target.targetBySubject.get(id) === level).length);
    if (counts.length < 2 || counts.some((n) => n < 5)) {
      return { status: 'not_available', reason: 'At least 5 independent subjects per outcome class are required for stratified five-fold nested cross-validation.' };
    }
  }

  const categorical = ['binary','multiclass','survival'].includes(protocol.outcomeType);
  const foldTarget = protocol.outcomeType === 'survival'
    ? new Map(subjects.map((subject) => [subject, target.targetBySubject.get(subject)?.event ?? 0]))
    : target.targetBySubject;
  if (protocol.outcomeType === 'survival') {
    const events = subjects.reduce((sum, subject) => sum + (target.targetBySubject.get(subject)?.event || 0), 0);
    if (events < 6) return { status: 'not_available', reason: 'At least 6 observed events are required for cross-validated survival prediction.' };
  }
  const outerFolds = deterministicFolds(subjects, foldTarget, Math.min(5, subjects.length), categorical);
  const lambdas = [0.01,0.1,1,10];
  const predictions = [];
  const foldSummaries = [];

  for (let foldIndex = 0; foldIndex < outerFolds.length; foldIndex += 1) {
    const testSubjects = outerFolds[foldIndex];
    const testSet = new Set(testSubjects);
    const trainSubjects = subjects.filter((subject) => !testSet.has(subject));
    const selected = selectPredictionFeatures(aggregatedByLayer, loadedLayers, trainSubjects, target, protocol, 12);
    if (!selected.length) continue;

    const innerFoldTarget = protocol.outcomeType === 'survival'
      ? new Map(trainSubjects.map((subject) => [subject, target.targetBySubject.get(subject)?.event ?? 0]))
      : target.targetBySubject;
    const innerFolds = deterministicFolds(trainSubjects, innerFoldTarget, Math.min(4, trainSubjects.length), categorical);
    // Preparing one inner-training dataset for each split (instead of once
    // for each penalty) is essential when thousands of features are uploaded.
    // All four tuning candidates use the *same*, properly isolated folds.
    const innerPrepared = [];
    for (const innerTest of innerFolds) {
      const innerTestSet = new Set(innerTest);
      const innerTrain = trainSubjects.filter((subject) => !innerTestSet.has(subject));
      if (innerTrain.length < 5 || !innerTest.length) continue;
      const innerSelected = selectPredictionFeatures(
        aggregatedByLayer, loadedLayers, innerTrain, target, protocol, 12
      );
      if (!innerSelected.length) continue;
      const trainMatrix = predictionMatrix(innerSelected, aggregatedByLayer, innerTrain, innerTrain, protocol);
      const testMatrix = predictionMatrix(innerSelected, aggregatedByLayer, innerTest, innerTrain, protocol);
      if (!trainMatrix.columns.length || trainMatrix.columns.length !== testMatrix.columns.length) continue;
      innerPrepared.push({
        training: trainMatrix.design,
        testing: testMatrix.design,
        trainTruth: innerTrain.map((subject) => target.targetBySubject.get(subject)),
        testTruth: innerTest.map((subject) => target.targetBySubject.get(subject))
      });
    }
    if (innerPrepared.length < 2) continue;
    let bestLambda = lambdas[0];
    let bestLoss = Infinity;
    for (const lambda of lambdas) {
      const losses = [];
      for (const split of innerPrepared) {
        const model = fitPredictiveModel(split.training, split.trainTruth, protocol, target.levels, lambda);
        if (!model) continue;
        const pred = predictModel(model, split.testing, protocol);
        losses.push(validationLoss(protocol, split.testTruth, pred, target.levels));
      }
      // Do not pick a penalty based on just one successful inner fold.
      if (losses.length !== innerPrepared.length) continue;
      const loss = mean(losses);
      if (Number.isFinite(loss) && loss < bestLoss) {
        bestLoss = loss;
        bestLambda = lambda;
      }
    }
    if (!Number.isFinite(bestLoss)) continue;

    const trainMatrix = predictionMatrix(selected, aggregatedByLayer, trainSubjects, trainSubjects, protocol);
    const testMatrix = predictionMatrix(selected, aggregatedByLayer, testSubjects, trainSubjects, protocol);
    if (!trainMatrix.columns.length || trainMatrix.columns.length !== testMatrix.columns.length) continue;
    const trainTruth = trainSubjects.map((subject) => target.targetBySubject.get(subject));
    const testTruth = testSubjects.map((subject) => target.targetBySubject.get(subject));
    const model = fitPredictiveModel(trainMatrix.design, trainTruth, protocol, target.levels, bestLambda);
    if (!model) continue;
    const pred = predictModel(model, testMatrix.design, protocol);
    testSubjects.forEach((subject, i) => predictions.push({ subjectId: subject, truth: testTruth[i], prediction: pred[i], fold: foldIndex + 1 }));
    foldSummaries.push({
      fold: foldIndex + 1,
      trainingSubjects: trainSubjects.length,
      testSubjects: testSubjects.length,
      selectedFeatures: selected.length,
      innerFeatureSelection: 'refitted on each inner-training fold',
      lambda: bestLambda,
      innerLoss: bestLoss
    });
  }

  // Fail closed: evaluating only the folds that converged can bias CV metrics.
  // Each eligible independent subject must have exactly one held-out prediction.
  const uniquePredictedSubjects = new Set(predictions.map((item) => item.subjectId));
  if (foldSummaries.length !== outerFolds.length ||
      predictions.length !== subjects.length ||
      uniquePredictedSubjects.size !== subjects.length) {
    return {
      status: 'not_available',
      reason: 'Incomplete outer cross-validation: every independent subject must have exactly one held-out prediction. Metrics are withheld when a fold fails.',
      attemptedFolds: outerFolds.length,
      completedFolds: foldSummaries.length,
      evaluatedSubjects: predictions.length,
      totalSubjects: subjects.length
    };
  }

  let metrics;
  if (protocol.outcomeType === 'binary') {
    const labels = predictions.map((item) => item.truth === target.levels[1] ? 1 : 0);
    const scores = predictions.map((item) => item.prediction);
    metrics = {
      auc: aucScore(labels, scores),
      accuracy: mean(labels.map((value, i) => (scores[i] >= 0.5 ? 1 : 0) === value ? 1 : 0)),
      // Proper quadratic probability score, evaluated ONLY on each held-out
      // subject's outer-fold probability. AUC ranks subjects; it does not
      // measure probability calibration. This is internal CV, not external
      // clinical validation or a calibration-slope estimate.
      brier: mean(labels.map((value, i) => (scores[i] - value) ** 2)),
      logLoss: validationLoss(protocol, predictions.map((item) => item.truth), scores, target.levels)
    };
  } else if (protocol.outcomeType === 'multiclass') {
    metrics = {
      accuracy: mean(predictions.map((item) => item.truth === item.prediction ? 1 : 0))
    };
  } else if (protocol.outcomeType === 'survival') {
    const truth = predictions.map((item) => item.truth);
    const risks = predictions.map((item) => Number(item.prediction));
    metrics = {
      cIndex: harrellCIndex(truth.map((item) => item.time), truth.map((item) => item.event), risks),
      events: truth.reduce((sum, item) => sum + item.event, 0)
    };
  } else {
    const truth = predictions.map((item) => Number(item.truth));
    const pred = predictions.map((item) => Number(item.prediction));
    const rmse = Math.sqrt(mean(truth.map((value, i) => (value - pred[i]) ** 2)));
    const baseline = mean(truth);
    const ssRes = truth.reduce((sum, value, i) => sum + (value - pred[i]) ** 2, 0);
    const ssTot = truth.reduce((sum, value) => sum + (value - baseline) ** 2, 0);
    metrics = { rmse, r2: ssTot > 0 ? 1 - ssRes / ssTot : null };
  }

  // Degenerate test sets (e.g. no comparable survival pairs) do not
  // establish predictive performance, even if the fold fits completed.
  const primaryMetric = protocol.outcomeType === 'binary' ? metrics.auc
    : protocol.outcomeType === 'multiclass' ? metrics.accuracy
      : protocol.outcomeType === 'survival' ? metrics.cIndex : metrics.rmse;
  if (!Number.isFinite(primaryMetric)) {
    return {
      status: 'not_available',
      reason: 'Cross-validated performance is undefined for these held-out outcomes; no metric is reported.',
      attemptedFolds: outerFolds.length,
      completedFolds: foldSummaries.length,
      evaluatedSubjects: predictions.length,
      totalSubjects: subjects.length
    };
  }

  return {
    status: 'ok',
    method: protocol.outcomeType === 'survival'
      ? 'nested cross-validation with training-only nuisance adjustment, univariate Cox feature screening and ridge-penalized Cox prediction'
      : 'nested cross-validation with training-only nuisance adjustment, feature selection and ridge regularisation',
    nuisanceAdjustment: {
      policy: 'fold-local',
      batch: true,
      covariates: protocol.covariateColumns || [],
      note: 'Nuisance adjustment, feature selection, missing-value fill, centering and scaling use each training split. The separate preprocessingPolicy records whether upstream operations are fold-isolated.'
    },
    outcomeType: protocol.outcomeType,
    subjects: subjects.length,
    outerFolds: foldSummaries.length,
    metrics,
    predictions,
    foldSummaries,
    caveat: 'Nested subject-level CV is internal screening, not external validation; performance is reported only for complete, one-prediction-per-subject out-of-fold coverage. Predictive feature selection is refitted inside each inner training fold. Upstream filtering is excluded from the supported as-supplied / integer-count predictive lane; unsupported scales may still undergo whole-cohort processing. External transportability and confounding remain unvalidated.'
  };
}

function explorationLayerResult(aggregated, exploration, layer) {
  const loadings = exploration.components?.[0]?.layerTopLoadings?.[layer]
    || exploration.components?.[0]?.topLoadings?.filter((item) => item.layer === layer)
    || [];
  const byFeature = new Map(loadings.map((item) => [item.feature, item.loading]));
  const rows = topVariableFeatures(aggregated, 100).map((feature) => ({
    feature,
    effect: byFeature.get(feature) ?? 0,
    effectScale: 'PC1 loading',
    pValue: null,
    qValue: null,
    variance: variance([...aggregated.values.get(feature).values()])
  })).sort((a,b) => Math.abs(b.effect) - Math.abs(a.effect) || b.variance - a.variance);
  const selected = rows.filter((row) => row.effect !== 0).slice(0,25);
  return {
    rows,
    selected,
    selectionRule: exploration.analysisMode === 'single_omic'
      ? 'top absolute loadings on single-omics PC1, capped at 25 features'
      : 'top absolute loadings on balanced multi-block PC1, capped at 25 features',
    effectScale: 'PC1 loading',
    contrast: exploration.analysisMode === 'single_omic' ? 'unsupervised within-layer variation' : 'unsupervised shared structure',
    mode: exploration.analysisMode === 'single_omic' ? 'exploratory-single-omic-pca' : 'exploratory-multiblock-pca',
    inferenceMethod: exploration.method,
    groupSizes: [],
    steps: aggregated.steps
  };
}

function subjectFeatureValues(aggregated, feature, options) {
  const source = aggregated.values.get(feature);
  const perSubject = new Map();
  for (const [sampleId, row] of aggregated.sampleMeta.entries()) {
    const value = source.get(sampleId);
    if (!Number.isFinite(value)) continue;
    if (!perSubject.has(row.subjectId)) perSubject.set(row.subjectId, []);
    perSubject.get(row.subjectId).push({ value, row });
  }

  const conditions = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.condition));
  const timepoints = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.timepoint));
  if (conditions.length !== 2) return { error: 'Current deterministic inference requires exactly two conditions.', conditions, timepoints };

  if (options.paired && !options.longitudinal) {
    const differences = [];
    for (const entries of perSubject.values()) {
      const a = entries.filter((x) => x.row.condition === conditions[0]).map((x) => x.value);
      const b = entries.filter((x) => x.row.condition === conditions[1]).map((x) => x.value);
      if (!a.length || !b.length) continue;
      differences.push(mean(b)-mean(a));
    }
    return {
      groups: [differences.map(() => 0), differences],
      differences,
      conditions,
      timepoints,
      contrast: `paired ${conditions[1]} − ${conditions[0]}`,
      mode: 'paired-permutation'
    };
  }

  if (options.longitudinal && timepoints.length >= 2) {
    const first = timepoints[0];
    const last = timepoints[timepoints.length - 1];
    const groups = [[],[]];

    if (timepoints.length > 2) {
      for (const entries of perSubject.values()) {
        const byTime = new Map();
        for (const entry of entries) {
          const time = numericTime(entry.row.timepoint);
          if (!Number.isFinite(time)) continue;
          if (!byTime.has(time)) byTime.set(time, []);
          byTime.get(time).push(entry.value);
        }
        const points = [...byTime.entries()].map(([time, values]) => ({ time, value: mean(values) }));
        const slope = linearSlope(points);
        if (!Number.isFinite(slope)) continue;
        const condition = entries[0].row.condition;
        const index = conditions.indexOf(condition);
        if (index >= 0) groups[index].push(slope);
      }
      return {
        groups,
        conditions,
        timepoints,
        contrast: `individual time slope: ${conditions[1]} vs ${conditions[0]}`,
        mode: 'longitudinal-slope'
      };
    }

    for (const entries of perSubject.values()) {
      const baseline = entries.filter((x) => x.row.timepoint === first).map((x) => x.value);
      const endpoint = entries.filter((x) => x.row.timepoint === last).map((x) => x.value);
      if (!baseline.length || !endpoint.length) continue;
      const condition = entries.find((x) => x.row.timepoint === last)?.row.condition || entries[0].row.condition;
      const index = conditions.indexOf(condition);
      if (index >= 0) groups[index].push(mean(endpoint) - mean(baseline));
    }
    return { groups, conditions, timepoints, contrast: `Δ${last}−${first}: ${conditions[1]} vs ${conditions[0]}`, mode: 'difference-in-differences' };
  }

  const targetTime = timepoints.length ? timepoints[timepoints.length - 1] : '';
  const groups = [[],[]];
  for (const entries of perSubject.values()) {
    const eligible = targetTime ? entries.filter((x) => x.row.timepoint === targetTime) : entries;
    if (!eligible.length) continue;
    const condition = eligible[0].row.condition;
    const index = conditions.indexOf(condition);
    if (index >= 0) groups[index].push(mean(eligible.map((x) => x.value)));
  }
  return { groups, conditions, timepoints, contrast: `${conditions[1]} vs ${conditions[0]}${targetTime ? ` at ${targetTime}` : ''}`, mode: 'two-group' };
}

function analyseLayer(aggregated, layer, options) {
  const rows = [];
  let contrast = '';
  let mode = '';
  let groupSizes = [0,0];
  const allConditions = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.condition));
  const multiGroup = allConditions.length > 2;
  const usePermutation = aggregated.features.length <= 500;
  const inferenceMethod = usePermutation
    ? 'deterministic permutation'
    : (multiGroup ? 'one-way ANOVA' : options.paired ? 'paired t-test' : 'Welch t-test');

  for (const feature of aggregated.features) {
    if (multiGroup) {
      if (options.longitudinal || options.paired) {
        return { error: 'Multi-group longitudinal/paired inference is not implemented.', rows: [], contrast: '', mode: '', groupSizes: [], selected: [], steps: aggregated.steps };
      }
      const data = independentMultiGroupValues(aggregated, feature);
      if (data.groups.some((group) => group.length < 2)) continue;
      const means = data.groups.map((group) => mean(group));
      const minMean = Math.min(...means);
      const maxMean = Math.max(...means);
      const effect = maxMean-minMean;
      contrast = `one-way condition effect across ${data.conditions.length} groups${data.targetTime ? ` at ${data.targetTime}` : ''}`;
      mode = 'multi-group-permutation-anova';
      groupSizes = data.groups.map((group) => group.length);
      rows.push({
        feature,
        effect,
        foldRatio: aggregated.scale === 'log2' ? Math.pow(2,effect) : null,
        effectScale: aggregated.scale,
        pValue: usePermutation
          ? multiGroupPermutationPValue(data.groups, `${layer}|${feature}|${contrast}`)
          : oneWayAnovaPValue(data.groups),
        qValue: null,
        groupMeans: Object.fromEntries(data.conditions.map((condition,index) => [condition, means[index]])),
        groupSizes: Object.fromEntries(data.conditions.map((condition,index) => [condition, data.groups[index].length])),
        nReference: null,
        nComparison: null
      });
      continue;
    }

    const data = subjectFeatureValues(aggregated, feature, options);
    if (data.error) return { error: data.error, rows: [], contrast: '', mode: '', groupSizes: [0,0], selected: [], steps: aggregated.steps };
    const [a,b] = data.groups;
    contrast = data.contrast;
    mode = data.mode;
    groupSizes = [a.length,b.length];
    if (!a.length || !b.length) continue;
    const effect = data.mode === 'paired-permutation' ? mean(data.differences) : mean(b) - mean(a);
    const pValue = data.mode === 'paired-permutation'
      ? (usePermutation
          ? pairedPermutationPValue(data.differences, `${layer}|${feature}|${contrast}`)
          : pairedTPValue(data.differences))
      : (usePermutation
          ? permutationPValue(a,b, `${layer}|${feature}|${contrast}`)
          : welchPValue(a,b));
    rows.push({
      feature,
      effect,
      foldRatio: aggregated.scale === 'log2' ? Math.pow(2,effect) : null,
      effectScale: aggregated.scale,
      pValue,
      qValue: null,
      nReference: a.length,
      nComparison: b.length
    });
  }
  bhAdjust(rows);
  rows.sort((a,b) => Math.abs(b.effect) - Math.abs(a.effect));
  const significant = rows.filter((row) =>
    Number.isFinite(row.qValue) &&
    row.qValue <= 0.10 &&
    (aggregated.scale !== 'log2' || Math.abs(row.effect) >= Math.log2(1.2))
  );
  const useConfirmedSet = significant.length >= 10;
  const selected = useConfirmedSet
    ? significant.slice(0, 50)
    : rows.slice().sort((a,b) => {
        const aq = Number.isFinite(a.qValue) ? a.qValue : 1;
        const bq = Number.isFinite(b.qValue) ? b.qValue : 1;
        if (aq !== bq) return aq - bq;
        return Math.abs(b.effect) - Math.abs(a.effect);
      }).slice(0, Math.min(25, rows.length));
  return {
    rows,
    selected,
    selectionRule: useConfirmedSet
      ? (aggregated.scale === 'log2'
          ? 'q ≤ 0.10 and |fold change| ≥ 1.2, capped at 50 features'
          : 'q ≤ 0.10 on the supplied scale, capped at 50 features')
      : `only ${significant.length} FDR-qualified feature(s); top ${selected.length} ranked features used for exploratory pathway mapping`,
    effectScale: aggregated.scale,
    contrast,
    mode,
    inferenceMethod,
    inferencePolicy: '≤500 features: seeded/exact permutation tests; >500 features: analytic Welch/paired t-test or one-way ANOVA for browser-scale performance',
    groupSizes,
    steps: aggregated.steps
  };
}

function subjectFeatureSummary(aggregated, feature, options) {
  const source = aggregated.values.get(feature);
  if (!source) return [];
  const perSubject = new Map();
  for (const [sampleId, row] of aggregated.sampleMeta.entries()) {
    const value = source.get(sampleId);
    if (!Number.isFinite(value)) continue;
    if (!perSubject.has(row.subjectId)) perSubject.set(row.subjectId, []);
    perSubject.get(row.subjectId).push({ value, row });
  }

  const timepoints = naturalOrder([...aggregated.sampleMeta.values()].map((row) => row.timepoint));
  const first = timepoints[0] || '';
  const last = timepoints[timepoints.length-1] || '';
  const out = [];
  for (const [subjectId, entries] of perSubject.entries()) {
    if (options.longitudinal && timepoints.length >= 2) {
      if (timepoints.length > 2) {
        const byTime = new Map();
        for (const entry of entries) {
          const time = numericTime(entry.row.timepoint);
          if (!Number.isFinite(time)) continue;
          if (!byTime.has(time)) byTime.set(time, []);
          byTime.get(time).push(entry.value);
        }
        const slope = linearSlope([...byTime.entries()].map(([time, values]) => ({ time, value: mean(values) })));
        if (!Number.isFinite(slope)) continue;
        out.push({ subjectId, condition: entries[0].row.condition, value: slope });
      } else {
        const baseline = entries.filter((x) => x.row.timepoint === first).map((x) => x.value);
        const endpoint = entries.filter((x) => x.row.timepoint === last).map((x) => x.value);
        if (!baseline.length || !endpoint.length) continue;
        const endpointRow = entries.find((x) => x.row.timepoint === last)?.row || entries[0].row;
        out.push({ subjectId, condition: endpointRow.condition, value: mean(endpoint)-mean(baseline) });
      }
    } else {
      const eligible = last ? entries.filter((x) => x.row.timepoint === last) : entries;
      if (!eligible.length) continue;
      out.push({ subjectId, condition: eligible[0].row.condition, value: mean(eligible.map((x) => x.value)) });
    }
  }
  return out;
}

function topVariableFeatures(aggregated, n = 20) {
  return aggregated.features
    .map((feature) => ({
      feature,
      variance: variance([...aggregated.values.get(feature).values()])
    }))
    .filter((x) => Number.isFinite(x.variance) && x.variance > 0)
    .sort((a,b) => b.variance-a.variance)
    .slice(0,n)
    .map((x) => x.feature);
}

function integrationCandidates(aggregated, layerResult, max = 30) {
  const ordered = [
    ...(layerResult?.selected || []).map((x) => x.feature),
    ...topVariableFeatures(aggregated, 20)
  ];
  return [...new Set(ordered)].slice(0,max);
}

function classifyCorrelationChange(rReference, rComparison) {
  if (!Number.isFinite(rReference) || !Number.isFinite(rComparison)) return 'not_estimable';
  const absRef = Math.abs(rReference);
  const absCmp = Math.abs(rComparison);
  if (Math.sign(rReference) !== Math.sign(rComparison) && absRef >= 0.3 && absCmp >= 0.3) return 'sign_reversal';
  if (absRef < 0.3 && absCmp >= 0.5) return 'gained_in_comparison';
  if (absRef >= 0.5 && absCmp < 0.3) return 'lost_in_comparison';
  if (Math.sign(rReference) === Math.sign(rComparison) && absCmp-absRef >= 0.3) return 'strengthened_in_comparison';
  if (Math.sign(rReference) === Math.sign(rComparison) && absRef-absCmp >= 0.3) return 'weakened_in_comparison';
  return 'modest_change';
}

function analyseCrossOmics(aggregatedByLayer, layers, loadedLayers, options) {
  const pairs = [];
  const layerPairs = [];
  for (let i = 0; i < loadedLayers.length; i += 1) {
    for (let j = i+1; j < loadedLayers.length; j += 1) layerPairs.push([loadedLayers[i], loadedLayers[j]]);
  }

  for (const [layerA, layerB] of layerPairs) {
    const featuresA = integrationCandidates(aggregatedByLayer[layerA], layers[layerA], 30);
    const featuresB = integrationCandidates(aggregatedByLayer[layerB], layers[layerB], 30);
    const conditions = naturalOrder([
      ...aggregatedByLayer[layerA].sampleMeta.values(),
      ...aggregatedByLayer[layerB].sampleMeta.values()
    ].map((row) => row.condition));
    if (conditions.length !== 2) continue;

    for (const featureA of featuresA) {
      const a = subjectFeatureSummary(aggregatedByLayer[layerA], featureA, options);
      const aMap = new Map(a.map((x) => [x.subjectId, x]));
      for (const featureB of featuresB) {
        const b = subjectFeatureSummary(aggregatedByLayer[layerB], featureB, options);
        const matched = b.filter((x) => aMap.has(x.subjectId)).map((x) => ({
          subjectId: x.subjectId,
          condition: x.condition,
          x: aMap.get(x.subjectId).value,
          y: x.value
        }));
        const ref = matched.filter((x) => x.condition === conditions[0]);
        const cmp = matched.filter((x) => x.condition === conditions[1]);
        if (ref.length < 4 || cmp.length < 4) continue;
        const stat = differentialCorrelationPair(
          ref.map((x) => x.x), ref.map((x) => x.y),
          cmp.map((x) => x.x), cmp.map((x) => x.y),
          'spearman'
        );
        pairs.push({
          layerA, featureA, layerB, featureB,
          reference: conditions[0],
          comparison: conditions[1],
          nReference: ref.length,
          nComparison: cmp.length,
          ...stat,
          // A Fisher-z approximation to Spearman's rho is not calibrated
          // for paired longitudinal slopes/repeated observations. Only
          // descriptive delta-r is published for such experiments.
          ...(options.longitudinal ? { pValue: null, z: null } : {}),
          pattern: classifyCorrelationChange(stat.rReference, stat.rComparison),
          qValue: null
        });
      }
    }
  }

  bhAdjust(pairs);
  pairs.sort((a,b) => {
    const aq = Number.isFinite(a.qValue) ? a.qValue : 1;
    const bq = Number.isFinite(b.qValue) ? b.qValue : 1;
    if (aq !== bq) return aq-bq;
    return Math.abs(b.deltaR)-Math.abs(a.deltaR);
  });
  return {
    method: options.longitudinal
      ? 'Descriptive cross-omics Spearman differences for subject trajectories; no browser p/q because longitudinal Fisher-z calibration is unverified.'
      : 'Exploratory candidate-selected Spearman correlation by condition; approximate Fisher z for independent groups and BH across the tested selected pairs (not an unbiased discovery test).',
    testedPairs: pairs.length,
    significantPairs: pairs.filter((x) => Number.isFinite(x.qValue) && x.qValue <= 0.10).length,
    pairs: pairs.slice(0,100)
  };
}

function pathwayField(pathway, path, fallback = null) {
  let value = pathway;
  for (const key of path) value = value?.[key];
  return value ?? fallback;
}

function normaliseReactomePathway(pathway) {
  return {
    id: pathway.stId || pathway.dbId || pathway.id || '',
    name: pathway.name || pathway.displayName || 'Unnamed pathway',
    species: pathway.species || pathway.speciesName || '',
    entitiesFound: Number(pathwayField(pathway,['entities','found'], pathway.entitiesFound ?? 0)),
    entitiesTotal: Number(pathwayField(pathway,['entities','total'], pathway.entitiesTotal ?? 0)),
    pValue: Number(pathwayField(pathway,['entities','pValue'], pathway.pValue ?? NaN)),
    fdr: Number(pathwayField(pathway,['entities','fdr'], pathway.fdr ?? NaN)),
    reactionsFound: Number(pathwayField(pathway,['reactions','found'], pathway.reactionsFound ?? 0)),
    raw: pathway
  };
}

function logChooseStable(n, k) {
  if (!Number.isFinite(n) || !Number.isFinite(k) || k < 0 || n < 0 || k > n) return -Infinity;
  return logGamma(n + 1) - logGamma(k + 1) - logGamma(n - k + 1);
}

function hypergeometricUpperTail(k, population, successes, draws) {
  const M = Math.max(0, Math.round(population));
  const K = Math.max(0, Math.min(M, Math.round(successes)));
  const n = Math.max(0, Math.min(M, Math.round(draws)));
  const observed = Math.max(0, Math.round(k));
  const maxX = Math.min(K, n);
  const minX = Math.max(observed, Math.max(0, n - (M - K)));
  if (minX > maxX || M <= 0) return 1;
  const logs = [];
  for (let x = minX; x <= maxX; x += 1) {
    logs.push(logChooseStable(K, x) + logChooseStable(M - K, n - x) - logChooseStable(M, n));
  }
  const maxLog = Math.max(...logs);
  const sum = logs.reduce((acc, value) => acc + Math.exp(value - maxLog), 0);
  return Math.max(0, Math.min(1, Math.exp(maxLog) * sum));
}

export function applyAssayUniverseBackground(selectedResult, universeResult, selectedSubmitted, universeSubmitted) {
  // Reactome's selected-list response is itself filtered by the selected
  // identifiers. Correcting BH only over that list understates multiplicity.
  // Define the statistical family from *all* pathways hit by the measured
  // assay universe, including pathways with k=0 in the selected feature set.
  const background = universeResult?.pathways || [];
  const selected = selectedResult?.pathways || [];
  const backgroundMap = new Map(background.map((row) => [row.id, row]));
  const selectedMap = new Map(selected.map((row) => [row.id, row]));
  const M = Math.max(0, universeSubmitted - Number(universeResult?.identifiersNotFound || 0));
  const n = Math.max(0, selectedSubmitted - Number(selectedResult?.identifiersNotFound || 0));
  const expectedBackground = Number(universeResult?.pathwaysFound ?? background.length);
  const expectedSelected = Number(selectedResult?.pathwaysFound ?? selected.length);
  const complete = Number.isFinite(expectedBackground) && Number.isFinite(expectedSelected)
    && background.length >= expectedBackground && selected.length >= expectedSelected;
  const validUniverse = M > 0 && n > 0 && n <= M && complete &&
    background.length > 0 && background.every((item) => item.id);
  const inconsistent = selected.some((pathway) => {
    const b = backgroundMap.get(pathway.id);
    return !b || Number(pathway.entitiesFound || 0) > Number(b.entitiesFound || 0);
  });
  const allPValues = [];
  if (validUniverse && !inconsistent) {
    for (const pathway of background) {
      const K = Math.min(M, Math.max(0, Number(pathway.entitiesFound || 0)));
      const k = Math.max(0, Number(selectedMap.get(pathway.id)?.entitiesFound || 0));
      allPValues.push({
        id: pathway.id,
        pValue: hypergeometricUpperTail(k, M, K, n),
        qValue: null
      });
    }
    bhAdjust(allPValues);
  }
  const adjusted = new Map(allPValues.map((row) => [row.id, row]));
  const rows = selected.map((pathway) => {
    const backgroundPathway = backgroundMap.get(pathway.id);
    const result = adjusted.get(pathway.id);
    const K = backgroundPathway && M > 0
      ? Math.min(M, Math.max(0, Number(backgroundPathway.entitiesFound || 0))) : null;
    const k = Number(pathway.entitiesFound || 0);
    return {
      ...pathway,
      assayUniversePValue: result?.pValue ?? null,
      assayUniverseFdr: result?.qValue ?? null,
      assayUniverseEntities: K === null ? null : {
        selectedHits: k,
        selectedMapped: n,
        universeHits: K,
        universeMapped: M
      }
    };
  });
  rows.sort((a,b) => {
    const aq = Number.isFinite(a.assayUniverseFdr) ? a.assayUniverseFdr : 1;
    const bq = Number.isFinite(b.assayUniverseFdr) ? b.assayUniverseFdr : 1;
    if (aq !== bq) return aq - bq;
    return b.entitiesFound - a.entitiesFound;
  });
  const reason = !complete ? 'Reactome response is paginated/truncated: entire tested pathway family was not retrieved.'
    : inconsistent ? 'Selected features are not a valid subset of the declared background pathway hits.'
    : !validUniverse ? 'Insufficient or incoherent measurable assay universe.'
    : null;
  return {
    ...selectedResult,
    pathways: rows,
    assayUniverse: {
      submitted: universeSubmitted,
      mapped: M,
      selectedSubmitted,
      selectedMapped: n,
      hypothesesTested: allPValues.length,
      multiplicityStatus: reason ? 'not_estimable' : 'full_assay_universe_family',
      warning: reason,
      method: 'Hypergeometric assay-universe over-representation; BH applied over every Reactome pathway hit by the measured universe (including pathways with zero selected hits), not just the enriched/returned selected pathways. Pathway selection is exploratory and biologically dependent.'
    }
  };
}

export async function reactomeOverRepresentation(ids, { projectToHuman = true, pageSize = 100 } = {}) {
  const unique = [...new Set(ids.map((id) => normaliseText(id)).filter(Boolean))];
  if (!unique.length) return { pathways: [], summary: null, identifiersNotFound: 0, pathwaysFound: 0, token: null };
  const endpoint = projectToHuman
    ? 'https://reactome.org/AnalysisService/identifiers/projection/'
    : 'https://reactome.org/AnalysisService/identifiers/';
  const response = await fetch(`${endpoint}?pageSize=${pageSize}&page=1`, {
    method: 'POST',
    headers: { 'Content-Type': 'text/plain' },
    body: ['#Identifiers', ...unique].join('\n')
  });
  if (!response.ok) throw new Error(`Reactome API returned HTTP ${response.status}`);
  const json = await response.json();
  return {
    pathways: (json.pathways || []).map(normaliseReactomePathway),
    summary: json.summary || null,
    identifiersNotFound: Number(json.identifiersNotFound ?? 0),
    pathwaysFound: Number(json.pathwaysFound ?? 0),
    token: json.summary?.token || null
  };
}

export function mergeReactomeResults(combined, perLayer) {
  const layerMaps = Object.fromEntries(Object.entries(perLayer).map(([layer,result]) => [
    layer,
    new Map((result?.pathways || []).map((pathway) => [pathway.id, pathway]))
  ]));
  return (combined?.pathways || []).map((pathway) => {
    const layerEvidence = {};
    let supportingLayers = 0;
    for (const layer of LAYERS) {
      const hit = layerMaps[layer]?.get(pathway.id);
      // Only the assay-universe-adjusted FDR can be counted as local
      // supporting evidence. Reactome's default FDR uses a different
      // reference universe and must never silently substitute for it.
      const fdr = hit && Number.isFinite(hit.assayUniverseFdr)
        ? hit.assayUniverseFdr : null;
      layerEvidence[layer] = hit ? {
        fdr,
        reactomeDefaultFdr: Number.isFinite(hit.fdr) ? hit.fdr : null,
        assayUniverseFdr: Number.isFinite(hit.assayUniverseFdr) ? hit.assayUniverseFdr : null,
        entitiesFound: hit.entitiesFound
      } : null;
      if (fdr != null && fdr <= 0.10) supportingLayers += 1;
    }
    return { ...pathway, layerEvidence, supportingLayers };
  }).sort((a,b) => {
    if (b.supportingLayers !== a.supportingLayers) return b.supportingLayers - a.supportingLayers;
    const af = Number.isFinite(a.assayUniverseFdr) ? a.assayUniverseFdr : Number.isFinite(a.fdr) ? a.fdr : 1;
    const bf = Number.isFinite(b.assayUniverseFdr) ? b.assayUniverseFdr : Number.isFinite(b.fdr) ? b.fdr : 1;
    if (af !== bf) return af-bf;
    return b.entitiesFound-a.entitiesFound;
  });
}

export function layerOverlapSummary(metadata, loadedLayers) {
  const subjectSets = Object.fromEntries(loadedLayers.map((layer) => [
    layer,
    new Set(metadata.filter((row) => row.omic === layer).map((row) => row.subjectId))
  ]));
  const layerSubjects = Object.fromEntries(loadedLayers.map((layer) => [layer, subjectSets[layer].size]));
  const pairwise = [];
  for (let i = 0; i < loadedLayers.length; i += 1) {
    for (let j = i+1; j < loadedLayers.length; j += 1) {
      const a = loadedLayers[i];
      const b = loadedLayers[j];
      const overlap = [...subjectSets[a]].filter((id) => subjectSets[b].has(id)).length;
      pairwise.push({ layerA: a, layerB: b, matchedSubjects: overlap, layerASubjects: subjectSets[a].size, layerBSubjects: subjectSets[b].size });
    }
  }
  const allSubjects = [...new Set(loadedLayers.flatMap((layer) => [...subjectSets[layer]]))];
  const allMatched = loadedLayers.length
    ? allSubjects.filter((id) => loadedLayers.every((layer) => subjectSets[layer].has(id))).length
    : 0;
  const coveragePatterns = {};
  for (const subject of allSubjects) {
    const pattern = loadedLayers.filter((layer) => subjectSets[layer].has(subject)).join('+') || 'none';
    coveragePatterns[pattern] = (coveragePatterns[pattern] || 0) + 1;
  }
  return { layerSubjects, pairwise, allMatched, totalSubjectsAnyLayer: allSubjects.length, coveragePatterns };
}

function auditBatchDesign(metadata, loadedLayers, protocol) {
  const perLayer = {};
  const blocking = [];

  for (const layer of loadedLayers) {
    const rows = metadata.filter((row) => row.omic === layer);
    const batches = naturalOrder(rows.map((row) => row.batch));
    const conditions = naturalOrder(rows.map((row) => row.condition));
    const timepoints = naturalOrder(rows.map((row) => row.timepoint));

    if (!batches.length) {
      const reasons = protocol.batchKnown === 'yes'
        ? ['protocol declares known technical batches but no batch labels are present in metadata']
        : [];
      perLayer[layer] = {
        status: reasons.length ? 'incomplete' : 'not_provided',
        batches: [],
        blockingReasons: reasons,
        note: reasons.length
          ? 'Batch-aware validation cannot be completed because the protocol declares known batches but metadata provide none.'
          : 'No technical batch labels were supplied for this layer.'
      };
      if (reasons.length) blocking.push({ layer, reasons });
      continue;
    }

    const missingBatchRows = rows.filter((row) => !row.batch).length;
    if (missingBatchRows > 0) {
      const reasons = [`batch labels are missing for ${missingBatchRows} of ${rows.length} assay row(s)`];
      perLayer[layer] = {
        status: 'incomplete',
        batches,
        batchCount: batches.length,
        missingBatchRows,
        blockingReasons: reasons,
        note: 'Partial batch annotation is not sufficient for an auditable deterministic analysis.'
      };
      blocking.push({ layer, reasons });
      continue;
    }

    if (batches.length === 1) {
      perLayer[layer] = {
        status: 'single_batch',
        batches,
        blockingReasons: [],
        note: 'A single technical batch is represented in this layer; no between-batch adjustment is required.'
      };
      continue;
    }

    const batchConditions = new Map(batches.map((batch) => [batch, new Set()]));
    const batchTimes = new Map(batches.map((batch) => [batch, new Set()]));
    const batchOutcomes = new Map(batches.map((batch) => [batch, new Set()]));
    for (const row of rows) {
      if (!row.batch) continue;
      if (row.condition) batchConditions.get(row.batch)?.add(row.condition);
      if (row.timepoint) batchTimes.get(row.batch)?.add(row.timepoint);
      if (row.outcome) batchOutcomes.get(row.batch)?.add(row.outcome);
    }

    const reasons = [];
    const anyBatchSpansConditions = [...batchConditions.values()].some((set) => set.size > 1);
    const anyBatchSpansTimes = [...batchTimes.values()].some((set) => set.size > 1);
    const anyBatchSpansOutcomes = [...batchOutcomes.values()].some((set) => set.size > 1);
    const categoricalOutcomes = naturalOrder(rows.map((row) => row.outcome));

    if (['groups','time','explore'].includes(protocol.objective) && conditions.length > 1 && !anyBatchSpansConditions) {
      reasons.push('condition is completely confounded with batch');
    }
    if (protocol.longitudinal && timepoints.length > 1 && !anyBatchSpansTimes) {
      reasons.push('timepoint is completely confounded with batch');
    }
    if (protocol.objective === 'outcome' && ['binary','multiclass'].includes(protocol.outcomeType) && categoricalOutcomes.length > 1 && !anyBatchSpansOutcomes) {
      reasons.push('categorical outcome is completely confounded with batch');
    }

    const status = reasons.length ? 'confounded' : 'multiple_batches_adjusted';
    const note = reasons.length
      ? `Inference is blocked because ${reasons.join(' and ')}.`
      : protocol.objective === 'outcome'
        ? 'Multiple technical batches are present and enter each outcome model directly as nuisance covariates.'
        : 'Multiple technical batches are present and are adjusted feature-wise by deterministic OLS residualisation before biological inference.';

    perLayer[layer] = {
      status,
      batches,
      batchCount: batches.length,
      blockingReasons: reasons,
      note
    };
    if (reasons.length) blocking.push({ layer, reasons });
  }

  return { perLayer, blocking };
}

// The endpoint is a property of the independent subject, not of the uploaded
// assay. Duplicated sample-sheet entries must never silently overwrite labels.
function assertConsistentSubjectEndpoints(metadata, protocol) {
  const targets = new Map();
  const groupLabels = new Map();
  for (const row of metadata) {
    if (protocol.objective === 'outcome') {
      let value = null;
      if (protocol.outcomeType === 'survival') {
        const hasTime = row.survivalTime !== '';
        const hasEvent = row.survivalEvent !== '';
        if (hasTime !== hasEvent) {
          throw new Error('Incomplete survival endpoint for subject ' + row.subjectId
            + ': survival_time and survival_event must be provided together.');
        }
        if (!hasTime) continue;
        const time = finiteNumber(row.survivalTime);
        const event = finiteNumber(row.survivalEvent);
        if (!(time > 0) || (event !== 0 && event !== 1)) {
          throw new Error('Invalid survival endpoint for subject ' + row.subjectId
            + ': require survival_time > 0 and survival_event coded 0 or 1.');
        }
        value = time + '|' + event;
      } else if (row.outcome !== '') {
        if (['continuous','count'].includes(protocol.outcomeType)) {
          const numeric = finiteNumber(row.outcome);
          if (!Number.isFinite(numeric) || (protocol.outcomeType === 'count' &&
            (numeric < 0 || !Number.isInteger(numeric)))) {
            throw new Error('Invalid numeric outcome for subject ' + row.subjectId
              + ': counts must be non-negative integers and all numeric outcomes finite.');
          }
          value = String(numeric);
        } else value = row.outcome;
      }
      if (value != null) {
        if (targets.has(row.subjectId) && targets.get(row.subjectId) !== value) {
          throw new Error('Conflicting outcome for subject ' + row.subjectId
            + ' across omics or visits. Each independent subject must have one consistent endpoint.');
        }
        targets.set(row.subjectId, value);
      }
    }

    // Independent-group assignments and longitudinal study arms are subject
    // properties. Paired designs are excluded: a subject may provide two conditions.
    const mustBeSingleArm = (protocol.objective === 'groups' &&
      protocol.designType === 'independent') ||
      (protocol.objective === 'time' && protocol.designType === 'repeated');
    if (mustBeSingleArm && row.condition) {
      if (groupLabels.has(row.subjectId) && groupLabels.get(row.subjectId) !== row.condition) {
        throw new Error('Conflicting condition for subject ' + row.subjectId
          + ': this design declares independent arms rather than crossover or paired observations.');
      }
      groupLabels.set(row.subjectId, row.condition);
    }
  }
}

export async function runDeterministicAnalysis({ files, metadataRows, columnMapping, protocol, dataTypes, identifierTypes = {}, useReactome = true, resolveIdentifiers = true }) {
  const inputManifest = {
    metadata: {
      name: files.metadata?.name || 'metadata',
      rows: metadataRows.length,
      fingerprint: 'fnv1a32:' + hashString(JSON.stringify(metadataRows)).toString(16).padStart(8, '0')
    }
  };
  const covariateColumns = Array.isArray(protocol.covariateColumns) ? protocol.covariateColumns.filter(Boolean) : [];
  const loadedLayers = LAYERS.filter((layer) => files[layer]);
  // A sample sheet may describe more measurements than the uploaded matrices.
  // Scope all cohort counts, contrasts, batches and reports to the analysed omics.
  const metadata = canonicalMetadata(metadataRows, columnMapping, covariateColumns)
    .filter((row) => loadedLayers.includes(row.omic));
  if (!metadata.length) throw new Error('No matching metadata rows for the uploaded omics matrices.');
  assertSampleSheetIdentifiers(metadata);
  const biologicalMetadata = metadata.filter((row) => !['blank','qc'].includes(canonicalSampleType(row.sampleType)));
  if (!biologicalMetadata.length) throw new Error('No biological metadata rows remain after excluding blank/QC injections.');
  const singleOmic = loadedLayers.length === 1;
  if (!loadedLayers.length) {
    throw new Error('At least one transcriptomics, proteomics or metabolomics matrix is required.');
  }
  if (protocol.designType === 'crossover') {
    throw new Error('Crossover designs require period/sequence-aware inference and are not yet implemented. No simplified paired analysis was run.');
  }
  assertConsistentSubjectEndpoints(biologicalMetadata, protocol);

  const conditions = naturalOrder(biologicalMetadata.map((row) => row.condition));
  const timepoints = naturalOrder(biologicalMetadata.map((row) => row.timepoint));
  const requiresConditionContrast = protocol.objective === 'groups' || protocol.objective === 'time';
  if (requiresConditionContrast && conditions.length < 2) {
    throw new Error('This objective requires at least two biological conditions; found ' + conditions.length + '.');
  }
  if (requiresConditionContrast && conditions.length > 2) {
    const supportedIndependentMultiGroup = protocol.objective === 'groups'
      && !protocol.longitudinal
      && protocol.designType === 'independent';
    const supportedLongitudinalMultiGroup = protocol.objective === 'time'
      && protocol.longitudinal
      && protocol.designType === 'repeated';
    if (!supportedIndependentMultiGroup && !supportedLongitudinalMultiGroup) {
      throw new Error('More than two conditions require either an independent multi-group comparison or a repeated-measures longitudinal condition × time model.');
    }
  }
  if (protocol.longitudinal && timepoints.length > 2) {
    const numericTimes = timepoints.map((value) => Number(String(value).match(/-?\d+(?:\.\d+)?/)?.[0]));
    if (numericTimes.some((value) => !Number.isFinite(value))) {
      throw new Error('Longitudinal studies with >2 time points require numeric or numeric-labelled time points (for example 0, 6, 12 or T0, T6, T12).');
    }
  }

  if (protocol.objective === 'outcome') {
    const outcomeType = protocol.outcomeType || 'continuous';
    if (timepoints.length > 1 && !protocol.outcomeTimepoint) {
      throw new Error('Outcome analysis with multiple omics time points requires an explicit outcomeTimepoint. Choose the visit/time whose molecular measurements enter the outcome model.');
    }
    if (protocol.outcomeTimepoint && timepoints.length && !timepoints.includes(protocol.outcomeTimepoint)) {
      throw new Error('Selected outcomeTimepoint is not present in the mapped metadata: ' + protocol.outcomeTimepoint);
    }
    if (outcomeType === 'survival') {
      const hasTime = biologicalMetadata.some((row) => Number.isFinite(finiteNumber(row.survivalTime)));
      const hasEvent = biologicalMetadata.some((row) => Number.isFinite(finiteNumber(row.survivalEvent)));
      if (!hasTime || !hasEvent) {
        throw new Error('Survival outcome analysis requires mapped survival_time and survival_event metadata columns.');
      }
    } else if (!biologicalMetadata.some((row) => row.outcome !== '')) {
      throw new Error('Outcome analysis requires a mapped outcome column with non-empty values.');
    }
  }

  const batchAudit = auditBatchDesign(biologicalMetadata, loadedLayers, protocol);
  if (batchAudit.blocking.length) {
    const details = batchAudit.blocking
      .map(({ layer, reasons }) => layer + ': ' + reasons.join('; '))
      .join(' | ');
    throw new Error('Technical batch confounding prevents identifiable biological inference. ' + details + '. Re-balance the design or provide data in which biological conditions/time points overlap technical batches.');
  }

  const overlap = layerOverlapSummary(biologicalMetadata, loadedLayers);
  const layers = {};
  const aggregatedByLayer = {};
  const predictionAggregatedByLayer = {};
  const predictionFoldIsolatedByLayer = {};
  const adjustments = {};

  for (const layer of loadedLayers) {
    const expected = metadata.filter((row) => row.omic === layer).map((row) => row.assayId);
    const text = await files[layer].text();
    inputManifest[layer] = {
      name: files[layer]?.name || layer,
      bytes: Number(files[layer]?.size || text.length),
      fingerprint: 'fnv1a32:' + hashString(text).toString(16).padStart(8, '0')
    };
    const matrix = matrixFromText(text, expected);
    const foldIsolatedPrediction = protocol.objective === 'outcome'
      ? prepareFoldIsolatedPredictionMatrix(matrix, layer, dataTypes[layer])
      : null;
    const allLayerMetadata = metadata.filter((row) => row.omic === layer);
    const layerMetadata = biologicalMetadata.filter((row) => row.omic === layer);
    const msPrepared = layer === 'metabolomics'
      ? applyMetabolomicsMsQc(matrix, allLayerMetadata, {
          blankMode: protocol.msBlankFilter || 'flag',
          blankFold: protocol.msBlankFold,
          qcRsdFilter: protocol.msQcRsdFilter !== false && protocol.msQcRsdFilter !== 'no',
          qcRsdThreshold: protocol.msQcRsdThreshold,
          driftCorrection: protocol.msDriftCorrection !== false && protocol.msDriftCorrection !== 'no',
          mnarStrategy: protocol.msMnarStrategy || 'none'
        })
      : { matrix, qc: null };
    const qcPrepared = prepareMatrixQc(msPrepared.matrix, layer, dataTypes[layer], layerMetadata);
    if (!qcPrepared.matrix.features.length) {
      throw new Error(layer + ': no features remain after modality-specific QC filtering.');
    }
    const processed = preprocessMatrix(qcPrepared.matrix, layer, dataTypes[layer]);
    const rawAggregated = aggregateTechnicalReplicates(processed, biologicalMetadata, layer);

    let aggregated;
    let adjustment;
    let analysisAggregated = rawAggregated;
    const directIndependentGroup = protocol.objective === 'groups'
      && !protocol.longitudinal
      && protocol.designType === 'independent';
    const longitudinalMixed = protocol.objective === 'time'
      && protocol.longitudinal
      && protocol.designType === 'repeated';

    if (protocol.objective === 'outcome' || directIndependentGroup || longitudinalMixed) {
      const nuisanceEncoder = buildCovariateEncoder(
        [...rawAggregated.sampleMeta.values()],
        covariateColumns,
        { includeBatch: true }
      );
      const residualizedForIntegration = residualizeAggregated(rawAggregated, covariateColumns);
      aggregated = residualizedForIntegration.aggregated;
      adjustment = {
        applied: nuisanceEncoder.columnNames.length > 1,
        method: protocol.objective === 'outcome'
          ? 'direct nuisance adjustment inside each outcome model'
          : longitudinalMixed
            ? 'direct nuisance adjustment inside each random-intercept longitudinal model; residualized copy used only for cross-omics correlations'
            : 'direct nuisance adjustment inside each group feature model; residualized copy used only for cross-omics correlations',
        columns: nuisanceEncoder.columnNames.slice(1),
        covariates: covariateColumns,
        note: nuisanceEncoder.columnNames.length > 1
          ? 'Batch and selected covariates enter each feature model directly.'
          : 'No varying batch or selected covariate required adjustment.'
      };
    } else {
      const adjusted = residualizeAggregated(rawAggregated, covariateColumns);
      aggregated = adjusted.aggregated;
      analysisAggregated = aggregated;
      adjustment = adjusted.adjustment;
    }
    aggregated.qc = {
      ...qcPrepared.qc,
      warnings: [...(qcPrepared.qc.warnings || []), ...(processed.audit?.warnings || [])],
      msQc: msPrepared.qc,
      preprocessingSteps: processed.steps,
      preprocessingAudit: processed.audit,
      pca: layerQcPca(rawAggregated),
      inferenceTier: layer === 'transcriptomics' && dataTypes[layer] === 'raw_counts'
        ? {
            level: 'screening',
            label: 'browser screening model',
            note: 'Raw RNA-seq counts are filtered and library-size normalized in-browser for deterministic screening. For publication-grade differential inference, use the provided DESeq2/limma reference R adapter with the same design and covariates.'
          }
        : {
            level: 'native',
            label: 'browser native model',
            note: 'Inference uses the declared processed scale and the explicit design matrix shown in the result.'
          }
    };
    aggregated.matrixShape = { features: matrix.features.length, retainedFeatures: qcPrepared.matrix.features.length, assays: matrix.assays.length, transposed: matrix.transposed };
    aggregated.replicateGroups = rawAggregated.replicateGroups;
    aggregatedByLayer[layer] = aggregated;
    predictionFoldIsolatedByLayer[layer] = Boolean(foldIsolatedPrediction);
    predictionAggregatedByLayer[layer] = foldIsolatedPrediction
      ? aggregateTechnicalReplicates(foldIsolatedPrediction, biologicalMetadata, layer)
      : rawAggregated;
    adjustments[layer] = adjustment;

    if (protocol.objective === 'outcome') {
      layers[layer] = analyseOutcomeLayer(analysisAggregated, layer, {
        outcomeType: protocol.outcomeType || 'continuous',
        covariateColumns,
        targetTimepoint: protocol.outcomeTimepoint || ''
      });
    } else if (directIndependentGroup) {
      layers[layer] = analyseIndependentAdjustedLayer(analysisAggregated, layer, { covariateColumns });
    } else if (longitudinalMixed) {
      layers[layer] = analyseLongitudinalMixedLayer(analysisAggregated, layer, { covariateColumns });
    } else if (protocol.objective !== 'explore') {
      layers[layer] = analyseLayer(analysisAggregated, layer, {
        longitudinal: protocol.longitudinal,
        paired: protocol.designType === 'paired'
      });
    }

    if (layers[layer]) {
      layers[layer].replicateGroups = rawAggregated.replicateGroups;
      layers[layer].matrixShape = { features: matrix.features.length, retainedFeatures: qcPrepared.matrix.features.length, assays: matrix.assays.length, transposed: matrix.transposed };
      layers[layer].qc = aggregated.qc;
      layers[layer].adjustment = adjustment;
    } else {
      aggregated.matrixShape = { features: matrix.features.length, assays: matrix.assays.length, transposed: matrix.transposed };
      aggregated.replicateGroups = rawAggregated.replicateGroups;
    }
  }

  let exploration = null;
  if (protocol.objective === 'explore') {
    exploration = analyseExploratoryIntegration(
      aggregatedByLayer,
      loadedLayers,
      50,
      protocol.partialOmicsExpected === 'yes'
    );
    if (exploration.error) throw new Error(exploration.error);
    for (const layer of loadedLayers) {
      const result = explorationLayerResult(aggregatedByLayer[layer], exploration, layer);
      result.replicateGroups = aggregatedByLayer[layer].replicateGroups || [];
      result.matrixShape = aggregatedByLayer[layer].matrixShape || { features: aggregatedByLayer[layer].features.length, assays: aggregatedByLayer[layer].sampleMeta.size, transposed: false };
      result.qc = aggregatedByLayer[layer].qc || null;
      result.adjustment = adjustments[layer];
      layers[layer] = result;
    }
  }

  const supervisedIntegration = singleOmic ? null : analyseSupervisedMultiblock(
    aggregatedByLayer,
    layers,
    loadedLayers,
    biologicalMetadata,
    protocol
  );
  const predictionUsesGlobalPreprocessing = loadedLayers.some((layer) =>
    predictionFoldIsolatedByLayer[layer] !== true
  );
  const predictiveOutcome = analysePredictiveOutcome(
    predictionAggregatedByLayer,
    loadedLayers,
    biologicalMetadata,
    protocol
  );
  if (predictiveOutcome?.status === 'ok') {
    predictiveOutcome.preprocessingLeakageRisk = predictionUsesGlobalPreprocessing;
    predictiveOutcome.rawIntensityScreening = loadedLayers.some((layer) =>
      (layer === 'proteomics' && dataTypes[layer] === 'lfq_intensity') ||
      (layer === 'metabolomics' && dataTypes[layer] === 'peak_area')
    );
    if (predictiveOutcome.rawIntensityScreening) {
      predictiveOutcome.intensityPanelAssumption = 'Per-sample median centering requires a consistent predefined assay feature panel, including on future external samples. Missing panels or drifting instrument responses invalidate direct portability.';
      predictiveOutcome.technicalQcCaveat = 'Prediction uses the prespecified assay-local intensity transform, not the global MS/pooled-QC drift, blank, RSD and MNAR processing used for exploratory differential analysis. Review technical QC and run an independent external sensitivity analysis before interpreting biomarker performance.';
    }
    predictiveOutcome.preprocessingPolicy = predictionUsesGlobalPreprocessing
      ? 'Some modalities retain whole-cohort QC, intensity transformation or drift correction upstream of folds: exploratory screening only.'
      : 'No cohort-level feature filtering, feature selection or fitted transformation before folds: only fixed transforms and per-assay normalization are applied (including LFQ/MS raw intensity where supported). The raw-intensity predictive branch intentionally excludes fitted pooled-QC and drift corrections; scientific sensitivity and external validation remain necessary.';
  }

  let crossOmics;
  if (singleOmic) {
    crossOmics = {
      method: 'Single-omics analysis: no cross-omics correlations or multiblock integration are calculated.',
      testedPairs: 0, significantPairs: 0, pairs: []
    };
  } else if (protocol.objective === 'explore') {
    crossOmics = {
      method: 'Differential cross-omics correlation is not applicable in the unsupervised branch; integration is performed by balanced multi-block PCA.',
      testedPairs: 0,
      significantPairs: 0,
      pairs: []
    };
  } else if (protocol.objective === 'outcome') {
    crossOmics = {
      method: 'Outcome-targeted feature models are fitted per omics layer; differential correlation is not used as the primary outcome model.',
      testedPairs: 0,
      significantPairs: 0,
      pairs: []
    };
  } else if (conditions.length > 2) {
    crossOmics = {
      method: 'Differential cross-omics correlations currently require exactly two conditions; omitted for the multi-group omnibus analysis.',
      testedPairs: 0,
      significantPairs: 0,
      pairs: []
    };
  } else if (protocol.designType === 'paired' && !protocol.longitudinal) {
    crossOmics = {
      method: 'Direct cross-omics correlation comparison is not run for paired designs in the current engine.',
      testedPairs: 0,
      significantPairs: 0,
      pairs: []
    };
  } else {
    crossOmics = analyseCrossOmics(aggregatedByLayer, layers, loadedLayers, { longitudinal: protocol.longitudinal });
  }

  const selectedIds = Object.fromEntries(
    loadedLayers.map((layer) => [layer, (layers[layer].selected || []).map((row) => row.feature)])
  );
  let identifierResolution = null;
  if (resolveIdentifiers) {
    identifierResolution = {};
    if (selectedIds.transcriptomics?.length) {
      identifierResolution.transcriptomics = await resolveGeneIdentifiers(selectedIds.transcriptomics, {
        identifierType: identifierTypes.transcriptomics || 'unknown',
        organism: protocol.organism || 'human'
      });
    }
    if (selectedIds.proteomics?.length) {
      identifierResolution.proteomics = await resolveProteinIdentifiers(selectedIds.proteomics, {
        identifierType: identifierTypes.proteomics || 'unknown',
        organism: protocol.organism || 'human'
      });
    }
    if (selectedIds.metabolomics?.length) {
      identifierResolution.metabolomics = await resolveMetaboliteIdentifiers(selectedIds.metabolomics);
    }
  }

  const reactomeIds = { ...selectedIds };
  for (const layer of loadedLayers) {
    if (!identifierResolution?.[layer]) continue;
    const map = new Map(identifierResolution[layer].mappings.map((item) => [item.original, item.resolved || item.original]));
    reactomeIds[layer] = selectedIds[layer].map((id) => map.get(id) || id);
  }
  const combinedIds = loadedLayers.flatMap((layer) => reactomeIds[layer] || []);

  let reactome = null;
  let reactomeError = null;
  if (useReactome && combinedIds.length) {
    try {
      const projectToHuman = true;
      const universeIds = Object.fromEntries(loadedLayers.map((layer) => [
        layer,
        [...new Set(aggregatedByLayer[layer].features.map((id) => normaliseText(id)).filter(Boolean))]
      ]));
      const combinedUniverseIds = [...new Set(loadedLayers.flatMap((layer) => universeIds[layer]))];
      const selectedEntries = await Promise.all([
        reactomeOverRepresentation(combinedIds, {projectToHuman, pageSize: 2000}),
        ...loadedLayers.map((layer) => reactomeOverRepresentation(reactomeIds[layer], {projectToHuman, pageSize: 2000}))
      ]);
      const universeEntries = await Promise.all([
        reactomeOverRepresentation(combinedUniverseIds, {projectToHuman, pageSize: 2000}),
        ...loadedLayers.map((layer) => reactomeOverRepresentation(universeIds[layer], {projectToHuman, pageSize: 2000}))
      ]);

      const combined = applyAssayUniverseBackground(
        selectedEntries[0],
        universeEntries[0],
        [...new Set(combinedIds)].length,
        combinedUniverseIds.length
      );
      const perLayer = Object.fromEntries(loadedLayers.map((layer,i) => [
        layer,
        applyAssayUniverseBackground(
          selectedEntries[i+1],
          universeEntries[i+1],
          [...new Set(reactomeIds[layer] || [])].length,
          universeIds[layer].length
        )
      ]));
      reactome = {
        combined,
        perLayer,
        consensus: mergeReactomeResults(combined, perLayer),
        backgroundPolicy: 'All Reactome pathways observable within the QC-retained assay universe; local hypergeometric tests and BH across the full pathway family.',
        backgroundCaveat: 'If either Reactome response is truncated or the background is inconsistent, assay-universe FDR is not estimable. Reactome default FDR is shown separately as an external-database statistic and is never substituted for assay-universe FDR.'
      };
    } catch (error) {
      reactomeError = error instanceof Error ? error.message : 'Reactome API request failed.';
    }
  }

  return {
    generatedAt: new Date().toISOString(),
    engine: {
      name: 'PMx Explain deterministic multi-omics engine',
      version: MULTIOMICS_ENGINE_VERSION,
      execution: 'browser/local deterministic JavaScript',
      analysisMode: singleOmic ? 'single_omic' : 'multiomics',
      analysedLayers: loadedLayers,
      externalServices: {
        ChEBI: resolveIdentifiers,
        Ensembl: resolveIdentifiers,
        UniProt: resolveIdentifiers,
        Reactome: useReactome
      }
    },
    inputManifest,
    protocol: { ...protocol, covariateColumns },
    metadataSummary: {
      subjects: new Set(biologicalMetadata.map((row) => row.subjectId)).size,
      samples: new Set(biologicalMetadata.map((row) => row.sampleId)).size,
      assays: new Set(biologicalMetadata.map((row) => row.assayId)).size,
      technicalQcAssays: new Set(metadata.filter((row) => !biologicalMetadata.includes(row)).map((row) => row.assayId)).size,
      conditions,
      timepoints,
      overlap,
      batchAudit: batchAudit.perLayer,
      adjustment: adjustments
    },
    layers,
    exploration,
    supervisedIntegration,
    predictiveOutcome,
    crossOmics,
    selectedIds,
    reactomeIds,
    identifierResolution,
    reactome,
    reactomeError
  };
}

export function resultToCsv(rows) {
  const header = ['feature','effect','effect_scale','ci95_low','ci95_high','exponentiated_effect','exponentiated_ci95_low','exponentiated_ci95_high','fold_ratio_if_log2','p_value','q_value','n','n_reference','n_comparison','model'];
  const lines = rows.map((row) => [
    row.feature,
    row.effect,
    row.effectScale ?? '',
    row.ciLow ?? '',
    row.ciHigh ?? '',
    row.exponentiatedEffect ?? '',
    row.exponentiatedCiLow ?? '',
    row.exponentiatedCiHigh ?? '',
    row.foldRatio ?? '',
    row.pValue ?? '',
    row.qValue ?? '',
    row.n ?? '',
    row.nReference ?? '',
    row.nComparison ?? '',
    row.model ?? ''
  ].join(','));
  return [header.join(','), ...lines].join('\n');
}
