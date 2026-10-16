import assert from 'node:assert/strict';
import {
  resolveMetaboliteIdentifier,
  resolveMetaboliteIdentifiers
} from '../src/lib/multiomics/deterministic.js';

function jsonResponse(value, status = 200) {
  return new Response(JSON.stringify(value), {
    status,
    headers: { 'Content-Type': 'application/json' }
  });
}

// HMDB -> ChEBI must use an explicit UniChem source-to-source mapping.
{
  const seen = [];
  const fetchFn = async (input) => {
    const url = String(input);
    seen.push(url);
    assert.match(url, /unichem\/rest\/src_compound_id\/HMDB0015623\/18\/7$/);
    return jsonResponse([{ src_compound_id: 'CHEBI:34967' }]);
  };
  const result = await resolveMetaboliteIdentifier('HMDB0015623', { fetchFn, identifierType: 'hmdb' });
  assert.equal(result.resolved, 'CHEBI:34967');
  assert.equal(result.status, 'resolved_external_id');
  assert.equal(result.method, 'unichem_hmdb_to_chebi');
  assert.equal(result.sourceDatabase, 'HMDB');
  assert.equal(result.canonicalDatabase, 'ChEBI');
  assert.equal(result.resolutionConfidence, 'unique_external_mapping');
  assert.equal(result.mappingUsedForIntegration, true);
  assert.equal(result.ambiguityForced, false);
  assert.equal(seen.length, 1);
}

// KEGG Compound is no longer routed through the retired UniChem KEGG source.
// A ChEBI search is accepted only when it returns one unique ChEBI accession.
{
  const fetchFn = async (input) => {
    const url = new URL(String(input));
    assert.equal(url.hostname, 'www.ebi.ac.uk');
    assert.match(url.pathname, /chebi\/backend\/api\/public\/es_search/);
    assert.equal(url.searchParams.get('term'), 'C07430');
    return jsonResponse({
      results: [{
        _source: {
          chebi_accession: 'CHEBI:8060',
          name: 'example',
          ascii_name: 'example',
          stars: 3
        }
      }]
    });
  };
  const result = await resolveMetaboliteIdentifier('C07430', { fetchFn, identifierType: 'kegg_compound' });
  assert.equal(result.resolved, 'CHEBI:8060');
  assert.equal(result.status, 'resolved_external_id');
  assert.equal(result.method, 'chebi_exact_kegg_xref_search');
  assert.equal(result.sourceDatabase, 'KEGG Compound');
  assert.equal(result.resolutionConfidence, 'unique_external_mapping');
}

// Multiple ChEBI hits for one KEGG identifier remain ambiguous instead of
// silently choosing the first search result.
{
  const fetchFn = async () => jsonResponse({
    results: [
      { _source: { chebi_accession: 'CHEBI:111', name: 'a', ascii_name: 'a' } },
      { _source: { chebi_accession: 'CHEBI:222', name: 'b', ascii_name: 'b' } }
    ]
  });
  const result = await resolveMetaboliteIdentifier('C99999', { fetchFn, identifierType: 'kegg_compound' });
  assert.equal(result.resolved, null);
  assert.equal(result.status, 'ambiguous');
  assert.equal(result.method, 'chebi_exact_kegg_xref_search');
  assert.deepEqual(result.candidates, ['CHEBI:111', 'CHEBI:222']);
  assert.equal(result.resolutionConfidence, 'ambiguous_not_used');
  assert.equal(result.mappingUsedForIntegration, false);
  assert.equal(result.ambiguityForced, false);
}

// Numeric PubChem CIDs are mapped only when the identifier type explicitly says PubChem.
{
  const fetchFn = async (input) => {
    const url = String(input);
    assert.match(url, /unichem\/rest\/src_compound_id\/3675\/22\/7$/);
    return jsonResponse([{ src_compound_id: 'CHEBI:8060' }]);
  };
  const result = await resolveMetaboliteIdentifier('3675', { fetchFn, identifierType: 'pubchem' });
  assert.equal(result.resolved, 'CHEBI:8060');
  assert.equal(result.method, 'unichem_pubchem_to_chebi');
  assert.equal(result.sourceDatabase, 'PubChem Compound');
  assert.equal(result.canonicalDatabase, 'ChEBI');
}

// InChIKey uses the structure-linked UniChem response and accepts ChEBI only
// when a single ChEBI accession is present.
{
  const inchikey = 'JZCPYUJPEARBJL-UHFFFAOYSA-N';
  const fetchFn = async (input) => {
    const url = String(input);
    assert.match(url, new RegExp(`unichem/rest/verbose_inchikey/${inchikey}$`));
    return jsonResponse([
      { name: 'chebi', name_label: 'chebi', src_id: 7, src_compound_id: ['CHEBI:34967'] },
      { name: 'hmdb', name_label: 'hmdb', src_id: 18, src_compound_id: ['HMDB0015623'] },
      { name: 'pubchem', name_label: 'pubchem', src_id: 22, src_compound_id: ['3675'] }
    ]);
  };
  const result = await resolveMetaboliteIdentifier(inchikey, { fetchFn, identifierType: 'inchikey' });
  assert.equal(result.resolved, 'CHEBI:34967');
  assert.equal(result.status, 'resolved_external_id');
  assert.equal(result.method, 'unichem_inchikey_to_chebi');
  assert.equal(result.sourceDatabase, 'InChIKey');
  assert.equal(result.resolutionConfidence, 'unique_external_mapping');
}

// Multiple structure-equivalent ChEBI identifiers are never silently collapsed to the first hit.
{
  const fetchFn = async (input) => {
    const url = String(input);
    assert.match(url, /unichem\/rest\/src_compound_id\/HMDB9999999\/18\/7$/);
    return jsonResponse([
      { src_compound_id: 'CHEBI:111' },
      { src_compound_id: 'CHEBI:222' }
    ]);
  };
  const result = await resolveMetaboliteIdentifier('HMDB9999999', { fetchFn, identifierType: 'hmdb' });
  assert.equal(result.resolved, null);
  assert.equal(result.status, 'ambiguous');
  assert.deepEqual(result.candidates, ['CHEBI:111', 'CHEBI:222']);
}

// Plain metabolite names continue to use conservative exact ChEBI name matching.
{
  const fetchFn = async (input) => {
    const url = new URL(String(input));
    assert.equal(url.hostname, 'www.ebi.ac.uk');
    assert.match(url.pathname, /chebi\/backend\/api\/public\/es_search/);
    assert.equal(url.searchParams.get('term'), 'linoleic acid');
    return jsonResponse({
      results: [{
        _source: {
          chebi_accession: 'CHEBI:17351',
          name: 'linoleic acid',
          ascii_name: 'linoleic acid',
          stars: 3
        }
      }]
    });
  };
  const result = await resolveMetaboliteIdentifier('linoleic acid', { fetchFn, identifierType: 'name' });
  assert.equal(result.resolved, 'CHEBI:17351');
  assert.equal(result.status, 'resolved');
  assert.equal(result.sourceDatabase, 'metabolite name');
  assert.equal(result.resolutionConfidence, 'exact_unique_match');
}

// Batch mapping reports a full deterministic coverage audit without forcing mappings.
{
  const fetchFn = async (input) => {
    const url = String(input);
    if (url.includes('/HMDB0000001/18/7')) return jsonResponse([{ src_compound_id: 'CHEBI:1' }]);
    if (url.includes('/HMDB0000002/18/7')) return jsonResponse([]);
    if (url.includes('/HMDB0000003/18/7')) return jsonResponse([
      { src_compound_id: 'CHEBI:3' },
      { src_compound_id: 'CHEBI:33' }
    ]);
    return jsonResponse([], 404);
  };
  const result = await resolveMetaboliteIdentifiers(
    ['HMDB0000001', 'HMDB0000002', 'HMDB0000003'],
    { fetchFn, identifierType: 'hmdb', maxQueries: 5 }
  );
  assert.equal(result.resolvedCount, 1);
  assert.equal(result.unresolvedCount, 2);
  assert.equal(result.strictlyUnresolvedCount, 1);
  assert.equal(result.ambiguousCount, 1);
  assert.equal(result.apiErrorCount, 0);
  assert.equal(result.queryLimitCount, 0);
  assert.equal(result.coverageFraction, 1 / 3);
  assert.equal(result.canonicalDatabase, 'ChEBI');
  assert.equal(result.sourceDatabaseCounts.HMDB, 3);
  assert.equal(result.statusCounts.resolved_external_id, 1);
  assert.equal(result.statusCounts.unresolved, 1);
  assert.equal(result.statusCounts.ambiguous, 1);
  assert.equal(result.methodCounts.unichem_hmdb_to_chebi, 3);
  assert.match(result.mappingPolicy, /Ambiguous mappings are never forced/i);
}

console.log('multiomics strict chemical identifier mapping: PASS');
