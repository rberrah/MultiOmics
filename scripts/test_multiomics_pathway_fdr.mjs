import assert from 'node:assert/strict';
import { applyAssayUniverseBackground } from '../src/lib/multiomics/deterministic.js';

function pathway(id, found, p = 1) {
  return {
    id,
    name: 'Test ' + id,
    species: 'Homo sapiens',
    entitiesFound: found,
    entitiesTotal: 200,
    fdr: p
  };
}
const background = Array.from({length:25},(_,i)=>pathway('R-TEST-'+i,5));
const selected = [pathway('R-TEST-0',5,1e-20)];
const full = applyAssayUniverseBackground(
  {pathways:selected,pathwaysFound:1,identifiersNotFound:0},
  {pathways:background,pathwaysFound:background.length,identifiersNotFound:0},
  10,100
);
assert.equal(full.assayUniverse.hypothesesTested,25);
assert.equal(full.assayUniverse.multiplicityStatus,'full_assay_universe_family');
const row = full.pathways[0];
assert(Number.isFinite(row.assayUniversePValue));
assert(Number.isFinite(row.assayUniverseFdr));
assert(row.assayUniversePValue > 0 && row.assayUniversePValue < 1);
assert(Math.abs(row.assayUniverseFdr / row.assayUniversePValue - 25) < 1e-4,
  'BH must include 24 zero-hit pathways, not only the selected enriched pathway');
assert(row.assayUniverseFdr > row.assayUniversePValue,
  'local assay-universe FDR cannot be the unadjusted Reactome p-value');

const paginated = applyAssayUniverseBackground(
  {pathways:selected,pathwaysFound:1,identifiersNotFound:0},
  {pathways:background,pathwaysFound:100,identifiersNotFound:0},
  10,100
);
assert.equal(paginated.assayUniverse.multiplicityStatus,'not_estimable');
assert.equal(paginated.pathways[0].assayUniverseFdr,null,
  'never output a partial-pathway FDR for truncated Reactome responses');

const invalidSubset = applyAssayUniverseBackground(
  {pathways:[pathway('R-TEST-0',7)],pathwaysFound:1,identifiersNotFound:0},
  {pathways:background,pathwaysFound:background.length,identifiersNotFound:0},
  10,100
);
assert.equal(invalidSubset.pathways[0].assayUniversePValue,null);
assert.equal(invalidSubset.pathways[0].assayUniverseFdr,null);
assert.equal(invalidSubset.assayUniverse.multiplicityStatus,'not_estimable');

const unknownSelected = applyAssayUniverseBackground(
  {pathways:[pathway('R-UNKNOWN',3)],pathwaysFound:1,identifiersNotFound:0},
  {pathways:background,pathwaysFound:background.length,identifiersNotFound:0},
  10,100
);
assert.equal(unknownSelected.pathways[0].assayUniverseFdr,null);
console.log('Reactome measured-universe ORA BH family, truncation and subset integrity: PASS');
