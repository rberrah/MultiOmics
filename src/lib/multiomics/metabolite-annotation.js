// @ts-nocheck
import { parseDelimited } from './deterministic.js';

// A user-provided map is an annotation claim, NOT an identity verification.
// It is only consumed by visualization and must never affect statistical fits.
export function parseMetaboliteAnnotations(csv) {
  const { headers, rows } = parseDelimited(csv);
  const norm = (name) => String(name || '').trim().toLowerCase().replace(/[\s-]+/g, '_');
  const featureCol = headers.find((h) => ['feature_id','feature','original_id','peak_id'].includes(norm(h)));
  const chebiCol = headers.find((h) => ['chebi_id','chebi','chemical_id'].includes(norm(h)));
  if (!featureCol || !chebiCol || featureCol === chebiCol)
    throw new Error('Annotation CSV: provide feature_id and chebi_id columns.');
  if (!rows.length) throw new Error('Annotation CSV is empty.');
  const unique = new Map();
  for (const row of rows) {
    const feature = String(row[featureCol] ?? '').trim();
    const raw = String(row[chebiCol] ?? '').trim();
    if (!feature || !/^CHEBI:[1-9]\d*$/i.test(raw))
      throw new Error('Annotation CSV: each feature must have an exact CHEBI:number ID. Offending row: ' + feature + ' / ' + raw);
    const chebi = raw.toUpperCase();
    if (unique.has(feature))
      throw new Error('Annotation CSV: duplicate feature_id ' + feature + '. Resolve contradictory annotations before import.');
    unique.set(feature, chebi);
  }
  return [...unique].map(([feature, chebi]) => ({
    feature, chebi, status: 'user_declared', provenance: 'user_supplied_csv'
  }));
}
