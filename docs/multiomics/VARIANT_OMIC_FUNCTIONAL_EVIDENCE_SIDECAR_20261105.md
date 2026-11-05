# Optional research sidecar: variant → omic discordance → functional hypotheses

**Prototype (2026-11-05), not a clinical or causal variant predictor.**
This addition is intentionally **outside** the frozen Research v1
statistical core and the unchanged `DOWNSTREAM_EXCHANGE_V1`
manifest/sample/feature/matrix contract. It supplies a research-only
Python function to build optional evidence annotations on original
source-qualified per-assay data.

## Biological reasoning with three distinct levels

```text
Original genome allele + build + transcript + evidence of consequence
  → assay-resolved molecular observations in the *same* subject/tissue/landmark
  → comparison to original appropriate control within each assay
  → biologically compatible hypotheses and required orthogonal experiments
  → Omics2PK (separately): compare clinical+PGx+PopPK vs actually measured
    pre-landmark omic panels on independent subjects' *measured* PK
```

**Observations are not explanations.** Valid mRNA abundance and
normal peptide abundance can coexist with absent enzymatic activity.
Conversely, normal mRNA with low protein can suggest abnormal translation
or instability, but an actual **protein degradation mechanism** requires
direct turnover/degradation or rescue experiments. Short-read gene totals
can coexist with loss of a **canonical functional splice isoform**. A
variant with high transcript and protein counts may be missense and
nonfunctional because it affects active-site chemistry, substrate
binding, trafficking, localization or partner interactions; classifying
it as posttranslational degradation would be an error.

| Qualifying original evidence | Hypothesis returned | Additional decisive experiments |
|---|---|---|
| Total RNA reference-like; protein decreased | Translation or protein instability **possible** | Pulse-chase/proteasome rescue, targeted full-length and variant peptide |
| Total RNA and protein reference-like; substrate-specific activity decreased | Catalysis, localisation or binding **possible** | Substrate-specific function, proteoform, localization |
| Total RNA reference-like; canonical splice RNA decreased | Aberrant splicing **possible** | Long-read, junction-specific RNA, independent functional activity |
| Stop gain; RNA total reference-like | NMD escape or isoform context **possible** | Same-transcript stop/junction geometry and allele-specific RNA |
| No mutation detected in qualified limited panel | No variant on *assayed sites* | More comprehensive genotype if needed; activity can still vary |
| Genotype untyped | No mechanism claim from an asserted variant | Original VCF / panel and independent functional evidence |
| An endogenous metabolite changes | Network/pathway proxy | Chemical confirmation, multiple enzymes, isotope flux as applicable |
| Protein abundance reference-like but thiopurine *cell* tolerance decreased | Variant may alter drug response independently of abundance | Orthogonal enzyme substrate assay, drug uptake and cell-viability controls; **never call cell PD an enzyme activity measurement** |

All proposed mechanisms are **noncausal hypotheses**.
These are not automated “pathogenic/benign” variant scores.
Functionally classified CYP3A5 *3/*6/*7 and known variants
are external literature priors, not clinical effects re-inferred from
RNA abundance or protein-level SHAP.

## Research interface

`scripts/variant_omic_evidence_sidecar.py` exports the
`reconcile_variant_omic_evidence(case)` function. Inputs require:

- Original gene, actual genotype assay coverage, exact variant and
  transcript/source URL **only if the variant was truly called**.
- For every assayed layer, a source feature, technology, explicit
  readout (gene total, canonical/aberrant RNA isoform, total/intact
  protein, metabolite, *orthogonal* substrate-specific function).
- A qualitative state (`reference_like`, `decreased`, `increased`,
  `absent`, `unknown`) already **validated inside its own assay**
  against a prequalified comparison group. Raw RNA count, MS intensity,
  and protein abundance are **not cross-assay comparable scales**.
- Matching original biological unit, tissue and collection landmark.
  Different donors/organs or predose-vs-postdose assays **abstain**,
  not join on matching gene name. Explicit `original_pairing_scope` is
  mandatory: `same_original_specimen` or
  `same_construct_variant_independent_assays`. The latter allows only
  **construct-level** reconciliation; it cannot be upgraded to same-donor
  clinical biology.
- Optional independent literature reference for a known functional
  variant. **Merely listing an evidence URL is no substitute for
  external source audit.**

Outputs retain observed states, mechanism hypotheses,
uncertainty/abstentions, experiments necessary to confirm mechanism,
and explicit `protein_decay_observed_directly=False` and
`in_vivo_drug_PK_or_PD_improvement_estimated=False`.

```bash
python -m unittest discover -s scripts -p 'test_variant_omic_evidence_sidecar.py' -v
```

Tests are synthetic and cover missense/normal RNA+low protein,
missense/normal RNA+protein but low function, splice-specific RNA,
stop/NMD-unknown, genotype untyped or panel-negative, different
biological subject/sample times, original feature/method validity and
metabolite nonspecificity.

## Original public-dataset evidence (computed in Omics2PK, not fabricated by this sidecar)

**CYP2C9 [original public n4,421 missense variant/abundance/activity benchmark](https://github.com/rberrah/Omics2PK/blob/feat/cyp3a5-functional-evidence-triage-20261103/docs/TWO_WAY_GENETIC_VARIANT_OMICS_TO_PK_PD_PUBLIC_VALIDATION_20261105.md):** 400 independently source-classified missense variants have *decreased activity* despite *wild-type-like abundance*. These assays share the variant identity but were measured in yeast (activity Click-seq) versus cultured human cells (abundance VAMP-seq), **not** paired patient tissues.

**NUDT15 [original public MaveDB n2,825 strict score-paired missense variants](https://github.com/rberrah/Omics2PK/blob/feat/cyp3a5-functional-evidence-triage-20261103/docs/NUDT15_PUBLIC_MAVEDB_VARIANT_STABILITY_THIOPURINE_PD_20261105.md):** 131 source variants show protein-stability score ≥0.8 but thiopurine cellular-tolerance score below 0.3; **another 222** show stability below 0.3 despite cellular tolerance ≥0.8, at exploratory *nonclinical* thresholds. The original 3,015 MaveDB rows include 170 incomplete and 20 nonstandard/uninterpretable protein-ID records. This is a **cellular PD assay**, not a specific hydrolase activity measurement, human treatment response or proof of protein decay.

**CYP3A5 [actual GEO source-metadata pairing](https://github.com/rberrah/Omics2PK/blob/feat/cyp3a5-functional-evidence-triage-20261103/docs/CYP3A5_GEO_HEPATOCYTE_RNA_GENOTYPE_SOURCE_CROSSWALK_20261105.md):** 60 original hepatocyte RNA samples have matching genotype-source sample titles within 96 genotype records. Original `rs776746` SNP calls and exon/isoform quantities are **not independently quantified yet**; published exon-QTL ≠ our own per-sample model. No tacrolimus clinical PK endpoints in these GEO data.

These source results motivate the sidecar and ensure its **hypotheses can be contradicted by real measurements**; the MultiOmics module itself currently has synthetic software-unit tests, not an original cohort-ready/clinically validated variant-effect classifier.

## Complementarity with Omics2PK

**MultiOmics:** describe and reconcile what the original independent
assays measure; generate *a mechanism hypothesis with graded
provenance*, not a % catalytic capacity inferred from SHAP or raw
expression.

**Omics2PK:** physiological drug, tissue and process-specific
qualification, variant codon and splicing annotation, then
**same-subject, outcome-blinded nested CV** to evaluate whether any
additional actual measured layer improves *PK/PD* beyond clinical and
fully assayed authentic PGx. Compare zero/single/multilayer options,
clinically appropriate errors and bias, and an external cohort.
A knowledge prior should be empirically falsifiable; the black-box
model can retain an unknown predictive peak without false identity.

An optional `variant_omic_evidence_sidecar_v1` should be a
**separate versioned attachment** to the exchange, not a modification
to `DOWNSTREAM_EXCHANGE_V1`, and should contain no participant
rows in the public repository. Site integration, authentic
cohort-level calibration and Omics2PK consumption will require a
separately qualified contract PR.

### Source and literature context

- [Suiter et al. 2020 NUDT15 paired protein stability and thiopurine
  cellular sensitivity](https://pmc.ncbi.nlm.nih.gov/articles/PMC7071893/).
- [Amorosi et al. 2021 CYP2C9 matched variant abundance and activity](https://doi.org/10.1016/j.ajhg.2021.07.001).
- [Omics2PK experimental two-way validation protocol](https://github.com/rberrah/Omics2PK/blob/feat/cyp3a5-functional-evidence-triage-20261103/docs/TWO_WAY_GENETIC_VARIANT_OMICS_TO_PK_PD_PUBLIC_VALIDATION_20261105.md).
