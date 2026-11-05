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
  not join on matching gene name.
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
