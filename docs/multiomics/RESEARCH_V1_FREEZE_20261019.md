# MultiOmics — Research v1 freeze

**Freeze date:** 19 October 2026  
**Software release:** `1.0.0` / **Research v1**  
**Intended immutable label:** `research-v1.0.0`

Research v1 is a freeze of the **research software, evidence contracts and explicit scientific boundaries**. It is not a medical-device release and it is not a declaration that every p-value or biological conclusion produced from arbitrary input data is valid.

## What is frozen

The Research v1 source line includes:

- standalone FR/EN, light/dark beginner and expert workflows;
- transcriptomic, proteomic and metabolomic import/QC plus single- and multi-omics exploration;
- a versioned, reproducible research ZIP with provenance and native-R tables where those tables were actually executed;
- fail-closed native-R confirmatory **preparation** for the supported group, repeated group×time and independent-survival contracts;
- explicit separation between browser exploratory numbers and native-R results;
- public-data qualification paths for Airway, Nutrimouse, Sacurine raw LC-MS, NKI survival and GSE191328 longitudinal data;
- the versioned MultiOmics → Omics2PK exchange contract without claiming that an omics score is a validated PK covariate.

The machine-readable scope is `research-v1-manifest.json`.

## Meaning of “confirmatory”

In Research v1, “confirmatory” means that the software can **prepare and display the complete output of a compatible native-R model under a strict eligibility contract**. It does not mean that MultiOmics has certified the study.

A native-R PASS never establishes on its own:

- adequate power or sample size;
- correct causal design;
- nominal FDR control in an arbitrary study;
- absence of confounding, MNAR missingness or model misspecification;
- biological replication;
- predictive or clinical utility.

Human scientific review therefore remains required for every publication-level conclusion.

## Evidence required for the freeze

The source revision may be called Research v1 only if its release CI verifies:

1. the standalone SvelteKit/JS suite, npm audit, build and browser workflow;
2. the frozen native-R contracts and runtime;
3. Airway source-exact replay;
4. real Nutrimouse MOFA2 and DIABLO qualification;
5. real 12-file Sacurine xcms/QC/blank evidence;
6. native Cox and longitudinal numerical contracts;
7. real external NKI/DMFS Cox and GSE191328 longitudinal wrapper-vs-direct qualification;
8. the downstream Omics2PK exchange contract;
9. release-version, citation and scientific-boundary consistency.

A failing or skipped required item means **the commit is not the Research v1 freeze**.

## Known limitations intentionally retained

### DESeq2 calibration — issue #44

The software reproduces genuine DESeq2 calculations, but null simulations produced study-level false discoveries above the nominal 5% target in several generators. BH arithmetic, simple size-factor estimation, Cook screening and independent filtering have been investigated and do not explain the complete signal in the tested settings; dispersion/Wald-tail behaviour remains generator-sensitive.

**Release decision:** keep native DESeq2 unchanged, retain the evidence, and prohibit any “FDR universally validated at 5%” claim. This known limitation does not invalidate the software execution contract, but it prevents a universal calibration claim.

### Raw LC-MS analytical scope — issue #47

The authentic 12-file Sacurine v3 study has already passed source authentication, native xcms execution, independent numerical QA and a prespecified parameter-sensitivity analysis. The observed dataset itself contains substantial QC warnings.

**Release decision:** raw-MS software execution is qualified. Chemical identification, standards, drift/carryover validation, interlaboratory accuracy and a universally valid feature-exclusion policy remain outside the frozen claim.

### Individual biological signatures → Omics2PK — issue #20

The exchange format and independent consumer are qualified, but no cross-platform signature is declared biologically transferable merely because a CSV satisfies the schema.

**Release decision:** technical interoperability is frozen; biological PK/PD covariate validation remains a separate Omics2PK research problem.

## Version and runtime

Research v1 pins the application release to **1.0.0** while retaining component-level engine versions for provenance. The reference scientific runtime remains R **4.5.2**, Bioconductor **3.22**, with the immutable Bioconductor container digest and package fingerprints recorded in `multiomics-engine/environment/reference-runtime.json`.

Node.js **22.17+** is supported for source use; Node 24 is recommended for current development.

## Immutable identity

The source file cannot safely contain its own Git commit SHA. Therefore the immutable Research v1 identity is the **final merged commit SHA** of the freeze PR, recorded externally in the release note and, where available, the `research-v1.0.0` Git tag. The repository must not call a later commit Research v1 unless it is a byte-for-byte descendant explicitly released as a new version.

## Change policy after freeze

After Research v1:

- scientific-method changes require a new version and qualification evidence;
- dependency/runtime changes require CI requalification;
- documentation-only corrections must not silently broaden scientific claims;
- PMx hosting remains a separate pinned-SHA integration step;
- Omics2PK, MIPD Engine and PMx Explain remain separate repositories and are not copied back into MultiOmics.

Research v1 is therefore a **stable research-software baseline with known, documented limitations**, not an endpoint for future methodological work.
