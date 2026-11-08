# RE1 Silencer Gates — Research Rationale and Study Design

**Scientific companion to the [main operating guide](README.md).**

Gabrielle DeCuir and Vrinda Sehgal · Johns Hopkins Applied Comparative Genomics · Fall 2026

This file retains its original filename so existing links and editor tabs remain useful. It now serves as the research-design document rather than a second installation guide. The operational status lives in [README.md](README.md); file provenance lives in [data/README.md](data/README.md).

> **Status as of 2026-11-03:** input discovery and a sampled source-assembly audit are complete. Conservation effects, model comparisons, gate-quality scores, and candidate rankings have not been estimated. Statements below about those outputs describe planned work.

## Contents

1. [Research question](#research-question)
2. [Why RE1 is relevant to the project](#why-re1-is-relevant-to-the-project)
3. [What the experiment measures](#what-the-experiment-measures)
4. [Populations and units of observation](#populations-and-units-of-observation)
5. [Hypotheses and planned comparisons](#hypotheses-and-planned-comparisons)
6. [Evolutionary retention as evidence](#evolutionary-retention-as-evidence)
7. [Model families and evaluation design](#model-families-and-evaluation-design)
8. [Gate-quality concepts awaiting definitions](#gate-quality-concepts-awaiting-definitions)
9. [Threats to validity](#threats-to-validity)
10. [Decision register](#decision-register)
11. [Deliverables and reporting standards](#deliverables-and-reporting-standards)
12. [References and team](#references-and-team)

## Research question

**Does evolutionary conservation predict the repression associated with human RE1 elements beyond the information contained in their REST motif scores?**

The proposed engineering use is to prioritize human regulatory sequences for subsequent experiments on transcriptional control. The immediate computational objective is narrower: establish trustworthy source measurements and coordinates, derive comparable features, and evaluate whether those features improve predictions under a defensible held-out design.

The project distinguishes three questions that require different evidence:

| Question | Evidence needed | Present status |
| --- | --- | --- |
| Are the source records, sequences, and coordinate systems understood? | File inspection, provenance, sequence comparisons, and explicit validation | Partially established; full Week 1 validation is outstanding |
| Does conservation add predictive information? | Defined outcomes, verified features, suitable controls, and held-out comparisons | Not tested |
| Would a selected element function as a useful genetic circuit component? | Defined gate metrics and experimental characterization in the intended context | Not demonstrated |

A positive answer to the first question does not establish either of the later answers.

## Why RE1 is relevant to the project

RE1 is a REST-binding regulatory sequence associated with transcriptional repression. The source study characterizes RE1 using a modified massively parallel reporter assay, MPRAduo. This provides a basis for asking how sequence and comparative genomic properties relate to measured reporter behavior.

The NOT-gate framing is an engineering analogy: a regulatory input associated with REST-mediated repression can reduce a transcriptional output. Calling a sequence a candidate NOT-gate component does not establish its transfer function, switching threshold, response time, dynamic range, robustness, or behavior in a neuronal system.

The project therefore treats neural biocomputing as a motivation for selecting experimental candidates. It does not claim that reporter assays in the four source cell lines are equivalent to validated operation in neurons or living computational systems.

## What the experiment measures

The local discovery work establishes the presence of construct-level count columns and processed `log2FoldChange`/`lfcSE` values. Native elements appear with five enhancer contexts: **En02, En09, En11, En19, and En21**. These contexts are part of the experiment, not accidental duplicate rows.

| Measurement distinction | Why it matters |
| --- | --- |
| Silencer element versus enhancer/silencer construct | The same native element has multiple experimental contexts |
| Cell type versus enhancer context | Variation between contexts must not be mistaken for cell specificity |
| Raw counts versus reported activity | Count handling and activity normalization are separate operations |
| Reported activity versus a binary repression call | A continuous value does not supply a validated positive/negative label by itself |
| Native sequence versus variant or scrambled construct | Different sequence populations cannot be pooled merely because IDs share text |
| Present source row versus usable measurement | Some rows have missing activity despite an existing element identifier |

The source `log2FoldChange` column is not automatically the final modeling target. Its interpretation, normalization, baseline, and relationship to the intended repression scale must be established. The existing loader performs no transformation to resolve these scientific choices.

## Populations and units of observation

### Native populations

| Source subset | Distinct elements per cell table | Construct rows per cell table | Intended role |
| --- | ---: | ---: | --- |
| Exact `project == "motif"` | 8,436 | 42,180 | Canonical-motif population |
| Exact `project == "noMotif"` | 4,430 | 22,150 | REST-peak-no-canonical-motif contrast population |

Names such as `motif,maf1` and `motif,rare` identify other source categories and are not included by the exact native `motif` filter. The full WGREST tables have 139,325 rows each, substantially more than the native population alone.

### Geographic scope of the assembly finding

The native interval numbers are supported as **hg19**, with unchanged zero-based, half-open coordinates yielding exact full-length matches in 47/47 sampled published oligos. The project target remains hg38. The source-to-target conversion is now complete: 12,815 verified, 47 lifted with sequence mismatches, and 4 failed lifts out of 12,866 inputs. Only the verified hg38 CSV enters downstream genomic analyses; excluded elements remain in a complete audit with repeat-bias diagnostics.

This assembly conclusion does not automatically apply to variant coordinate encodings, the benchmark library, or the separate noncanonical library. Their provenance must be examined on its own terms.

### Intended modeling unit remains a decision

The scaffold's tidy loader expects one row per element and cell type. The source assay also contains an enhancer-context dimension. Before preparing a modeling table, the study must establish whether to retain that dimension, choose a justified context, or use an explicitly defined summary.

No averaging or selection rule has been implemented. In particular, five contexts must not be turned into five independent copies of an element when estimating generalization to new sequences.

## Hypotheses and planned comparisons

| Proposed comparison | Scientific question | Required control or qualification |
| --- | --- | --- |
| Motif score alone versus motif plus conservation | Does conservation add information beyond the binding motif? | Consistent outcome, population, folds, and preprocessing |
| Motif score alone versus motif plus cofactor occupancy | Does verified occupancy information improve prediction? | Cell-type and assembly compatibility of occupancy sources |
| Motif-based models versus a gapped k-mer model | Does broader sequence information add predictive value? | Training-only model selection and compatible held-out evaluation |
| Retention classes versus repression | Is observed evolutionary retention associated with activity? | Motif score and other justified covariates; missing evidence remains explicit |
| K562 versus SK-N-SH | How consistent are measurements or predictions across the two contexts? | Explicit usable-element intersection, matching enhancer definition, and comparable activity direction |
| Native versus sequence-variant behavior | Which elements tolerate changes in sequence? | Verified reference/variant relationships and a defined robustness metric |

The comparison is predictive or associative unless the design establishes a stronger causal claim. Evolutionary retention may correlate with motif strength, genomic context, selection, assay coverage, and other properties. An observed association would not by itself show that age causes stronger repression.

Potential covariates discussed in the original project plan include basal activity, transposable-element overlap, and GC content. These are candidate design considerations; their exact definitions, source files, inclusion rules, and fitted coefficients do not yet exist in the repository.

## Evolutionary retention as evidence

The target assemblies are rheMac10, mm39, canFam6, bosTau9, and monDom5. The implementation classifies an explicitly supplied species-state matrix using the deepest observed retention:

| Deepest retained evidence | Configured label | Interpretation |
| --- | --- | --- |
| Opossum | `therian_retention` | Retention observed across this deeper comparison |
| Mouse, dog, or cow | `eutherian_retention` | Retention observed across at least one of these comparisons |
| Macaque | `primate_retention` | Retention observed in this comparison, without deeper positive evidence |

The labels intentionally avoid precise evolutionary dates. Mouse, dog, and cow are not ordered as successive ages within this rule table.

Absence and missingness also require distinct interpretations:

- `lifted_no_motif` records a mapped and scored sequence without established retention.
- `failed_lift` records a mapping failure, not demonstrated evolutionary absence.
- `lifted_unscorable` records a mapped sequence whose motif status remains unknown.
- With no positive evidence, complete scored absence yields `no_observed_outgroup_motif`; incomplete evidence yields `unresolved`.

A best motif hit anywhere inside a lifted 200 bp interval is not sufficient proof of orthology to the human motif. Mapping quality, target interval length, orientation, and positional correspondence remain part of the scientific validation problem.

## Model families and evaluation design

### Planned model families

| Entry point | Intended features | Implementation status |
| --- | --- | --- |
| `models.baseline` | Motif score | Stub |
| `models.motif_conservation` | Motif score plus verified conservation features | Stub |
| `models.motif_cofactor` | Motif score plus verified cofactor occupancy | Stub |
| `models.gapped_kmer_svm` | Sequence through lsgkm subprocesses | Stub |

Each currently raises `NotImplementedError`. A stub is a visible boundary around an unverified scientific choice, not a successful model that happens to return no result.

### Dependence-aware evaluation

The implemented fold constructor links all chromosomes that share a CD-HIT cluster, including transitive connections, and assigns whole components to folds. This guards against splitting an existing chromosome or reported sequence cluster between training and testing. It does not prove that every similar sequence pair was discovered by greedy clustering.

All measurements, contexts, and derived features attached to an element must respect its held-out assignment. Model selection, preprocessing, imputation, feature selection, and any learned normalization must be designed so that test information is not used during training.

The shared evaluator computes Spearman and AUROC by cell type and fold from supplied predictions. AUROC requires an externally established binary label. A single-class fold has undefined AUROC, and constant/insufficient values can make Spearman undefined; the evaluator reports these as missing with a reason.

### Sign and orientation of the response

The planned Week 1 check requires repressing elements to occupy the **left** tail on the source activity scale. The evaluator's documented scoring direction is larger predicted values for stronger repression or a greater probability of the positive class. These are different representations and must be connected by an explicit, justified modeling contract.

Preserve the source measurement, verify its interpretation against trusted annotations or controls, and document any later transformed target. Distribution shape alone must not be treated as proof that biological labels or sign are correct. No sign correction has been applied in the current repository.

## Gate-quality concepts awaiting definitions

| Concept | Intended question | Definition still needed |
| --- | --- | --- |
| Dynamic range | How far apart are the on and off outputs? | States, units, normalization, and comparison scale |
| Off-state leakiness | How much expression remains in the off state? | Reference expression and normalization |
| Separation | How distinguishable are on/off measurements? | Replicate structure, variability model, and statistic |
| Variant robustness | How stable is behavior under sequence changes? | Variant population, reference mapping, effect scale, and aggregation |

All four functions in `src/metrics.py` remain stubs. No plausible-looking values, provisional rankings, or undocumented formulas stand in for these definitions.

## Threats to validity

| Issue | Consequence if ignored | Current response |
| --- | --- | --- |
| hg19 source interpreted as hg38 | Wrong reference sequence and genomic annotations | Explicit source/target constants and empirical assembly evidence |
| Compound element names | Dropped, duplicated, or incorrectly assigned sequences | 3,074 explicit publisher aliases resolved with one-to-one checks and source-name retention |
| Enhancer contexts collapsed without justification | Outcome definition changes silently | No aggregation implemented |
| Missing activity confused with absence of an element | Population and intersection counts become misleading | Source presence and usable activity must be reported separately |
| Mapping failure called motif absence | Retention classes become biased by mappability | Distinct failure and absence states |
| Incompatible motif scoring settings | Published threshold is applied to a different score scale | Compatibility still requires verification |
| Genomic interval length treated as sequence verification | Incorrect assembly or orientation can pass a length-only check | Published-sequence comparison is a separate evidence layer |
| Features or related sequences leak across folds | Predictive performance can be overstated | Component grouping implemented; future model fitting must honor it |
| Reporter behavior generalized directly to neurons | Candidate utility is overstated | Final claims remain restricted to measured contexts and proposed follow-up |

## Decision register

| Decision | Present state | Evidence needed to close it |
| --- | --- | --- |
| Native source assembly | hg19 supported by exact sequence comparisons | Sampled source audit complete; full mapped set now compared to hg38 library sequences |
| Analysis target | hg38 conversion complete | 12,815 exact-sequence-verified elements; all exclusions audited |
| Canonical sequence-name reconciliation | Complete | All 3,074 missing direct joins resolved through explicit Table S2 aliases |
| Enhancer handling | Undecided | A documented analysis unit consistent with the assay |
| Missing-measurement policy | Undecided | Counts, causes, and effects on eligibility/intersections |
| Sign and repression labels | Unverified | Trusted source definitions, controls, and call criteria |
| Minimum K562/SK-N-SH usable overlap | Unspecified | An explicit criterion chosen for the intended comparison |
| Cross-species `min_match` | Mouse-pilot value 0.5 from a four-value sweep | 227/390 lift; span/failure diagnostics retained; other species need reassessment |
| CD-HIT identity/coverage | Operational defaults of 0.95/0.95 | A documented sensitivity/validation plan before scientific use |
| Gate metrics and model specifications | Unimplemented | Verified definitions and an evaluation plan |

## Deliverables and reporting standards

The work proceeds from input evidence to scientific interpretation:

1. **Input record:** exact files, source URLs, timestamps, sizes, hashes, schemas, and population definitions.
2. **Validation record:** expected versus observed values for the six Week 1 checks, with pass/fail and explicit prerequisites.
3. **Feature record:** verified coordinate conversions, retained failures, motif settings, conservation-track provenance, and species states.
4. **Evaluation record:** frozen element folds, explicit outcomes, predictions, undefined metrics, and comparable model runs.
5. **Candidate record:** only after the earlier stages, transparent ranking criteria and experimental limitations.

No biological result should be presented as established merely because a wrapper executed successfully or a synthetic test passed. Conversely, a failed validation is an informative result that must be retained rather than repaired by silently changing thresholds.

## References and team

The primary experimental source is Mouri K, et al. *Whole-genome functional characterization of RE1 silencers using a modified massively parallel reporter assay*. Cell Genomics. 2023;3(1):100234. [DOI](https://doi.org/10.1016/j.xgen.2022.100234), [full article and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903721/), [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171).

Additional background references retained from the original project plan:

| Reference | Role in the project's background |
| --- | --- |
| Johnson R, et al. *Evolution of the vertebrate gene regulatory network controlled by the transcriptional repressor REST*. Mol Biol Evol. 2009;26(7):1491–1507. [DOI](https://doi.org/10.1093/molbev/msp058) | REST regulatory-network evolution |
| Kagan BJ, et al. *In vitro neurons learn and exhibit sentience when embodied in a simulated game-world*. Neuron. 2022;110(23):3952–3969. [DOI](https://doi.org/10.1016/j.neuron.2022.09.001) | Broader motivation for biological computing |
| Smirnova L, et al. *Organoid intelligence (OI): the new frontier in biocomputing and intelligence-in-a-dish*. Front Sci. 2023;1:1017235. [DOI](https://doi.org/10.3389/fsci.2023.1017235) | Broader motivation for biological computing |

These background citations do not validate this repository's proposed circuit components or supply results for its unfinished analyses.

**Project team:** Gabrielle DeCuir and Vrinda Sehgal.
