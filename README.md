# RE1 Silencer Gates

## Findings

- **Exact motif scoring:** all **7,005 of 7,005** published canonical-site scores match after correcting the motif offset to **(90, 110)**.
- **m-value enrichment:** at **20.86**, diagnostic negative-activity calls are enriched **5.24×** in HepG2, **4.74×** in SK-N-SH, **6.95×** in K562 and **1.11×** in GM12878.
- **Enhancer context changes repression:** the SK-N-SH motif-versus-noMotif gap is **0.74** with En02 and **−0.01** with En21; score–repression correlations span **0.57 to 0.28** across **20** contexts.
- **Chain-only classification mixes distinct evidence:** the **390** mouse-pilot elements resolve into **20 retention candidates, 18 nonretention candidates, 63 unresolved and 289 human-weak sites**.

**Comparative genomics and reporter assays for prioritizing RE1 silencers as candidate transcriptional NOT-gate components.**

Gabrielle DeCuir and Vrinda Sehgal · Johns Hopkins Applied Comparative Genomics

## The question and current state

RE1 is a DNA binding site for the transcriptional repressor REST. This project asks whether sequence conservation adds predictive information about human RE1-mediated repression beyond the REST motif score. The source MPRAduo experiment measures reporter activity with a silencer paired to an enhancer. The engineering goal is a defensible candidate list for experimental testing.

The completed analyses establish sequence-verified human coordinates, a mouse comparison pilot, validated motif scores, and a continuous activity dataset. Human assembly conversion retained **12,815 of 12,866** native elements. The analysis dataset preserves **267,440** source observations at **element × cell type × enhancer** resolution, including missing measurements. The mouse classification track is **complete and frozen at pilot scope**; conservation features, predictive models and gate-performance analyses remain open. [Human conversion](results/tables/liftover_hg19_to_hg38.md) · [Analysis dataset](results/tables/analysis_dataset.md)

The assayed cells are GM12878, HepG2, K562 and SK-N-SH. The dataset contains canonical-site elements (`motif`), REST-peak elements without the canonical motif (`noMotif`), and negative controls (`NegCtrl`). All measured results below come from the linked reports; the [report index](results/README.md) distinguishes current evidence from superseded interpretations.

## Results

### Validated motif scoring and the m-value relationship

Scoring was checked before biological interpretation: **7,005 of 7,005 published canonical-site scores match exactly** after correcting the motif offset to **(90, 110), zero-based inclusive**. The fixed m-value is **20.86**. All **8,436** canonical sites were scored; **1,431** have no published canonical-position hit for direct comparison and were scored with uncensored FIMO rather than assigned zero. [Scoring validation](results/tables/scoring_validation.md)

Above-threshold sites show approximately fivefold enrichment for **diagnostic negative-activity calls** in HepG2 and SK-N-SH, and closer to sevenfold in K562; GM12878 changes little. The diagnostic is negative `log2FoldChange` with supplied `padj ≤ 0.01` in any measured enhancer context.

| Cell type | Diagnostic calls / tested below 20.86 | Fraction below | Diagnostic calls / tested at or above 20.86 | Fraction at or above | Diagnostic enrichment (fold) |
| --- | ---: | ---: | ---: | ---: | ---: |
| HepG2 | 433 / 6,207 | 0.07 | 635 / 1,736 | 0.37 | 5.24 |
| SK-N-SH | 435 / 6,207 | 0.07 | 577 / 1,736 | 0.33 | 4.74 |
| K562 | 71 / 6,207 | 0.01 | 138 / 1,736 | 0.08 | 6.95 |
| GM12878 | 197 / 6,207 | 0.03 | 61 / 1,736 | 0.04 | 1.11 |

Source: the [fixed above/below split](results/tables/scoring_validation.md#fixed-abovebelow-split-diagnostic-calls-only). Enrichment ratios are calculated from those integer counts in [threshold_enrichment.md](results/tables/threshold_enrichment.md). The response rises broadly, with **no clear abrupt step exactly at 20.86**.

Continuous analysis is authorized by the recorded override in `data/processed/scoring_validation/gate.json`; scoring and sign checks remain required. The [scoring report](results/tables/scoring_validation.md) documents the decision.

### Repression depends strongly on enhancer context

**An RE1's effective dynamic range is a property of the silencer–enhancer pair rather than of the silencer alone.** In SK-N-SH, the motif-versus-noMotif median gap is **0.74** against En02 and **−0.01** against En21, in log2 activity units. The gap is `median(noMotif activity) − median(motif activity)`; a positive value means the motif population has lower activity. Averaging enhancers would hide this difference. [Context measurements](results/tables/analysis_dataset.md)

Within canonical elements, Spearman correlation between FIMO score and matched repression ranges from **0.57 to 0.28** across the **20** cell/enhancer combinations. The highest-ranked contexts are:

| Rank | Cell / enhancer | Measured motif elements | NegCtrl baseline | Motif-versus-noMotif gap | Spearman ρ |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | SK-N-SH / En02 | 7,602 | 0.02 | 0.74 | 0.57 |
| 2 | SK-N-SH / En11 | 6,627 | 0.72 | 0.53 | 0.54 |
| 3 | HepG2 / En02 | 7,602 | 0.10 | 0.86 | 0.53 |
| 4 | HepG2 / En11 | 6,627 | 0.65 | 0.76 | 0.51 |
| 5 | K562 / En02 | 7,602 | −0.02 | 0.56 | 0.51 |
| 6 | HepG2 / En19 | 7,763 | 2.37 | 0.44 | 0.51 |
| 7 | K562 / En11 | 6,627 | 1.19 | 0.43 | 0.50 |
| 8 | GM12878 / En09 | 7,890 | 1.89 | 0.18 | 0.48 |

Baseline and gap are in log2 activity units; correlation is unitless. The [complete ranked table](results/tables/analysis_dataset.md#ranked-contexts) includes every context, baseline spread, pair counts and Pearson correlations. GM12878 has a positive continuous score response despite its weak binary diagnostic enrichment; the endpoints answer different questions.

For sequence-response analysis, carry forward the strongest context within each cell: SK-N-SH/En02, HepG2/En02, K562/En02 and GM12878/En09. Also retain En19 separately in every cell for active-enhancer follow-up: En02 is the source study's background-control enhancer and has a near-zero baseline. [Recommendation and baseline evidence](results/tables/analysis_dataset.md#context-recommendation)

### Chain lifting alone contaminates ortholog classification

The mouse pilot exposes distinct sources of contamination in chain-based present/absent calls:

- The element can lift while its motif region does not: **8** pilot elements have element-only mappings. [Independent motif track](results/tables/pilot_mouse_orthologs.md)
- A returned interval can bridge missing source bases or target insertions through the motif. The evidence audit records **112** disrupted-motif alignments. [Alignment evidence](results/tables/pilot_mouse_evidence.md)
- The human canonical site may already be below the binding-score threshold: **289** pilot elements are human-weak. [Human starting-site scores](results/tables/pilot_mouse_evidence.md)

The evidence-aware partition retains all **390** pilot elements:

| Group | Elements | Interpretation and use |
| --- | ---: | --- |
| Retention candidates | 20 | Strong human and corresponding mouse cores with intact, eligible correspondence; prioritize candidate review. |
| Nonretention candidates | 18 | Strong human core and an eligible corresponding mouse core below threshold; review as a sequence contrast. |
| Unresolved | 63 | Human-strong cases with insufficient motif correspondence or scoring evidence, including disrupted or regional alignments; retain the uncertainty. |
| Human-weak | 289 | Human canonical core below 20.86; retain for activity modeling but exclude from the strong-human retention contrast. |

Source: [corrected pilot comparison](results/tables/pilot_mouse_evidence.md). There are no source-unscorable elements in this reported partition. A naive two-way call would mix missing alignment, gapped correspondence and weak starting motifs into a purported lineage-specific class, inflating evidence for human-specific silencing.

Flank anchoring recovered **9** eligible regions from **157** unmapped elements and found **1** above-threshold motif in whole-region scans (**0** in midpoint scans), so it rescued little. [Flank-rescue report](results/tables/pilot_mouse_rescue.md)

## Dataset and repression sign

**Repression = median NegCtrl activity in the same cell and enhancer − element activity. Positive means stronger repression relative to that matched baseline.** A positive raw `log2FoldChange` can still represent repression. Negative repression means activity above the matched baseline; missing activity remains missing.

The [dataset report](results/tables/analysis_dataset.md#repression-sign-and-example) gives this measured example for `REST_chr10:100059172-100059372_Ref` in HepG2/En11:

```text
matched baseline: 0.65; raw activity: 0.34; repression: +0.30
```

Repression is calculated before rounding; displayed measurements use two decimals, with full precision in the linked reports. The dataset keeps supplied activity and padj, population, verified hg38 coordinates where available, canonical-core FIMO score, the threshold indicator and explicit missingness flags. It asserts source row count and exact element/cell/enhancer key preservation; no context is synthesized and no activity is imputed. Controls and unverified native intervals remain present with missing working coordinates. Scores for noMotif and NegCtrl are missing because a canonical site is not defined, not because they were scored below threshold. [Source accounting and missingness](results/tables/analysis_dataset.md)

## Design decisions and why

| Decision | Reason and evidence |
| --- | --- |
| Resolve hg19, then lift to hg38 with sequence verification | Published coordinates match hg19; complete library oligos are checked against the lifted hg38 reference. Only exact verified intervals enter genomic joins. Original IDs and source coordinates remain traceable; never parse working coordinates from IDs. [Assembly](results/tables/coordinate_assembly.md), [conversion](results/tables/liftover_hg19_to_hg38.md) |
| Lift the motif window independently of the full element | Element mapping can succeed without the motif. Separate tracks expose that failure and retain uncertain cases. [Motif-window report](results/tables/pilot_mouse_orthologs.md) |
| Inspect chain blocks and use flanks for regional anchoring | Bounding intervals conceal internal gaps. Flanks locate candidate syntenic regions; chain blocks resolve base-level motif correspondence. [Evidence audit](results/tables/pilot_mouse_evidence.md), [rescue](results/tables/pilot_mouse_rescue.md) |
| Keep element × cell × enhancer observations and matched controls | Baselines and motif effects differ across contexts. A universal zero or an enhancer average would change the scientific question. [Dataset](results/tables/analysis_dataset.md) |
| Use continuous conservation as the modeling feature | Continuous phyloP values can test added predictive information without turning alignment uncertainty into an evolutionary-age label. Track selection and extraction are unrun; the frozen classification track remains a documented result, not the target variable. [Dataset boundary](results/tables/analysis_dataset.md#coordinates-motif-scoring-and-missingness) |

Coordinates use zero-based, half-open BED intervals. The motif offset uses inclusive endpoints and is oriented with the library sequence; it must be converted explicitly when writing BED. LiftOver `minMatch` measures source mapping coverage, not sequence identity or REST binding. [Coordinate and mapping contract](results/tables/pilot_mouse_orthologs.md)

## Pipeline and implementation status

Paths are relative to the repository root. **Complete** means implemented and run for the stated scope. **Helper only** means tested mechanics exist but no production analysis is claimed. The mouse branch is complete and frozen.

| Stage, in execution order | Status | Module / entry point | Make target | Principal output | Report |
| --- | --- | --- | --- | --- | --- |
| Source discovery and schema inspection | Complete; inventory retained | Recorded inspection; source readers in `src/assembly_check.py` | None | Cached GEO tables and schema inventory | [Discovery inventory; state superseded](results/tables/week1_discovery.md) |
| Resolve source assembly | Complete | `src/assembly_check.py` | `assembly-check` | `results/tables/coordinate_assembly_comparisons.csv` | [Assembly](results/tables/coordinate_assembly.md) |
| hg19 → hg38 and sequence verification | Complete | `src/lift.py` | `lift` | `data/processed/elements_hg38.csv` | [Human conversion](results/tables/liftover_hg19_to_hg38.md) |
| Select canonical chr8 pilot and verify references | Complete; frozen pilot scope | `src/pilot.py` | `pilot` | `data/processed/pilot_chr8.bed` | [Pilot setup](results/tables/pilot_setup.md) |
| Element lift to mm39 and mapping sweep | Complete; frozen | `src/pilot_mouse.py` | `pilot-mouse` | `data/processed/pilot_mouse/final/lifted.mm39.bed` | [Mouse lift](results/tables/pilot_mouse_lift.md) |
| Independent motif-window lift and three-way mapping calls | Complete; frozen | `src/pilot_orthologs.py` | `pilot-orthologs` | `data/processed/pilot_mouse/orthologs.mm39.csv` | [Motif coverage](results/tables/pilot_mouse_orthologs.md) |
| Oriented mouse sequence extraction | Complete; frozen | `src/pilot_orthologs.py` | `pilot-orthologs` | `data/processed/pilot_mouse/orthologs.mm39.fa` | [Extraction](results/tables/pilot_mouse_orthologs.md) |
| Mouse FIMO and raw threshold calls | Complete; frozen | `src/pilot_fimo.py`, `src/motifs.py` | `pilot-fimo` | `data/processed/pilot_mouse/fimo_scores.csv` | [Raw FIMO](results/tables/pilot_mouse_fimo.md) |
| Flank-anchored regional rescue | Complete; frozen pilot result | `src/rescue.py` | `pilot-rescue` | `data/processed/pilot_mouse_rescue/rescue.resolved.csv` | [Regional evidence](results/tables/pilot_mouse_rescue.md) |
| Human strength, chain correspondence and comparison eligibility | Complete; frozen | `src/pilot_evidence.py`, `src/alignment_evidence.py`, `src/classify.py` | `pilot-validate` | `data/processed/pilot_mouse_evidence/calls.csv` | [Final pilot comparison](results/tables/pilot_mouse_evidence.md) |
| Full canonical scoring and activity/sign audit | Complete; continuous override recorded | `src/scoring_validation.py` | `scoring-validation` | `data/processed/scoring_validation/canonical_scores.csv`, `gate.json` | [Validation and limitation](results/tables/scoring_validation.md) |
| Context-matched continuous dataset | Complete | `src/dataset.py` | `dataset` | `data/processed/analysis_dataset.csv`, `analysis_dataset_contexts.csv` | [Dataset and context ranking](results/tables/analysis_dataset.md) |
| Remote conservation features and dataset join | Helper only; production extraction not started | `src/conservation.py`; track URLs unset | None | No production feature table | No production report |
| Additional species: rheMac10, canFam6, bosTau9, monDom5 | Analysis not started; shared lift/rescue helpers exist | `src/orthologs.py`, `src/rescue.py`; mouse orchestration is not generalized | None | No validated comparative table | No report |
| Cofactor ChIP-seq features | Not started | No production module or selected feature set | None | None | No report |
| Chromosome/CD-HIT grouped folds | Helper only; production folds not constructed | `src/folds.py` | None | `results/folds/assignments.csv` is a configured destination, not a completed result | No production report |
| Motif-only, motif+conservation, motif+cofactor, gapped-k-mer models | Training not started; explicit stubs | `src/models.py` | None | No fitted models or held-out predictions | No report |
| Gate dynamic range, leakiness and separation | Not started; explicit stubs | `src/metrics.py` | None | No gate metrics | No report |
| Variant robustness | Not started; explicit stub | `src/metrics.py::variant_robustness` | None | No robustness analysis | No report |

`pilot-validate` rebuilds pilot selection, element lifting, motif lifting and extraction, then runs mouse scoring, rescue and evidence publication in order. Human conversion is an explicit prerequisite, not silently rerun. The continuous dataset consumes validated scores and activities independently of mouse mapping eligibility. `classify.py` supplies conservative state mechanics.

The legacy `src/mpra_schema.py` helper assumes complete element/cell pairs, and the generic model evaluator requires binary labels and lacks enhancer keys. Neither is used by `dataset.py`. Production modeling must preserve context and missingness and adapt evaluation to the continuous outcome.

## Open analyses and candidate priorities

Prioritize retention candidates for sequence-and-activity review, nonretention candidates for matched contrasts, and unresolved cases for correspondence review. Human-weak sites remain part of repression modeling and a useful score contrast. These priorities follow the [pilot evidence categories](results/tables/pilot_mouse_evidence.md).

The central modeling work is to select and validate conservation tracks, join measured features to the continuous dataset, construct leakage-safe folds, and compare the stated model families within separate enhancer contexts. Additional species and cofactor features remain open. Gate metrics require defined ON/OFF conditions and measurement scales; variant robustness requires validated variant-to-element mapping.

## Reproducing and checking the work

Create the declared environment and run software checks from the repository root:

```bash
conda env create -f environment.yml
conda activate re1-gates
make test
make lint
```

With cached inputs and executables available, the analysis commands are:

```bash
make assembly-check
make lift
make pilot-validate
make scoring-validation
make dataset
```

Use `make TARGET PYTHON=/path/to/python` to select an interpreter. Reference acquisition is explicit and separate: `src/fetch_lift_inputs.py` fetches human conversion inputs, `src/fetch_mouse_scoring_inputs.py` fetches mouse scoring inputs, and `scripts/fetch_genomes.sh` manages manual genome downloads. See the [data catalog](data/README.md) for source locations and provenance. The scoring audit also requires its cached author result tables and score evidence; the analysis targets do not fetch missing inputs.

Software tests cover coordinate contracts, wrappers, mapping uncertainty, score handling, context-specific baselines, missingness, source-key preservation and fold mechanics. External-tool integration tests use synthetic fixtures and can skip when executables are unavailable; a skip is not successful tool verification. Tests of stubs establish explicit failure behavior, not completed models.

The [dataset replay receipt](results/tables/analysis_dataset_reproducibility.json) records byte-identical outputs from paired full builds. The [report index](results/README.md) also links the conversion and pilot replay receipts and describes their scope. These are recorded comparisons with their own inputs and implementations, not guarantees across software versions. Source hashes in manifests can change after documentation or module organization edits; recorded biological outputs must not be presented as a fresh pipeline run.

## Repository organization

| Location | Responsibility |
| --- | --- |
| `src/` | Importable analysis, validation, acquisition helpers and pipeline entry points; explicit stubs are documented as such. |
| `scripts/` | Operational shell commands and standalone rendering of saved results; scientific validation belongs in `src/`. |
| `tests/` | Unit contracts and external-tool integration fixtures. |
| `references/` | Published sequence/score fixtures and provenance. |
| `data/raw/`, `data/processed/` | Local source data, verified sequences and derived observation/evidence tables; do not commit artifacts. |
| `results/tables/` | Reported measurements, limitations and reproducibility receipts. |
| `results/figures/` | Saved mapping/FIMO flow and scoring-validation figures; the flow shows raw mapping/scan states, not the final comparison groups. |
| `results/README.md` | Current versus superseded report index and organization decisions. |

Generated data and most results are ignored by Git. The root README is self-contained for the measured findings; linked local artifacts require the documented caches or reproduction commands in a fresh checkout. The research-design companion is background, not the authority for implementation status.

## Limitations

- **Published binary labels remain unreproduced and unused.** Exact motif scores and identical GEO activities support continuous analysis; they do not repair the missing matched-enhancer Figure 2B labels. Supplied padj is not interpreted as that silencer FDR.
- **MPRA is episomal.** Reporter effects do not establish behavior in chromosomal integration, endogenous chromatin or coupled circuits.
- **These are cell lines, not mature neurons.** SK-N-SH does not establish behavior in neurons, organoids or neural computing systems.
- **The output is candidate evidence, not demonstrated gates.** Switching behavior, dynamic range, leakiness and robustness are unmeasured here.
- **Alignment uncertainty is not biological absence.** A disrupted motif or failed chain projection needs its own interpretation; the ancestral state is not identified by the mouse comparison.
- **Generalization is untested.** The mouse classification pilot does not establish other-species performance, and descriptive context rankings are not held-out model validation.
- **Conservation is observational.** Added predictive value, if demonstrated, would not establish that evolutionary retention causes repression.

## Sources and attribution

The assay and library sequences come from Mouri et al., *Whole-genome functional characterization of RE1 silencers using a modified massively parallel reporter assay*, with processed activities from GEO and supplementary author tables. UCSC supplies references and chains; JASPAR supplies the REST motif used by FIMO. Local reports document the actual inputs, validations and limitations.

[Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903721/) · [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171) · [Author analysis](https://github.com/tewhey-lab/duoREST) · [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) · [JASPAR](https://jaspar.elixir.no/)
