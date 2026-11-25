# Results and report index

Use [analysis_dataset.md](tables/analysis_dataset.md) for continuous activity and context selection, [scoring_validation.md](tables/scoring_validation.md) for score/sign validation and its binary-label limitation, and [pilot_mouse_evidence.md](tables/pilot_mouse_evidence.md) for the frozen mouse comparison result. These address different endpoints and should be read together.

“Current” means the report remains the evidence for its stated scope. It does not assert that the edited source has been rerun. “Superseded” identifies an interpretation or implementation-state claim replaced by the named report; source observations and recorded receipts are retained. Filenames are stable evidence references, not a project schedule.

## Scientific reports and supporting tables

| Report | Coverage | Status | Replacement or interpretive authority |
| --- | --- | --- | --- |
| [week1_discovery.md](tables/week1_discovery.md) | GEO directory, archive schemas, population inventory and unresolved source assembly. | **Superseded for assembly and implementation state**; useful source inventory. Its claims that assembly is unresolved and downstream stages are unbuilt do not describe the repository. | [coordinate_assembly.md](tables/coordinate_assembly.md) resolves the assembly; [liftover_hg19_to_hg38.md](tables/liftover_hg19_to_hg38.md) documents conversion; [analysis_dataset.md](tables/analysis_dataset.md) defines the production observation table. |
| [coordinate_assembly.md](tables/coordinate_assembly.md) | Exact published oligo comparisons at unchanged hg19/hg38 coordinates. | **Current**, scoped to sampled native populations. | No replacement; the conversion report is the full native-library handoff. |
| [coordinate_assembly_comparisons.csv](tables/coordinate_assembly_comparisons.csv) | Individual sampled comparisons and reference-query evidence supporting assembly resolution. | **Current supporting table**. | Interpreted by coordinate_assembly.md; not independent genome-wide validation. |
| [liftover_hg19_to_hg38.md](tables/liftover_hg19_to_hg38.md) | Full native-library sequence reconciliation, lift, exact verification, exclusions and repeat-overlap audit. | **Current**. | No replacement; authoritative coordinate admission record. |
| [pilot_setup.md](tables/pilot_setup.md) | Canonical chr8 selection, reference checks and sequence regression. | **Current for the frozen pilot**. Its heading is retained as part of the saved artifact. | No replacement for selection; final comparison interpretation is in pilot_mouse_evidence.md. |
| [pilot_mouse_lift.md](tables/pilot_mouse_lift.md) | Element-only mapping sweep and width diagnostics. | **Current mapping measurement**; insufficient as a motif-retention classification. | pilot_mouse_orthologs.md adds motif coverage; pilot_mouse_evidence.md governs biological comparisons. |
| [pilot_mouse_orthologs.md](tables/pilot_mouse_orthologs.md) | Corrected motif window, independent element/motif tracks, three-way mapping calls and sequence extraction. | **Current mapping/extraction result**. Raw ABSENT/MAPPED states are not final retention states. | [pilot_mouse_evidence.md](tables/pilot_mouse_evidence.md) supplies human strength and correspondence eligibility. |
| [pilot_mouse_fimo.md](tables/pilot_mouse_fimo.md) | Raw best-hit scores in mapped mouse windows, threshold calls and browser-check rows. | **Current raw scores**; interpreting `absent` as loss is **superseded**. Blank browser observations are not completed manual validation. | [pilot_mouse_evidence.md](tables/pilot_mouse_evidence.md) governs direct motif comparisons. |
| [pilot_mouse_rescue.md](tables/pilot_mouse_rescue.md) | Corrected neither-track-lifted cohort; flank sweep, anchor acceptance, whole-region and midpoint scans. | **Current regional evidence**; a regional hit does not establish original-motif correspondence. | pilot_mouse_evidence.md integrates regional uncertainty. The report itself replaces the unsupported ABSENT_CONFIRMED terminology. |
| [pilot_mouse_evidence.md](tables/pilot_mouse_evidence.md) | Human starting-site strength, base-level chain evidence and the retention/nonretention/unresolved/human-weak partition. | **Current, complete and frozen at mouse-pilot scope**. | Final authority for pilot comparison eligibility; no evidence of ancestral state or demonstrated binding. |
| [scoring_validation.md](tables/scoring_validation.md) | Exact published score agreement, activity sign, diagnostic threshold bins, missingness and unreproduced Figure 2B labels. | **Current audit and limitation**. Its original binary-label blocking verdict is explicitly overridden for continuous analysis; its diagnostic labels remain unverified as published silencer calls. | The override banner and gate.json specify authorization. [analysis_dataset.md](tables/analysis_dataset.md) supersedes enhancer-aggregated summaries as the production outcome definition, not the audit findings. |
| [threshold_enrichment.md](tables/threshold_enrichment.md) | Explicit above/below diagnostic enrichment ratios calculated from the scoring audit’s integer counts. | **Current derived summary**, with the same unverified-binary-label limitation. | No replacement; scoring_validation.md is the source and interpretive authority. |
| [analysis_dataset.md](tables/analysis_dataset.md) | Source-preserving continuous dataset, context-matched control baselines, repression sign, missingness and ranked context correlations. | **Current production dataset specification and measured results**. | No replacement. Conservation extraction and gate behavior are not claimed. |

## Reproducibility receipts

Receipts record measured repeated runs, not promises about any installation. Output hashes were checked during this organization pass without rerunning biological pipelines. Their saved source hashes retain the exact implementation identity of the recorded run.

| Receipt | Coverage | Status and supersession |
| --- | --- | --- |
| [liftover_hg19_to_hg38_reproducibility.json](tables/liftover_hg19_to_hg38_reproducibility.json) | Human conversion output comparison. | **Current saved-artifact receipt**; listed output hashes match. No replacement. |
| [pilot_mouse_reproducibility.json](tables/pilot_mouse_reproducibility.json) | Standalone pilot/element-lift replay. | **Superseded as the complete pilot replay receipt** by pilot_validation_reproducibility.json. Its paired runs matched each other, but saved hg38/mm39 twoBitInfo log files no longer match this receipt. The element-lift measurements remain documented. |
| [pilot_mouse_orthologs_reproducibility.json](tables/pilot_mouse_orthologs_reproducibility.json) | Corrected motif-lift and extraction subset of the full replay. | **Current saved-artifact subset receipt**; listed hashes match. Points to pilot_validation_reproducibility.json; replaces its archived counterpart. |
| [pilot_mouse_rescue_reproducibility.json](tables/pilot_mouse_rescue_reproducibility.json) | Corrected rescue outputs and scan artifacts from the full replay. | **Current saved-artifact subset receipt**; listed hashes match. Points to pilot_validation_reproducibility.json; replaces its archived counterpart. |
| [pilot_validation_reproducibility.json](tables/pilot_validation_reproducibility.json) | Complete corrected mouse pilot, mapping, FIMO, rescue and evidence publication. | **Current saved-artifact receipt**; listed output hashes match. Supersedes the narrower pilot replay as the reproducibility reference for the final comparison result. |
| [analysis_dataset_reproducibility.json](tables/analysis_dataset_reproducibility.json) | Paired full builds of the continuous dataset, context summary, manifest and report. | **Current saved-artifact receipt**; listed output hashes match. No replacement. |

The source headers, MPRA schema organization and scoring-audit location have changed. Accordingly, the pilot, scoring and dataset manifests do **not** certify the edited implementation: their code-input hashes differ while their recorded output hashes still match. In particular, `python -m src.pilot_evidence --check` correctly detects that implementation drift until an explicit pipeline replay. This documentation pass does not rewrite manifests, erase findings or imply a new biological run.

## Figures and preserved evidence

`figures/pilot_mouse_mapping_fimo_flow.{png,svg,pdf}` depicts raw mapping availability followed by raw mouse-window FIMO calls. Read it with pilot_mouse_orthologs.md and pilot_mouse_fimo.md; it does not depict the final evidence-aware comparison partition.

`figures/scoring_validation.{png,svg,pdf}` depicts the diagnostic score-bin endpoint described in scoring_validation.md. It is not a plot of verified published silencer fractions or the continuous matched-baseline outcome.

`archive/` preserves baseline snapshots and their identity manifests. They are superseded as the interpretation of current outputs by pilot_mouse_evidence.md or the scoring override, according to the snapshot's scope. Archive files and their findings remain unchanged. `folds/` is a destination for grouped assignments; no production fold assignment is present.

## Bounded organization decisions

Analysis and scientific validation belong in `src/`, including importable functions and command-line orchestration. `scripts/` contains operational shell utilities and standalone renderers of saved results. No package hierarchy has been restructured, and no classification code has been deleted.

### Moves applied

| Original location | Canonical location | Reason and compatibility |
| --- | --- | --- |
| `scripts/validate_scoring_activity.py` | `src/scoring_validation.py` | This is reusable scientific scoring/sign validation, with tests and a gate manifest, rather than a standalone presentation script. `make scoring-validation` invokes the canonical module. The original module command remains a thin compatibility entry point. |
| `src/io_mpra.py` | `src/mpra_schema.py` | The implementation validates an explicitly mapped, complete element/cell schema. It is not the production context-preserving MPRA reader. The canonical module's docstring states the missingness and enhancer limitations; `io_mpra.py` retains public compatibility imports. Its tests are named `test_mpra_schema.py`. |

### Naming review and deferred renames

These names were checked against implementation and consumers. Proposed names describe existing scope; they are not claims of generalization. Keeping frozen entry points avoids unnecessary import, command and provenance churn.

| Module | Scope that its name can obscure | Proposed clearer name | Decision and reason |
| --- | --- | --- | --- |
| `src/pilot.py` | Canonical pilot selection plus reference verification helpers used outside the pilot. | `pilot_setup.py`; a separate genome-verification module would require a larger split. | **Retain**. Mixed responsibilities make a filename-only rename incomplete; add the frozen-track header and document scope. |
| `src/pilot_mouse.py` | Specifically the element-lift minMatch sweep, not the whole mouse analysis. | `pilot_mouse_lift.py` | **Retain**. Stable frozen CLI and imports; Makefile and pipeline table identify the exact stage. |
| `src/pilot_orthologs.py` | Mouse-specific dual-track mapping and extraction, not general confirmed orthology. | `pilot_mouse_orthologs.py` | **Retain**. Embedded mouse paths and cross-module consumers remain; a rename would not generalize them. |
| `src/pilot_fimo.py` | Mouse mapped-window scanning and raw threshold calls. | `pilot_mouse_fimo.py` | **Retain**. Frozen entry point is referenced by evidence code and recorded commands; scope is explicit in headers and documentation. |
| `src/pilot_evidence.py` | Mouse-specific source-strength and correspondence audit. | `pilot_mouse_evidence.py` | **Retain**. Preserves the evidence CLI and recorded manifest references. |
| `src/orthologs.py` | Generic interval lifting and oriented extraction; it does not certify orthology. | `interval_projection.py` | **Retain**. Public helper imports remain stable; no need to change the mechanics for this documentation task. |
| `src/rescue.py` | Target-parameterized flank-region evidence, not confirmed rescue of a missing motif. | `flank_rescue.py` | **Retain**. Existing CLI already takes the target assembly; scope and uncertainty are clarified without interface churn. |
| `src/classify.py` | Conservative ortholog-state rules, not an activity classifier or validated evolutionary-age model. | `ortholog_states.py` | **Retain**. Shared enum imports are stable, and the frozen-track header distinguishes this result from continuous modeling. |
| `src/lift.py` | Human assembly conversion plus shared liftOver helpers. | `human_lift.py` for the orchestration only | **Retain**. Renaming the whole file would hide its shared consumers; splitting it is outside the bounded change. |
| `src/fetch_mouse_scoring_inputs.py` | Explicitly mouse-specific input acquisition. | No rename needed. | **Retain**. The name is accurate; parameterized helpers elsewhere do not make this downloader general. |
| `scripts/plot_pilot_mouse_flow.py` | Presentation-only rendering from saved mouse mapping and FIMO CSVs. | No rename needed. | **Retain in scripts/**. No scientific scoring or classification is performed here. |

The remaining `src/` names match their implemented responsibilities: `assembly_check`, `alignment_evidence`, `motifs`, `dataset`, `conservation`, `folds`, `models`, `metrics`, `downloads`, `fetch_lift_inputs`, `config` and `_utils`. “Models” and “metrics” are explicitly stubbed in code and in the pipeline table; naming them does not claim training or gate results. The model evaluator also retains a binary-label/element-cell contract and must be adapted before use with continuous enhancer-specific outcomes.

The ortholog-track modules carry a short complete-and-frozen header scoped to the mouse pilot, including the shared scoring, alignment and state helpers. Shared helpers remain reusable; no other species is declared validated. Makefile comment groups follow acquisition, human conversion, frozen classification, scoring validation, continuous data, unimplemented production stages and software checks. Cache cleanup remains separate from analysis.

### Verification

The complete test suite passed after the moves: **192 passed, no skips**, with the available genomics tools on PATH. Ruff passed. Existing tests exercise the canonical module imports; compatibility imports and the recorded scoring command were also checked without running analysis. No data were staged or committed. Scientific outputs were not regenerated for organizational changes.
