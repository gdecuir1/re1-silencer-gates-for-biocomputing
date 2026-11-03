# RE1 Silencer Gates

**A reproducible computational genomics project for evaluating human RE1 elements as candidate transcriptional NOT-gate components.**

Gabrielle DeCuir and Vrinda Sehgal · Johns Hopkins Applied Comparative Genomics · Fall 2026

**Current stage:** data discovery and source-assembly verification. The processing scaffold is implemented; the full Week 1 validation harness and scientific models are not yet implemented. Documentation status: **2026-11-03**.

> **Read this before using any coordinates.** The native `motif` and `noMotif` intervals in GEO GSE196171 use **hg19**. The project's intended analysis assembly is **hg38**. These are separate settings, and changing an assembly label does not convert coordinates. The original source intervals have not yet been lifted to hg38.

## Contents

1. [Project question and scope](#project-question-and-scope)
2. [Documentation map](#documentation-map)
3. [Current status and verified findings](#current-status-and-verified-findings)
4. [Installation and working environment](#installation-and-working-environment)
5. [Commands that exist today](#commands-that-exist-today)
6. [Repository layout](#repository-layout)
7. [Data model and configuration](#data-model-and-configuration)
8. [Assembly and coordinate conventions](#assembly-and-coordinate-conventions)
9. [Implemented modules and their contracts](#implemented-modules-and-their-contracts)
10. [Classification and missing evidence](#classification-and-missing-evidence)
11. [Cross-validation and model evaluation](#cross-validation-and-model-evaluation)
12. [Testing and reproducibility](#testing-and-reproducibility)
13. [Remaining work and completion criteria](#remaining-work-and-completion-criteria)
14. [Troubleshooting](#troubleshooting)
15. [Sources and attribution](#sources-and-attribution)

## Project question and scope

The central question is whether evolutionary retention and sequence conservation help predict RE1-mediated repression **after accounting for motif score**. The source experiment measures reporter activity for human RE1 elements across GM12878, HepG2, K562, and SK-N-SH. The planned comparative analysis examines retention in rhesus macaque, mouse, dog, cow, and opossum.

The eventual deliverable is a defensible prioritization of candidate regulatory parts for experimental testing. The project does not yet demonstrate functional logic gates, neuronal computing performance, causal effects of evolutionary age, or validated rankings of candidate elements. The available measurements and future model predictions must remain distinguishable.

See the [research design](README%20%282%29.md) for the hypothesis, proposed comparisons, interpretation limits, and experimental context.

## Documentation map

Each README has a defined role. Use the linked source of truth instead of maintaining independent copies of the same facts.

| Document | Purpose | Read it when you need to… |
| --- | --- | --- |
| [README.md](README.md) | Main operating guide | Understand the repository, install the environment, locate code, and run supported checks |
| [README (2).md](README%20%282%29.md) | Research rationale and study design | Explain the question, distinguish hypotheses from results, and assess interpretation limits |
| [data-README.md](data-README.md) | Week 1 acquisition and validation guide | Follow the dependency order, understand unresolved input decisions, and assess readiness |
| [data/README.md](data/README.md) | Authoritative data catalog and download ledger | Locate actual files, inspect source schemas, and recover URLs, timestamps, sizes, and hashes |
| [src/config.py](src/config.py) | Executable configuration | Inspect constants, path defaults, unset mappings, and classification rules |

Local evidence reports are generated artifacts, not guaranteed to exist in a fresh checkout:

| Report | What it establishes |
| --- | --- |
| [week1_discovery.md](results/tables/week1_discovery.md) | GEO listing, seven table schemas, row counts, and the original discovery findings |
| [coordinate_assembly.md](results/tables/coordinate_assembly.md) | Deterministic sequence comparisons distinguishing hg19 from hg38 |
| [coordinate_assembly_comparisons.csv](results/tables/coordinate_assembly_comparisons.csv) | Individual assembly comparisons, query URLs, and response hashes |

## Current status and verified findings

### Implementation status

| Component | Status | Important boundary |
| --- | --- | --- |
| GEO discovery | Completed | Native construct-level activity and interval columns were inspected directly |
| Native source assembly | Resolved empirically as hg19 | Based on a deterministic sequence sample; not an all-element reference audit |
| hg19 → hg38 conversion | Not performed | Required before applying the configured hg38 downstream workflows |
| Generic MPRA loading and validation | Implemented | Raw GEO tables still require scientifically justified preparation |
| FIMO, liftOver, bedtools, conservation wrappers | Implemented | Production inputs and several interpretation choices remain outstanding |
| Ortholog-state classification | Implemented and tested | Classifies supplied states; does not infer states from missing data |
| Chromosome/cluster fold construction | Implemented and tested | No real-data fold assignments have been produced |
| Shared held-out evaluator | Implemented | Does not train models or prove that predictions were held out |
| Four model-training entry points | Explicit stubs | Raise `NotImplementedError` |
| Four gate-quality metrics | Explicit stubs | Raise `NotImplementedError` |
| Full six-check Week 1 validation | Not implemented | `src/validate.py`, `make validate`, and the requested validation report do not exist yet |
| Manual genome acquisition target | Not implemented | `make data-genomes` does not exist yet |

### Findings supported by the downloaded inputs

| Observation | Verified value | Scope |
| --- | ---: | --- |
| Rows per WGREST cell-type table | 139,325 | Includes controls, variants, and five enhancer contexts |
| Distinct native canonical elements | 8,436 | Exact `project == "motif"` subset |
| Distinct native no-motif contrast elements | 4,430 | Exact `project == "noMotif"` subset |
| Enhancer contexts per native element | 5 | En02, En09, En11, En19, En21 |
| Native coordinate definitions shared across cell tables | 12,866 | All four cell tables agree on these interval records |
| Direct Table S2 sequence-name matches | 9,792 | 5,362 canonical + 4,430 no-motif; all matched sequences are 200 bp |
| Canonical names still unreconciled with Table S2 | 3,074 | Not imputed and not treated as reference-sequence mismatches |
| Exact reference matches in the assembly audit | hg19: 47/47; hg38: 10/47 | Full 200 bp sequences at unchanged interval numbers |

These are acquisition and validation findings. They are not conservation effects, predictive performance estimates, or experimental gate-quality results.

## Installation and working environment

### New environment

Run these commands from the repository root:

```bash
conda env create -f environment.yml
conda activate re1-gates
python --version
python -c "import src.config"
make test
make lint
```

The environment specifies Python **3.11**. `make env` is an alternative to the first command. It installs software dependencies; it does not acquire project datasets.

### Existing environment

```bash
conda env update -f environment.yml
conda activate re1-gates
make test
make lint
```

Python dependencies and command-line genomics packages are listed in [environment.yml](environment.yml). Do not assume that the system Python, the notebook kernel, and the activated Conda environment are the same interpreter.

### External executables

| Executable | Used by | Expected source |
| --- | --- | --- |
| `fimo` | Motif scoring | MEME package in the environment |
| `cd-hit-est` | DNA clustering | CD-HIT package in the environment |
| `liftOver` | Interval conversion | UCSC liftOver package in the environment |
| `bedtools` | Oriented reference-sequence extraction | bedtools package in the environment |
| `curl` | Explicit small-file downloads | Available on the system PATH; not declared as a Conda dependency here |
| `twoBitToFa` | Planned genome extraction workflow | UCSC package listed in the environment; acquisition workflow not yet implemented |
| `gkmtrain`, `gkmpredict` | Future gapped k-mer model | Separate lsgkm installation; model function is currently a stub |

```bash
command -v fimo cd-hit-est liftOver bedtools curl
```

For a different configured interpreter, use `make test lint PYTHON=python`. A local `.venv` may exist in a developer's checkout, but it is ignored and is not a portable replacement for `environment.yml`.

## Commands that exist today

| Command | Action | Project data required? | Network access? |
| --- | --- | --- | --- |
| `make env` | Create the Conda environment | No | Software package downloads |
| `make test` | Run pytest | No | None required by the tests |
| `make lint` | Run Ruff | No | None |
| `make clean` | Remove pytest/Ruff caches and Python caches under source, tests, and notebooks | No | None |
| `make assembly-check` | Recompute the assembly evidence report from local inputs and cached UCSC responses | Yes | None |
| `python -m src.assembly_check --fetch` | Obtain missing small assembly-evidence files and run the audit | GEO archive must already be present | Explicit, bounded requests |

**A fresh checkout cannot run the assembly audit from code alone.** The audit requires the GEO archive at `data/raw/geo/GSE196171/GSE196171_RAW.tar`. Its `--fetch` option fetches Table S2 and missing sequence-window responses, but does not fetch that GEO archive. See the [acquisition guide](data-README.md#acquisition-and-replay-procedures).

The following names describe planned deliverables, not available commands: `make validate` and `make data-genomes`.

## Repository layout

```text
.
├── README.md                         Main operating guide
├── README (2).md                     Research rationale and design
├── data-README.md                    Week 1 acquisition and validation guide
├── environment.yml                  Python 3.11 and dependency specification
├── pyproject.toml                   Python version, pytest, and Ruff configuration
├── Makefile                         Supported developer/audit commands
├── .gitignore                       Active Git exclusions
├── gitignore                        Retained legacy file; Git does not apply it
├── data/
│   ├── README.md                    Canonical catalog and append-only download ledger
│   ├── raw/                         Original downloads and cached responses; ignored
│   └── processed/                   Prepared inputs; currently placeholders only
├── src/
│   ├── __init__.py
│   ├── config.py
│   ├── _utils.py
│   ├── downloads.py
│   ├── assembly_check.py
│   ├── io_mpra.py
│   ├── motifs.py
│   ├── conservation.py
│   ├── orthologs.py
│   ├── classify.py
│   ├── folds.py
│   ├── models.py
│   └── metrics.py
├── tests/                           Synthetic unit and optional executable tests
├── notebooks/
│   └── 01_explore.ipynb              Unexecuted exploration starter
└── results/
    ├── tables/                      Discovery and assembly reports; ignored
    ├── figures/                     No scientific figures produced yet
    └── folds/                       Fold definitions may be versioned; none produced yet
```

Under `data/` and `results/`, the active ignore policy allows README files and `.gitkeep` placeholders, plus files under `results/folds/`. Large genomic formats are additionally ignored throughout the repository. No raw dataset is required to run the test suite. An allowed fold-output directory does not mean that real fold definitions have already been generated or committed.

## Data model and configuration

### Constants and their roles

| Setting | Value | Meaning |
| --- | --- | --- |
| `ELEMENT_LENGTH` | `200` | Expected native element-sequence length |
| `MOTIF_SPAN` | `(91, 111)` | Configured zero-based, inclusive motif offsets |
| `M_VALUE` | `20.86` | Configured score threshold; scoring compatibility must still be checked |
| `CELL_TYPES` | GM12878, HepG2, K562, SK-N-SH | Canonical cell labels in prepared inputs |
| `JASPAR_ID` | `MA0138.2` | Required REST matrix identifier by default |
| `GENOME_BUILD` | `hg38` | Intended analysis target |
| `GEO_NATIVE_GENOME_BUILD` | `hg19` | Empirically resolved assembly of the native source sets |
| `EXPECTED_ELEMENT_COUNT` | `8436` | Default loader expectation for the canonical set |
| `TARGET_SPECIES` | rheMac10, mm39, canFam6, bosTau9, monDom5 | Cross-species target assemblies |

Path defaults in `src/config.py` derive from the file's repository location using `pathlib`. The assembly audit also defines its own evidence/report paths relative to these configured directories. No developer-specific absolute path is needed in the source.

### Prepared MPRA table contract

The loader's canonical output has one row per `(element_id, cell_type)`:

| Field | Required? | Contract |
| --- | --- | --- |
| `element_id` | Yes | Nonempty, whitespace-free string; unique within a cell type |
| `sequence` | Yes | Exactly 200 bp of valid IUPAC DNA; consistent across cell types |
| `cell_type` | Yes | One of the four canonical labels |
| `repression` | Yes | Finite numeric value with an explicitly established interpretation |
| `chrom`, `start`, `end`, `strand` | Together, when coordinates are supplied | Valid oriented intervals; source/target assembly must be established outside the loader |

The default validator requires all four cell types for every canonical element. It rejects duplicate element/cell pairs, incomplete measurements, inconsistent metadata, malformed intervals, and incorrect lengths. It does not infer an assembly from interval numbers.

The downloaded GEO tables do **not** directly satisfy this contract: their IDs describe enhancer/silencer constructs, they contain five contexts per native element, some measurements are missing, and they lack sequence and strand columns. The source label `SK.N.SH` also differs from the canonical `SK-N-SH` label.

### Column mappings

| Configuration | Wide input | Long input |
| --- | --- | --- |
| `MPRA_LAYOUT` | `"wide"` | `"long"` |
| `MPRA_COLUMNS` | Map `element_id`, `sequence`, and optional metadata to source headers | Also map `cell_type` and `repression` |
| `MPRA_REPRESSION_COLUMNS` | Map each of the four cell types to its measurement column | Not used |
| `MPRA_SEPARATOR` | Tab by default | Tab by default |

Unset mappings deliberately fail. Do not populate them by renaming the raw construct ID to an element ID or treating `log2FoldChange` as fully validated repression. Enhancer-context handling, sequence reconciliation, missingness, sign, and assembly conversion must be established first.

## Assembly and coordinate conventions

### Evidence for hg19

The audit joins native GEO intervals to the `WG_REST` library in Table S2, then selects one directly matched oligo per available chromosome/project group before examining reference sequences. It tests 24 canonical groups and 23 no-motif groups, allowing exact forward or reverse-complement equality.

| Candidate reference | Exact matches at original coordinates | Matches with start −1 | Matches with start +1 |
| --- | ---: | ---: | ---: |
| hg19 | 47/47 | 0/47 | 0/47 |
| hg38 | 10/47 | 1/47 | 0/47 |

Every hg19 match is forward-oriented. The ten matches shared with hg38 are nondiscriminating regions, not evidence that the entire table can be read as hg38. The audit resolves the source assembly for the native sets, while full-population sequence checks and variant-coordinate provenance remain separate work.

### Coordinate systems must be translated explicitly

| Representation | Convention |
| --- | --- |
| Source native intervals supported by the audit | hg19, zero-based, half-open `[start, end)` |
| Intended downstream human intervals | hg38, zero-based, half-open, after explicit conversion and validation |
| `MOTIF_SPAN` in config | Zero-based, both endpoints inclusive: 21 bases |
| FIMO `start` and `stop` | One-based, inclusive, retained as emitted by FIMO |
| UCSC audit requests | 202 bases: the source 200 bp interval plus one diagnostic base at each end |

For the current configured motif offsets, conservation queries use `[start + 91, start + 112)` on a forward-oriented sequence and `[end - 112, end - 91)` on a reverse-oriented sequence. The relevant strand is the orientation of the element sequence relative to the genome, not the strand of a FIMO hit.

The assembly audit does not validate motif placement. The configured motif offsets and one-based FIMO positions must be compared in a common coordinate system during the future motif-position check; an exact 200 bp sequence match alone cannot establish the motif's location.

## Implemented modules and their contracts

| Module / entry point | Input | Output / behavior | Boundary |
| --- | --- | --- | --- |
| [`downloads.download_small_file`](src/downloads.py) | HTTPS URL and destination below `data/` | Bounded download or provenance-checked cache | Explicit `fetch=True`; 10 MB default cap, maximum 100 MB; not a genome downloader |
| [`assembly_check`](src/assembly_check.py) | GEO archive, Table S2, cached or requested UCSC windows | Markdown report and comparison CSV | Offline by default; native sequence sample only |
| [`io_mpra.load_mpra`](src/io_mpra.py) | Local delimited file and explicit mappings | Validated tidy DataFrame | No subset selection, context aggregation, sign change, or imputation |
| [`motifs.score_sequences`](src/motifs.py) | DNA FASTA and local MEME-format matrix | Best hit per sequence, retaining missing hits | No JASPAR download or genome extraction |
| [`conservation.query_conservation`](src/conservation.py) | Oriented 200 bp intervals and two track URLs | Exact motif-window means and coverage counts | Track URLs/builds must be established; no full bigWig download |
| [`orthologs.lift_elements`](src/orthologs.py) | Full source intervals, chain, explicit `min_match` | Mapped and unmapped records with failure reasons | Target lengths may change; no arbitrary multimapping selection |
| [`orthologs.extract_sequences`](src/orthologs.py) | Target intervals and local target FASTA | Strand-aware sequences | Checks skipped/truncated outputs; bedtools may create an adjacent FASTA index |
| [`classify.classify_element` / `classify_matrix`](src/classify.py) | Explicit per-species states | Retention class with missing evidence preserved | No automatic orthology or absence inference |
| [`folds.build_folds`](src/folds.py) | Unique element IDs, chromosomes, and sequences | Deterministic fold CSV | Entire chromosome/cluster components stay together |
| [`models.evaluate_folds`](src/models.py) | Supplied held-out observations, predictions, and labels | Spearman and AUROC by fold and cell type | No training or label-threshold inference |

FIMO runs in text mode with a p-value threshold of 1 so weak best scores are not removed by its default significance filter. Output is streamed to disk and reduced in chunks. The best score wins, with deterministic tie breaking. `retains_motif` uses `score >= 20.86`; unscorable sequences retain missing scores and missing retention calls. Provide the intended MEME matrix and verify its background/scoring compatibility before interpreting the cutoff scientifically.

Conservation checks HTTP byte-range support before opening remote bigWigs. It uses actual per-base values for the 21 bp window, reports the number of covered bases, and leaves all-missing windows as missing. Missing coverage is not converted to zero conservation.

## Classification and missing evidence

Each of the five target species must have an explicit state:

| State | Interpretation |
| --- | --- |
| `retained` | Motif retention has been established for the mapped element |
| `lifted_no_motif` | Mapping and scoring succeeded, but no retained motif was established |
| `failed_lift` | No mapped interval; biological absence is not established |
| `lifted_unscorable` | Mapped interval exists, but motif presence cannot be evaluated |

The deepest observed retention wins according to `EVOLUTION_RULES`:

| Evidence | Assigned class |
| --- | --- |
| Opossum retention | `therian_retention` |
| Mouse, dog, or cow retention, with no deeper retained species | `eutherian_retention` |
| Macaque retention, with no deeper retained species | `primate_retention` |
| No retention and some failed/unscorable evidence | `unresolved` |
| No retention and complete scored absence in all targets | `no_observed_outgroup_motif` |

These labels describe observed evidence, not divergence dates or proof of human specificity. A retained motif can support a class even when another species is unresolved; the missing evidence remains visible. The matrix output retains the input states and adds `evolutionary_class` and `evidence_complete`.

## Cross-validation and model evaluation

### Fold construction

`build_folds` clusters 200 bp DNA sequences using CD-HIT-EST on both strands, then connects chromosomes that share any cluster. Transitive connections matter: if a cluster connects two chromosomes and another connects one of them to a third, all three remain in one fold.

| Setting | Default | Interpretation |
| --- | ---: | --- |
| Global identity | 0.95 | Operational clustering choice, not a verified biological cutoff |
| Coverage of both sequences | 0.95 | Operational alignment-coverage requirement |
| Number of folds | 5 | Requires at least five independent components |

Whole components are allocated deterministically, with balancing by element count. Too few independent components raises an error. Every chromosome and every reported CD-HIT cluster is confined to one fold. CD-HIT uses representative-based greedy clustering; this does not establish an exhaustive all-pairs similarity guarantee.

The output records `element_id`, `chrom`, `cluster_id`, `component_id`, `fold`, and clustering thresholds. The default destination is `results/folds/assignments.csv`. All cell-type measurements belonging to an element must later inherit that element's fold.

### Evaluation and sign

`evaluate_folds` requires `element_id`, `cell_type`, `fold`, `observed`, `prediction`, and explicit binary `positive` labels. It rejects duplicated element/cell predictions, nonfinite numeric values, and elements assigned to multiple held-out folds. Undefined Spearman correlations and single-class AUROCs are reported as NaN with status fields.

**Direction is a pending scientific contract.** The planned source-data validation requires repression to lie in the left tail of the activity scale. The evaluator currently documents predictions oriented so that larger values mean stronger repression/more likely positive. Those conventions must be reconciled explicitly before modeling: preserve the measured activity, define any modeling transformation separately, and give `observed` and `prediction` compatible directions. No sign transformation or positive-label definition has been validated or applied yet.

The following functions intentionally raise `NotImplementedError`:

- Models: `baseline`, `motif_conservation`, `motif_cofactor`, `gapped_kmer_svm`.
- Gate metrics: `dynamic_range`, `off_state_leakiness`, `separation`, `variant_robustness`.

## Testing and reproducibility

The last recorded implementation verification on 2026-11-03 passed **74 tests**, with Ruff clean. This is a historical verification result, not a substitute for rerunning checks after changes.

| Test area | What it verifies |
| --- | --- |
| Configuration | Fixed constants, portable path resolution, and notebook output hygiene |
| MPRA loading | Required columns, counts, lengths, cell labels, duplicate pairs, missingness, and mappings |
| Classification | State distinctions and all 1,024 combinations of four states across five species |
| Folds | Transitive components, determinism, disjointness, and invalid membership handling |
| Conservation | Window coordinates, byte-range behavior, missing coverage, and handle cleanup |
| External wrappers | Expected output schemas, record preservation, and skipped-record detection |
| Assembly audit/downloads | Exact orientation checks, offset diagnostics, deterministic selection, offline cache behavior, and provenance |
| Models | Known synthetic evaluation cases, undefined metrics, and explicit stub failures |
| Optional integration checks | Installed FIMO, liftOver, bedtools, and CD-HIT-EST using temporary synthetic fixtures |

```bash
make test
make lint
# Run only tests that do not require external genomics executables:
python -m pytest -m "not integration"
```

Missing executables cause their integration tests to skip; a passing suite with skips does not establish those tools work locally. The assembly audit is a separate real-input check. The starter notebook has no committed outputs and intentionally requires prepared real data before it can run.

## Remaining work and completion criteria

The detailed sequence is in the [Week 1 guide](data-README.md). The immediate dependencies are:

1. Reconcile the 3,074 canonical sequence names without guessing aliases or introducing duplicates.
2. Convert native hg19 intervals to hg38, preserving source IDs, original coordinates, failures, orientation, and sequence comparisons.
3. Establish enhancer-context handling, complete/missing measurement policies, repression sign, and repression-call definitions.
4. Acquire and verify the remaining small inputs; prepare the separately triggered genome-download workflow.
5. Implement the six-check `src/validate.py` harness and `make validate`, with expected/observed values and explicit failures.
6. Only after input validation, create real cross-species features, frozen folds, models, and candidate-ranking metrics.

Current code availability is not evidence that these steps have been run on the biological data.

## Troubleshooting

| Symptom | Meaning | Appropriate next step |
| --- | --- | --- |
| `Fill column mapping in config` | Loader schema mapping is intentionally unset | Prepare the scientifically appropriate table, then map its real headers |
| Duplicate element/cell pairs | Enhancer contexts, replicate rows, or duplicates remain | Preserve the distinction and establish the intended analysis unit |
| Expected element count or completeness failure | Input does not satisfy the declared population | Report the observed discrepancy; do not lower the expectation to obtain a pass |
| Missing assembly-audit cache | Offline replay lacks an evidence file | Confirm the GEO archive exists; use the explicit small-file fetch option when appropriate |
| Cache lacks matching provenance / has changed | Local bytes or source identity disagree with the ledger | Investigate the file; do not rewrite historical hashes to bless it |
| ZIP reader rejects `mmc3.zip` | The retained PMC response is an HTML challenge | Use the verified `Table_S2_publisher.zip`, not the failed response |
| Required executable not on PATH | Wrong environment or missing command-line tool | Activate/update the specified environment and verify the executable |
| pyBigWig lacks remote support | Its build does not support the required HTTP access | Resolve the environment build; do not download whole tracks as a fallback |
| Too few independent fold components | Chromosome/cluster constraints preclude the requested split | Report the grouping constraint; do not split linked elements to force folds |
| Local report link is missing after checkout | Reports and downloaded evidence are ignored artifacts | Acquire the documented inputs and regenerate the supported report |

## Sources and attribution

The scientific input source is Mouri et al., *Whole-genome functional characterization of RE1 silencers using a modified massively parallel reporter assay*, *Cell Genomics* 3(1):100234. [Article and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903721/), [DOI](https://doi.org/10.1016/j.xgen.2022.100234), [GEO GSE196171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171), [author reference repository](https://github.com/tewhey-lab/duoREST).

Operational references: [UCSC sequence API](https://genome.ucsc.edu/goldenPath/help/api.html), [FIMO](https://meme-suite.org/meme/doc/fimo.html), [pyBigWig](https://github.com/deeptools/pyBigWig), [bedtools getfasta](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html), and [CD-HIT](https://github.com/weizhongli/cdhit).

Use the exact source URLs, timestamps, and file hashes in the [data ledger](data/README.md) when documenting a local run. Keep study measurements, repository verification results, and future scientific conclusions separately attributable.
