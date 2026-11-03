# Week 1 — Data Acquisition and Validation Guide

**Operational plan for getting the inputs on disk and establishing what they contain.**

[Main guide](README.md) · [Research design](README%20%282%29.md) · [Authoritative catalog and provenance](data/README.md)

This file is the Week 1 guide. It is **not** a second download ledger. All actual file-level acquisition records belong in [data/README.md](data/README.md), where existing records must remain intact.

> **Current handoff:** native GEO intervals are hg19; the analysis target is hg38. Discovery and the sampled assembly audit are complete. Sequence-name reconciliation, coordinate conversion, the remaining input acquisitions, and the full six-check validation harness remain outstanding.

## Contents

1. [Stage-by-stage status](#stage-by-stage-status)
2. [What is available locally](#what-is-available-locally)
3. [Required dependency order](#required-dependency-order)
4. [Acquisition and replay procedures](#acquisition-and-replay-procedures)
5. [Preparing a defensible MPRA input](#preparing-a-defensible-mpra-input)
6. [Specification for the six-check validation harness](#specification-for-the-six-check-validation-harness)
7. [Small inputs still required](#small-inputs-still-required)
8. [Large-download target specification](#large-download-target-specification)
9. [Provenance and file-handling rules](#provenance-and-file-handling-rules)
10. [Week 1 completion checklist](#week-1-completion-checklist)

## Stage-by-stage status

| Stage | Original objective | Current state |
| --- | --- | --- |
| A — Discovery | Inspect GEO supplements before constructing downstream workflows | Complete; seven tables inspected, all schemas recorded |
| Assembly follow-up | Identify the coordinate assembly | Complete for a deterministic native-element sample: hg19 supported |
| B — Small inputs | Acquire matrix, chains, author reference repository, and supplements | Partial: Table S2 acquired for assembly verification; other requested acquisitions remain pending |
| C — Large inputs | Build a resumable, manually invoked genome-download target | Not implemented; no whole-genome downloads have run |
| D — Validation harness | Produce six expected-versus-observed checks and figures | Not implemented; `make validate` is not yet available |
| E — Provenance | Record every actual download automatically | Implemented for the small downloader; all current raw acquisitions/responses are recorded |

The source-assembly issue is resolved sufficiently to identify the required conversion step. It does not eliminate the other input and interpretation dependencies.

## What is available locally

The catalog snapshot on 2026-11-03 contains **101 downloaded files/responses totaling 33,071,763 bytes**, excluding README files and zero-byte placeholders. This includes an intentionally retained failed HTML response; a provenance record does not by itself mean a file is valid scientific input.

| Local artifact group | Role | Usability |
| --- | --- | --- |
| GEO raw archive | Six compressed source tables | Opened and inspected without extracting duplicate files |
| GEO noncanonical K562 table | Separate processed experiment | Opened; coordinate columns are entirely missing |
| GEO inventory and metadata snapshots | Discovery provenance | Retained as downloaded |
| Publisher Table S2 ZIP | Synthesized oligo sequences | Valid archive; the `WG_REST` library is used by the assembly audit |
| PMC response named `mmc3.zip` | Failed supplement fetch | HTML challenge, not a usable ZIP |
| 94 UCSC JSON responses | 47 sampled loci × two assemblies | Each response contains a bounded 202 bp sequence window |

The detailed [data catalog](data/README.md) records filenames, row counts, source schemas, and checksums. The `data/processed/` directory currently contains only its placeholder. No finished MPRA modeling table, real-data fold CSV, or scientific figure is present.

## Required dependency order

```text
Inspect immutable source files and their schemas
    ↓
Distinguish element IDs, construct IDs, contexts, and native/variant populations
    ↓
Resolve source assembly and reconcile published sequence names
    ↓
Convert native hg19 intervals to hg38 and validate sequences/orientation
    ↓
Establish activity sign, enhancer-context handling, missingness, and repression labels
    ↓
Acquire/verify motif, mapping, sequence, and annotation inputs
    ↓
Run the six-check validation harness and retain every failure
    ↓
Only then construct features, frozen folds, and scientific models
```

Some acquisitions can be independent, but a completed download does not authorize interpreting a dataset before its assembly and schema are established. In particular, hg38-to-mammal chain files cannot be applied directly to the original hg19 interval numbers.

## Acquisition and replay procedures

### 1. Verify the software environment

From the repository root:

```bash
conda activate re1-gates
python -c "import src.config"
make test
make lint
```

Use the installation instructions in [README.md](README.md#installation-and-working-environment) if the environment does not yet exist. These commands do not acquire biological data.

### 2. Restore the small GEO archive when it is absent

The assembly audit does not download its GEO archive automatically. The following is an **explicit acquisition command**, not a setup hook. It uses the existing bounded downloader and appends provenance internally:

```bash
python - <<'PY'
from src import config
from src.downloads import download_small_file

url = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE196nnn/GSE196171/suppl/GSE196171_RAW.tar"
)
destination = config.RAW_DIR / "geo" / "GSE196171" / "GSE196171_RAW.tar"
download_small_file(url, destination, fetch=True, max_bytes=100_000_000)
PY
```

The observed archive size is 31,252,480 bytes. That observation is not a permanent guarantee about a future response: the downloader enforces the supplied size cap, and the audit must still be able to read the expected archive members. If a file exists without matching provenance, the downloader raises rather than silently accepting or replacing it.

### 3. Run the audit from an existing cache

```bash
make assembly-check
```

This is the normal replay path. It requires the GEO archive, the publisher Table S2 ZIP, and the cached UCSC responses. It uses no network access and rewrites the generated assembly Markdown/CSV reports from the evidence.

### 4. Fetch missing small audit evidence explicitly

```bash
python -m src.assembly_check --fetch
```

This can acquire Table S2 and missing UCSC windows. It does not download the GEO archive, whole genomes, chain files, or bigWigs. Each sequence query requests the source 200 bp interval plus one flanking base on each side. The implementation spaces uncached UCSC requests by one second.

The output is [coordinate_assembly.md](results/tables/coordinate_assembly.md) plus [coordinate_assembly_comparisons.csv](results/tables/coordinate_assembly_comparisons.csv). An unresolved assembly comparison produces a nonzero exit after writing its evidence report.

## Preparing a defensible MPRA input

### Preserve the source populations

Use exact `project` labels when identifying native subsets:

| Label | Population | Expected distinct native elements |
| --- | --- | ---: |
| `motif` | Canonical-motif set | 8,436 |
| `noMotif` | No-canonical-motif contrast set | 4,430 |

The other project labels are not aliases for these native subsets. Preserve the original source files and labels before deriving any new table.

### Preserve enhancer context until a rule is justified

Each native element has five rows per cell table, one per enhancer context. `id` identifies a construct; `silencer` identifies the candidate element. The five contexts are En02, En09, En11, En19, and En21.

The generic loader's one-row-per-element/cell contract cannot be satisfied by simply relabeling `id`, dropping duplicates, or averaging values. A context-selection or context-modeling policy must be documented and justified first. No such rule is currently implemented.

### Reconcile sequences explicitly

Table S2 contains 24,000 rows in its `WG_REST` library. Restricting to that library avoids oligo names reused in other experiments. Exact native-name joins recover 5,362 of 8,436 canonical sequences and all 4,430 no-motif sequences.

The 3,074 unmatched canonical names need a source-backed mapping. They are not established biological absences, sequence failures, or permission to lower the expected canonical count. Preserve the unresolved names and the evidence for any eventual reconciliation.

### Convert coordinates with an audit trail

Keep original hg19 chromosome/start/end values alongside any new hg38 coordinates. Record chain identity, mapping options, orientation changes, target interval lengths, and unmapped records. Validate target sequences against the appropriate source sequence evidence; do not assume that a successful lift preserves the 200 bp span or the exact motif location.

The existing liftOver wrapper accepts a supplied chain and explicit `min_match`; it does not choose the correct assembly pair or a scientifically justified mapping threshold for the caller.

### Establish activity interpretation

The native GEO tables include `log2FoldChange` and `lfcSE`, but not a finished binary repression-call column. Some activity values are missing. The native SK-N-SH source label is `SK.N.SH`.

Prepare a documented label normalization, missingness policy, and response definition. Keep the measured source activity distinct from any transformed model target. The planned source-validation rule places repression in the left tail; do not flip values merely because a visually convenient distribution would result.

## Specification for the six-check validation harness

**This section is a specification.** It does not claim that `src/validate.py`, its figures, or `results/tables/week1_validation.md` already exist.

Every check must report its input population, expectation, observed value, decision, and any unresolved prerequisite. A missing prerequisite must not be presented as a pass. Thresholds must not be adjusted to rescue a failed assertion.

| Check | Expected evidence | Required observed output | Current readiness |
| --- | --- | --- | --- |
| 1. Population | 8,436 canonical elements; report the 4,430 contrast set if present | Distinct element counts and construct/context counts separately | Native counts verified in discovery; final prepared-table population still pending |
| 2. Length | 200 bp native elements | Full sequence-length distribution and invalid/missing counts | 9,792 direct Table S2 matches are 200 bp; full canonical reconciliation pending |
| 3. Motif position | Configured motif location, compared in a common coordinate convention | Ten genomic extractions, actual best-hit positions, orientation, and pass/fail | Target conversion, extraction inputs, and matrix still pending |
| 4. Sign | Repressing observations occupy the left tail on the source scale | Per-cell distributions, figures, label/control evidence, and explicit failure if inverted | Trusted sign/repression-call interpretation not established |
| 5. Cell coverage | All four cell types | Counts by cell and K562/SK-N-SH overlap for explicitly defined eligible measurements | All four source tables exist; usable-overlap rule and acceptance criterion pending |
| 6. m-value | FIMO score threshold 20.86 | Number at/above threshold and fraction independently called repressing | Full sequences, matrix/scoring compatibility, and repression labels pending |

### Motif-position reporting

The current configuration uses zero-based inclusive offsets `(91, 111)`; FIMO emits one-based inclusive positions. The report must show the convention used for each value and convert before comparison. The assembly audit did not run FIMO and does not validate the configured motif position.

### Sign reporting

Show the source measurement definition and the evidence identifying repressing observations. A tail shape alone cannot establish which biological class it represents. If a trusted sign check fails, fail explicitly and retain the evidence; do not silently negate a column.

### Cell-intersection reporting

Distinguish at least the presence of an element ID from availability of an eligible finite measurement in the chosen context. An intersection of source IDs is not the same as an intersection of usable measurements. The project's minimum acceptable K562/SK-N-SH overlap has not yet been specified; do not invent it.

### Expected outputs once implemented

| Artifact | Destination | Current state |
| --- | --- | --- |
| Printed and saved six-check report | `results/tables/week1_validation.md` | Not implemented |
| Repression-distribution figures | `results/figures/` | Not generated |
| Reproducible entry point | `make validate` | Not implemented |
| Supporting per-record evidence | Explicitly documented paths under `results/` | To be defined with the harness |

## Small inputs still required

| Input | Intended use | Required verification before use |
| --- | --- | --- |
| JASPAR MA0138.2 in MEME format | FIMO scoring | Matrix identifier, format, strands/background, and compatibility with the score threshold |
| Source-to-target human chain | Native hg19 → project hg38 conversion | Correct assembly direction, chain provenance, and mapping behavior |
| hg38 → rheMac10 chain | Macaque mapping | Exact assemblies and direction |
| hg38 → mm39 chain | Mouse mapping | Exact assemblies and direction |
| hg38 → canFam6 chain | Dog mapping | Exact assemblies and direction |
| hg38 → bosTau9 chain | Cow mapping | Exact assemblies and direction |
| hg38 → monDom5 chain | Opossum mapping | Exact assemblies and direction |
| Author `duoREST` reference checkout | Read-only interpretation reference | Record revision; place under a newly ignored `vendor/` directory before cloning; do not execute it |
| Remaining required supplements, especially Table S11 | Cofactor/accession interpretation | Confirm the actual sheet/schema and accession meaning |
| phyloP and phastCons URLs | Remote conservation queries | hg38 compatibility, exact track identity, and HTTP byte-range support |

Table S2 is already present because it was necessary for assembly resolution. Reading selected author files through connected tools did not create a local `vendor/` checkout. The active `.gitignore` does not yet include a `vendor/` rule; that must be added as part of any future clone workflow.

## Large-download target specification

`data-genomes` is a planned Makefile target that must be triggered manually. It must not become a dependency of `env`, `test`, `lint`, `assembly-check`, or an import-time operation.

| Requirement | Intended behavior |
| --- | --- |
| Genome inventory | Identify hg38 `.2bit` and the five target assemblies explicitly |
| Size preview | Print verified remote byte sizes and a running/total byte estimate before transfer |
| Resume | Support interrupted transfers, for example through curl resume behavior |
| Integrity | Compare against UCSC MD5 files where available; distinguish local hashes from verified publisher matches |
| Extraction tool | Resolve and document the `twoBitToFa` executable path at runtime |
| Provenance | Append a record for each actual acquired file, without rewriting earlier records |
| Activation | Download only when the user explicitly invokes the target |

Current paths in `config.py` describe expected target FASTA locations. They do not establish that the files exist, that the assemblies are available at a guessed URL, or that a transfer has happened. Remote sizes must be checked when the target is implemented; the earlier approximate total is not a verified byte budget.

## Provenance and file-handling rules

1. Keep original downloaded bytes immutable. Derived tables belong in separate processed/output paths.
2. Route supported small downloads through `download_small_file` so recording cannot be forgotten.
3. Record URL, UTC completion time, byte size, and local MD5; separately identify any publisher-checksum comparison.
4. Validate content, not just the filename or HTTP status. The retained PMC HTML challenge is a concrete example of HTTP success without a valid dataset.
5. Read archive members directly when extraction is unnecessary. If files are later extracted, record their parent archive and member path as derivation provenance rather than pretending they were separately downloaded.
6. Keep raw files and generated tables/figures ignored. Never use `git add -f` to bypass the intended exclusions.
7. Do not run downloaded reference code or introduce automatic large downloads.

Before handing off changes:

```bash
git status --short
git diff --cached --name-only
git check-ignore data/raw/geo/GSE196171/GSE196171_RAW.tar
```

README documentation may be tracked; data payloads must not be staged. The full provenance ledger is authoritative even when this guide contains a concise summary.

## Week 1 completion checklist

- [x] Enumerate the GEO supplementary directory and record actual sizes.
- [x] Open all seven observed tables and report complete schemas.
- [x] Establish that native activity and coordinate columns exist.
- [x] Resolve native source assembly by exact published-sequence comparisons.
- [x] Record current downloads and small sequence responses with provenance.
- [ ] Reconcile all canonical sequence identities.
- [ ] Convert and validate native intervals for the hg38 workflow.
- [ ] Establish enhancer-context handling, missingness, sign, and repression labels.
- [ ] Acquire and verify the remaining requested small inputs.
- [ ] Implement the resumable, separately triggered genome-download target.
- [ ] Obtain required full-genome inputs only when that target is deliberately invoked.
- [ ] Implement and run all six expected-versus-observed validation checks.
- [ ] Produce the validation Markdown report and distribution figures.
- [ ] Confirm the actual eligible K562/SK-N-SH intersection against an explicit criterion.
- [ ] Verify no data payloads are staged before the final Week 1 handoff.

The completed discovery and assembly reports remain evidence for their specific questions. They do not substitute for the unfinished six-check Week 1 report.
