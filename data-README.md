# Week 1 — Data Acquisition and Validation Guide

**Operational plan for getting the inputs on disk and establishing what they contain.**

[Main guide](README.md) · [Research design](README%20%282%29.md) · [Authoritative catalog and provenance](data/README.md)

This file is the Week 1 guide. It is **not** a second download ledger. All actual file-level acquisition records belong in [data/README.md](data/README.md), where existing records must remain intact.

> **Current handoff (2026-11-08):** native source intervals are hg19; verified downstream intervals are hg38. Sequence reconciliation and conversion are complete: **12,815 verified, 47 lifted with sequence mismatches, and 4 failed lifts out of 12,866 inputs**. Use `data/processed/elements_hg38.csv`. Remaining acquisitions and the full six-check Week 1 harness are separate work.

## Week 2 pilot setup and mouse lift

**Completed on chr8:** 390 canonical-motif elements, all 200 bp, with **10/10 exact hg38 sequence checks**. The pilot BED is `data/processed/pilot_chr8.bed`; contrast elements are excluded. See [pilot_setup.md](results/tables/pilot_setup.md).

### Dependencies are scoped to the consumer

| Step | Required genomes | Entry point |
| --- | --- | --- |
| Pilot selection, BED writing, sequence regression | hg38 | `make pilot` |
| Mouse minMatch sweep and final lift | hg38 and mm39 | `make pilot-mouse` |
| Full multi-species lift (future) | hg38 source plus rheMac10, mm39, canFam6, bosTau9, monDom5 targets | `verify_step("multi_species")`; full workflow not yet implemented |

`src/pilot.py` retains the same integrity checks: twoBitInfo, every sequence name/length, sequence count, total bases, whole-file publisher MD5, and file stability during verification. Only the dependency list changes. A required genome failure stops that step; unused genomes are not inspected. The current pilot and mouse commands do not touch, repair, or download bosTau9 or monDom5; the existing `scripts/fetch_genomes.sh` process owns those downloads.

```bash
# Offline; requires only hg38.
make pilot
# Offline; rebuilds the pilot, then requires hg38 and mm39.
make pilot-mouse
```

For a fresh checkout, `python -m src.pilot --fetch-metadata` explicitly acquires only the pilot's small hg38 metadata and the official `twoBitInfo` executable if absent. Genome transfers never run through these commands. The environment declares twoBitInfo; `data/raw/tools/twoBitInfo` is the recorded local fallback.

Pilot selection uses exact `project == "motif"`, `PILOT_CHROM = "chr8"`, and the hg38 columns in `elements_hg38.csv`. IDs remain immutable even when their names contain historical coordinates; coordinates are never reconstructed from names. The count guard remains 150–600 (roughly 300–400 expected), widths must be exactly 200 bp, and BED6 preserves strand. The regression samples ten rows from the written BED with seed `20261108`; `twoBitToFa -noMask` must reproduce the library sequence exactly in the written orientation. Only a successful regression publishes the pilot.

### Measured mouse policy

[The mouse report](results/tables/pilot_mouse_lift.md) contains the four measured thresholds (0.95, 0.5, 0.2, 0.1), rates, medians, and width deviations. **`MINMATCH_CROSS_SPECIES = 0.5`** was selected after the sweep: **227/390 (58.21%)** lift, with 163 target spans within 160–240 bp and 64 outside. The final 163 unmapped elements comprise 130 “Deleted in new” and 33 “Partially deleted in new” tool categories.

Lowering 0.5 to 0.2 gains 21 mappings but loses one as duplicated, a net increase of 20; only two additional spans fall within the width band. Lowering 0.2 to 0.1 gains five mappings, all outside the band. The choice is a mouse-pilot compromise, not a validated optimum for every species. minMatch is remapped-base coverage, not nucleotide identity. Span stability does not establish orthology or REST-motif retention.

`src.lift.lift_intervals(source_bed, source_build, target_build, chain, output_dir, min_match=...)` is the generic BED6 wrapper. It preserves both outputs and asserts exhaustive, mutually exclusive input-ID accounting. The hg19-to-hg38 conversion uses it with no minMatch override and retains its original sequence checks and 95% success guard. The cross-species sweep varies only minMatch and leaves `-multiple` disabled; it does not inherit a same-species success-rate guard or require human/mouse DNA identity.

Every sweep run is preserved under `data/processed/pilot_mouse/minmatch_<value>/`. The final directory contains `lifted.mm39.bed`, `unmapped.hg38.bed`, `liftOver.log`, and `element_audit.csv` with one row per pilot ID. All mapped intervals are retained, including the 64 width flags. Tool failure wording is recorded verbatim and is not itself evidence of lineage-specific loss. The chain download is recorded in the append-only provenance ledger.

### Motif-window coverage and ortholog FASTA

`make pilot-orthologs` runs the preceding pilot/mouse steps, then separately lifts the motif window at the **same minMatch 0.5** and extracts sequences only for jointly mapped elements. Genome dependencies remain hg38/mm39; no other species genome is inspected or downloaded. [The complete report](results/tables/pilot_mouse_orthologs.md) includes the full length distribution.

The current result is **217 MAPPED + 15 UNDETERMINED + 158 ABSENT = 390**. Exactly 217 records are written to `data/processed/pilot_mouse/orthologs.mm39.fa`; `orthologs.mm39.csv` beside it contains all 390 inputs. Of the MAPPED elements, **60** have unpadded widths outside 160–240 bp. Extracted lengths range from **136 to 3,285 bp**; **0** are shorter than 50 bp, **0** contain N, and **0** require clipping at chromosome boundaries.

| Element lifted | Motif window lifted | Call | Pilot count | Downstream rule |
| --- | --- | --- | ---: | --- |
| Yes | Yes | MAPPED | 217 | Include padded element sequence for future rescoring |
| Yes | No | UNDETERMINED | 10 | Exclude mouse; cannot distinguish motif absence from missing alignment |
| No | Yes | UNDETERMINED | 5 | User-approved motif-only policy; preserve motif coordinates but exclude mouse |
| No | No | ABSENT | 158 | Alignment-based candidate loss; not a proven biological absence |

The motif-only policy was explicitly approved after observing the fourth logical outcome. Both UNDETERMINED subtypes retain their separate `call_reason`; neither can enter the FASTA or silently become ABSENT. `exclude_mouse_downstream=True` and `include_in_fasta=False` make that rule explicit in the sidecar.

Coordinate convention is inherited unchanged from config: `MOTIF_SPAN=(91,111)` is **zero-based, both endpoints inclusive**. For an element `[s,e)`, the 21 bp motif BED is `[s+91,s+112)` on plus strand and `[e-112,e-91)` on minus strand. These are configured offsets, not a new FIMO validation of motif placement. A passing minMatch 0.5 window need not align all 21 source bases; MAPPED records mapping availability, not motif retention.

Both successful tracks must agree on chromosome and strand, with the mapped motif window contained in its mapped element. All 217 passed. Extraction uses the **lifted element interval**, padded 15 bp on each genomic side, with BED6 strand and `twoBitToFa -noMask`. Chromosome-end padding would be clipped explicitly and recorded per side, never concealed. Each extracted length is asserted equal to padded_end − padded_start. Sequences under 50 bp would be flagged and retained to keep FASTA count equal to MAPPED count.

The sidecar records both hg38 source tracks, both mm39 mapped tracks, independent failure reasons, the unpadded width-band flag, padded extraction coordinates, actual padding, sequence length, short-sequence flag, and N-base count. Missing sequence measurements remain empty for excluded/absent elements. The original 21 bp motif BED, both motif lift outputs, and command log are retained under `pilot_mouse/motif_window/`; extraction diagnostics and an element-track replay are under `pilot_mouse/ortholog_extraction/`. The replay must match the previous final element lift exactly before joining tracks.

These three calls are mapping-availability states. They do not replace the later motif-score/evolutionary classification, and UNDETERMINED must not be encoded as a negative motif score. Outputs and printed reports are deterministic; replay evidence is in [pilot_mouse_orthologs_reproducibility.json](results/tables/pilot_mouse_orthologs_reproducibility.json).

## Contents

1. [Stage-by-stage status](#stage-by-stage-status)
2. [What is available locally](#what-is-available-locally)
3. [Required dependency order](#required-dependency-order)
4. [Acquisition and replay procedures](#acquisition-and-replay-procedures)
5. [Preparing a defensible MPRA input](#preparing-a-defensible-mpra-input)
6. [Verified hg38 handoff](#verified-hg38-handoff)
7. [Specification for the six-check validation harness](#specification-for-the-six-check-validation-harness)
8. [Small inputs still required](#small-inputs-still-required)
9. [Large-download target specification](#large-download-target-specification)
10. [Provenance and file-handling rules](#provenance-and-file-handling-rules)
11. [Week 1 completion checklist](#week-1-completion-checklist)

## Stage-by-stage status

| Stage | Original objective | Current state |
| --- | --- | --- |
| A — Discovery | Inspect GEO supplements before constructing downstream workflows | Complete; seven tables inspected, all schemas recorded |
| Assembly follow-up | Identify the coordinate assembly | Complete for a deterministic native-element sample: hg19 supported |
| B — Small inputs | Acquire matrix, chains, author reference repository, and supplements | Partial: Table S2 acquired for assembly verification; other requested acquisitions remain pending |
| C — Large inputs | Build a resumable, manually invoked genome-download target | `data-genomes` remains pending; hg38 was explicitly authorized, downloaded, and checksum-verified for conversion |
| D — Validation harness | Produce six expected-versus-observed checks and figures | Not implemented; `make validate` is not yet available |
| E — Provenance | Record every actual download automatically | Small and conversion-input downloaders append provenance internally; all five new conversion references are recorded |

Source-assembly identification and the sequence-verified conversion are complete. Other input and activity-interpretation dependencies remain.

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

Resolved on 2026-11-08: the remaining 3,074 names occur as explicit tokens in parenthesized, semicolon-separated Table S2 oligo names. `src/lift.py` expands that exact syntax, requires one source row per element, and preserves the full publisher name and join method. All 12,866 native oligos now have 200 bp library sequences; no reference DNA was imputed.

### Convert coordinates with an audit trail

Keep original hg19 chromosome/start/end values alongside any new hg38 coordinates. Record chain identity, mapping options, orientation changes, target interval lengths, and unmapped records. Validate target sequences against the appropriate source sequence evidence; do not assume that a successful lift preserves the 200 bp span or the exact motif location.

The existing liftOver wrapper accepts a supplied chain and explicit `min_match`; it does not choose the correct assembly pair or a scientifically justified mapping threshold for the caller.

### Establish activity interpretation

The native GEO tables include `log2FoldChange` and `lfcSE`, but not a finished binary repression-call column. Some activity values are missing. The native SK-N-SH source label is `SK.N.SH`.

Prepare a documented label normalization, missingness policy, and response definition. Keep the measured source activity distinct from any transformed model target. The planned source-validation rule places repression in the left tail; do not flip values merely because a visually convenient distribution would result.

## Verified hg38 handoff

**Completed 2026-11-08:** 12,815 verified, 47 lifted with sequence mismatches, and 4 failed lifts out of 12,866 inputs. The verified CSV contains **8,400 canonical** and **4,415 no-motif** elements. Original input expectations remain 8,436 and 4,430, respectively; losses are recorded rather than hidden by changing those expectations.

### Reproduce the conversion

Activate the project environment so `liftOver` and `twoBitToFa` are on PATH. From the repository root:

```bash
# Offline: inputs must already be present.
make lift

# Explicit acquisition only when references are missing. This downloads hg38.
python -m src.fetch_lift_inputs --genome
```

The downloader prints a cumulative byte budget, resumes interrupted transfers, and appends each completed download to `data/README.md`. hg38 must match UCSC's published MD5. Chain and RepeatMasker MD5s document local identity; no publisher checksum is claimed for them. Source GEO and Table S2 archives must already exist. The other species genomes are not fetched by this command.

`src/lift.py` accepts an optional `--input` CSV with `element_id,chrom,start,end,sequence,source_build` and optional `strand,project`. `source_build` must explicitly be `hg19`; missing/invalid sequences, duplicate IDs, fractional coordinates, or non-200-bp intervals fail before BED writing. This prevents accidentally lifting an already-hg38 table a second time.

### Outputs and coordinate contract

| Artifact | Purpose |
| --- | --- |
| `data/processed/elements_hg38.csv` | Only exact sequence-verified elements; downstream coordinate source |
| `data/processed/liftover_hg19_to_hg38/input_elements.csv` | Full original population, published oligos, explicit alias evidence, adapter diagnostics |
| `input.hg19.bed` in that directory | Original 200 bp BED6 records, one per element |
| `lifted.hg38.bed` | Unmodified mapped output, including intervals that subsequently fail sequence verification |
| `unmapped.hg19.bed` | Unmodified unmapped output, including liftOver's per-record reasons |
| `lifted.hg38.fa` | `twoBitToFa -bed=... -noMask` output for every mapped record |
| `element_audit.csv` | Every input exactly once; bucket, source/target coordinates, library/extracted DNA, mismatch details, repeat overlap |
| `liftOver.log`, `twoBitToFa.log` | Actual argv, stdout, stderr |
| `results/tables/liftover_hg19_to_hg38.md` | Counts, checks, mismatch breakdown, repeat comparison, input/tool MD5s |
| `results/tables/liftover_hg19_to_hg38_reproducibility.json` | Byte identity check across two real-data runs |

All intervals use zero-based, half-open coordinates. The verified CSV's `chrom,start,end,strand` are **hg38**; `sequence` is the uppercase library genomic core, with `strand` adjusted to describe its orientation relative to hg38. `hg19_chrom,hg19_start,hg19_end,hg19_strand` preserve the input. When no strand was published, `strand_supplied=False` records that the BED's original plus strand was only a mapping convention. `source_oligo_name` and `sequence_join` preserve the library evidence.

### Strict verification and exclusion rules

- liftOver uses its default `-minMatch=0.95`, with no override and no multiple-mapping option. An observed mapping success below 95% stops output publication; this operational guard is separate from the alignment parameter.
- Library and extracted DNA are uppercased, and both orientations are tested. Only an exact unambiguous A/C/G/T match of all 200 bases passes.
- Current published oligos are all 200 bp: **no adapters were stripped**. For future adapter-bearing inputs, explicit `--adapter-5p`/`--adapter-3p` settings or one unambiguous shared-flank split can remove adapters. Automatic inference requires at least two distinct long oligos, identical flanks of at least eight bases, and exactly one split leaving 200 bp. It never uses reference DNA to select a favorable subwindow. Ambiguous long oligos fail and remain counted.
- Mismatch diagnostics distinguish changed lengths, one-base differences, 2–49 differences, >=50/200 differences, ambiguous bases, and unresolved adapters. The high-divergence label is descriptive, not a homology inference. No mismatch is rescued by shifting, trimming the genome interval, or relaxing equality.
- Bucket membership must be exhaustive and mutually exclusive. A failed rerun removes the previous consumable CSV and writes a failure report, preventing stale verified data from appearing current.

### Repeat-bias interpretation

The audit intersects **original hg19 intervals** with the full hg19 RepeatMasker track. Any overlap of at least one base counts, on either strand, including simple and low-complexity repeats. Each element contributes once to an inside/outside denominator. The report separates mapping failure from total exclusion, and reports both populations separately as well as pooled.

The descriptive rule for a large absolute failure-rate difference is **at least five percentage points**, fixed before examining outcomes. The report also gives the risk ratio. A small absolute difference cannot establish that evolutionary-class bias is absent: repeat overlap is a proxy, and sequence-mismatch exclusions also matter.

### Downstream use

Use the verified hg38 CSV for conservation, ENCODE overlaps, and hg38-to-species chains. Keep raw hg19 artifacts only for traceability and offline reproduction of this one conversion. Activity/context handling, sign conventions, repression calls, and motif-position validation remain separate prerequisites; the coordinate table does not resolve them.

## Specification for the six-check validation harness

**This section is a specification.** It does not claim that `src/validate.py`, its figures, or `results/tables/week1_validation.md` already exist.

Every check must report its input population, expectation, observed value, decision, and any unresolved prerequisite. A missing prerequisite must not be presented as a pass. Thresholds must not be adjusted to rescue a failed assertion.

| Check | Expected evidence | Required observed output | Current readiness |
| --- | --- | --- | --- |
| 1. Population | 8,436 canonical elements; report the 4,430 contrast set if present | Distinct element counts and construct/context counts separately | Native counts verified in discovery; final prepared-table population still pending |
| 2. Length | 200 bp native elements | Full sequence-length distribution and invalid/missing counts | All 12,866 native library sequences are 200 bp; lifted length changes and mismatches are retained in the conversion audit |
| 3. Motif position | Configured motif location, compared in a common coordinate convention | Ten genomic extractions, actual best-hit positions, orientation, and pass/fail | Target conversion and exact-sequence extraction complete for retained elements; matrix/FIMO motif-position validation pending |
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

The human hg19-to-hg38 chain is already downloaded and validated through the completed conversion. The remaining comparative chains are separate inputs.

| Input | Intended use | Required verification before use |
| --- | --- | --- |
| JASPAR MA0138.2 in MEME format | FIMO scoring | Matrix identifier, format, strands/background, and compatibility with the score threshold |
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
- [x] Reconcile all canonical sequence identities through explicit publisher aliases.
- [x] Convert and validate native intervals for the hg38 workflow; retain every failure/mismatch.
- [ ] Establish enhancer-context handling, missingness, sign, and repression labels.
- [ ] Acquire and verify the remaining requested small inputs.
- [ ] Implement the resumable, separately triggered genome-download target.
- [x] Obtain and checksum-verify hg38 after explicit user authorization.
- [ ] Complete the remaining species genome acquisitions with explicit invocation and provenance.
- [ ] Implement and run all six expected-versus-observed validation checks.
- [ ] Produce the validation Markdown report and distribution figures.
- [ ] Confirm the actual eligible K562/SK-N-SH intersection against an explicit criterion.
- [ ] Verify no data payloads are staged before the final Week 1 handoff.

The completed discovery and assembly reports remain evidence for their specific questions. They do not substitute for the unfinished six-check Week 1 report.
