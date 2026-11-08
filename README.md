# RE1 Silencer Gates

**Comparative genomics for prioritizing human RE1 silencers as candidate transcriptional NOT-gate components.**

Gabrielle DeCuir and Vrinda Sehgal · Johns Hopkins Applied Comparative Genomics

## What this is

This repository connects published reporter measurements, verified human sequences, and comparative genome mappings to ask which RE1 elements merit experimental follow-up. RE1 is a binding site for the transcriptional repressor REST. The source experiment, MPRAduo, measures regulatory activity in an enhancer–silencer reporter context. The engineering motivation is to identify candidates whose repression could be useful in transcriptional circuits. Reporter activity alone does not demonstrate a working gate. [Source study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903721/)

**The completed workflow reaches sequence extraction for the mouse pilot.** It has resolved the source assembly, produced **12,815 sequence-verified hg38 elements**, selected **390 canonical chr8 elements**, and extracted **217 mouse sequences** with motif-window mapping tracked separately. These are data-validation and mapping results; conservation effects, predictive performance, and gate behavior have not been established. [Human conversion](results/tables/liftover_hg19_to_hg38.md) · [Pilot](results/tables/pilot_setup.md) · [Mouse ortholog sequences](results/tables/pilot_mouse_orthologs.md)

> **Coordinate contract:** published native coordinates are hg19. Downstream genomic work uses `data/processed/elements_hg38.csv`, whose retained sequences were verified against hg38. Original hg19 coordinates remain as traceability fields. Element IDs retain their original names and may contain hg19 coordinate text; never reconstruct working coordinates from an ID.

### Contents

- [The question](#the-question)
- [Approach and pipeline](#approach-and-pipeline)
- [Results](#results)
- [Design decisions and why](#design-decisions-and-why)
- [Verification](#verification)
- [Reproducing](#reproducing)
- [Repository layout](#repository-layout)
- [Data sources](#data-sources)
- [References](#references)
- [Scope and known limitations](#scope-and-known-limitations)
- [Team](#team)

## The question

**Do evolutionary retention and sequence conservation help predict human RE1-mediated repression beyond the information in REST motif score?**

The proposed comparison is a motif-only baseline against models that also use conservation, cross-species motif evidence, or cofactor features. Evaluation should hold out related genomic and sequence groups together. Any resulting ranking would prioritize candidates for experiments, rather than establish their circuit performance or the causal effect of evolutionary age.

The source cell types are **GM12878, HepG2, K562, and SK-N-SH**. The native populations are canonical-motif elements (`project == "motif"`) and REST-peak elements without the canonical motif (`project == "noMotif"`). The contrast population is retained in the verified human dataset but excluded from the canonical-element pilot. [Population and cell-table evidence](results/tables/coordinate_assembly.md) · [Pilot selection](results/tables/pilot_setup.md)

## Approach and pipeline

### Observation units and coordinate conventions

The GEO archive contains processed construct-level activity and genomic intervals, despite its `RAW.tar` filename. A construct row is not an independent genomic element: enhancer contexts and cell types can share an element. In the native WGREST tables, `silencer` identifies the element, while `id` identifies the construct; `log2FoldChange` and `lfcSE` are published activity fields. Collapsing identical coordinate metadata for sequence work does **not** aggregate activity. The full schema inventory is available in the [local report directory](results/tables/) and indexed by the [data catalog](data/README.md).

The coordinate workflow uses zero-based, half-open BED intervals. Human elements are **200 bp**. `MOTIF_SPAN=(91, 111)` uses zero-based, **inclusive** endpoints, giving a **21 bp** window. For an element `[s,e)`, the motif BED is `[s+91,s+112)` on the plus strand and `[e-112,e-91)` on the minus strand. This is the implemented coordinate convention, not an independent FIMO confirmation of the motif's location. [Conversion contract](results/tables/liftover_hg19_to_hg38.md) · [Motif-window contract](results/tables/pilot_mouse_orthologs.md)

### Implemented stages, in dependency order

Paths below are relative to the repository root. Each automated data target prints its report and writes the linked Markdown artifact.

| Stage | Input | Method and parameters | Principal outputs | Make target | Evidence |
| --- | --- | --- | --- | --- | --- |
| Inspect published inputs | GEO supplementary archive and processed tables | Inspect actual archive members, schemas, activity fields, coordinates, and native populations | Downloaded GEO inputs; schema and row-count inventory | No dedicated target; inspection is recorded, not an automated discovery command | GEO inventory in [report directory](results/tables/), indexed in [data catalog](data/README.md) |
| Resolve source assembly | Native GEO intervals, Table S2 oligos, cached UCSC reference queries | Deterministically sample before reference comparison; exact full-length equality in either orientation at unchanged coordinates | `results/tables/coordinate_assembly_comparisons.csv` | `make assembly-check` | [Assembly evidence](results/tables/coordinate_assembly.md) |
| Convert and verify human elements | Reconciled library sequences, hg19 intervals, hg19→hg38 chain, hg38 reference, hg19 RepeatMasker | Default `minMatch=0.95`; exact sequence verification; retain every failure and compare repeat-associated attrition | `data/processed/elements_hg38.csv`; BEDs, FASTA, and complete audit in `data/processed/liftover_hg19_to_hg38/` | `make lift` | [Conversion report](results/tables/liftover_hg19_to_hg38.md) |
| Select pilot | Verified hg38 CSV and hg38 reference | Verify required genome; select canonical chr8 elements; assert interval widths; extract a fixed random sample for exact sequence regression | `data/processed/pilot_chr8.bed`; evidence in `data/processed/pilot_setup/` | `make pilot` | [Pilot report](results/tables/pilot_setup.md) |
| Evaluate mouse mapping | Pilot BED, hg38→mm39 chain, hg38/mm39 references | Sweep `minMatch` over 0.95, 0.5, 0.2, 0.1; retain failures and width flags; repeat configured value for final output | `data/processed/pilot_mouse/sweep.csv`; per-threshold directories; `final/lifted.mm39.bed`, `final/unmapped.hg38.bed`, `final/element_audit.csv` | `make pilot-mouse` | [Mouse lift report](results/tables/pilot_mouse_lift.md) |
| Resolve motif coverage and extract mouse DNA | Pilot, final element-track outputs, same chain and references | Lift motif windows independently at 0.5; join by ID; validate concordance; extract MAPPED element spans with 15 bp padding per side | `data/processed/pilot_mouse/orthologs.mm39.fa`, `orthologs.mm39.csv`; both motif lift outputs and extraction evidence | `make pilot-orthologs` | [Ortholog sequence report](results/tables/pilot_mouse_orthologs.md) |

`make pilot-orthologs` invokes `make pilot-mouse`, which invokes `make pilot`. The human conversion is an explicit upstream handoff; downstream targets do not relift hg19. All these Make data targets are offline. Acquisition is separate.

### Genome dependencies

A gate verifies only the references used by that step, with the same integrity checks regardless of scope. Missing or changing unrelated downloads do not block the pilot.

| Step | Required references | Operational status |
| --- | --- | --- |
| Pilot selection, BED writing, human regression | hg38 | Implemented and run |
| Mouse sweep, motif-window mapping, mouse extraction | hg38 and mm39 | Implemented and run |
| Full comparative workflow | hg38 plus rheMac10, mm39, canFam6, bosTau9, monDom5 | Dependency policy defined; full analysis not run |

See [pilot dependency policy](results/tables/pilot_setup.md) and [`verify_step`](src/pilot.py). Verification never downloads or repairs a genome. Independent genome acquisition remains owned by the manual download script.

## Results

### Source assembly: hg19, not hg38

Published oligos from Table S2 were compared with reference DNA extracted at the **unchanged published interval numbers**. Sampling selected a directly name-matched element from each available chromosome/population group, ordered by start and ID before querying the reference. Both orientations were tested by exact equality; no alignment tolerance or motif score entered the decision.

| Candidate assembly | Sampled oligos | Exact matches at published coordinates |
| --- | ---: | ---: |
| hg19 | 47 | 47 |
| hg38 | 47 | 10 |

This supports hg19 for the sampled native populations. It is not an assembly audit of every variant or library in the publication. Diagnostic flanking-coordinate comparisons are preserved in the [assembly report](results/tables/coordinate_assembly.md) and [individual comparison CSV](results/tables/coordinate_assembly_comparisons.csv).

The project converts to hg38 because its intended mammalian conservation tracks, ENCODE integration, and target-species chain workflow use hg38 as the common reference. Keeping the activity table in hg19 would require repeated coordinate reconciliation at later joins. The conversion is therefore performed once, with original coordinates preserved and a sequence-verified hg38 table as the downstream input. Exact conservation tracks and ENCODE accessions still need to be selected and validated.

### Human conversion and retention

Sequence reconciliation recovered **9,792 exact-name matches** and **3,074 explicit aliases** from Table S2. Aliases are published tokens in oligo names, not names inferred from nearby coordinates. All **12,866** source intervals and library cores were **200 bp**; no adapters needed stripping. [Conversion report](results/tables/liftover_hg19_to_hg38.md)

The conversion keeps liftOver's same-species default, `minMatch=0.95`, and does not enable `-multiple`. It retains mapped and unmapped outputs, extracts hg38 DNA using `twoBitToFa -noMask`, uppercases library DNA, and compares both orientations by full-length exact equality. Only verified sequences enter the downstream CSV; the audit preserves mismatches, failed lifts, source coordinates, and reasons.

| Population | Input | Lifted + verified | Lifted + mismatch | Failed to lift |
| --- | ---: | ---: | ---: | ---: |
| All native elements | 12,866 | 12,815 | 47 | 4 |
| Canonical motif | 8,436 | 8,400 | 32 | 4 |
| REST-peak/no-canonical-motif contrast | 4,430 | 4,415 | 15 | 0 |

**Conversion attrition is 0.40%**, calculated as `(47 + 4) / 12,866 × 100`, rounded from the report's exclusion counts. This is loss from the verified dataset, not evidence of evolutionary loss. The bucket sums are asserted for the complete input and each population. [Exhaustive accounting](results/tables/liftover_hg19_to_hg38.md)

| Sequence mismatch category | Elements |
| --- | ---: |
| Length differs | 13 |
| One-base difference | 28 |
| Multiple-base differences | 6 |
| High divergence/unrelated, ambiguous bases, or unresolved adapters | 0 in each category |

The failed lifts comprise **2 “Deleted in new”** and **2 “Partially deleted in new”** records. No interval was resized to make it pass. [Mismatch and failure details](results/tables/liftover_hg19_to_hg38.md)

Repeat-associated failure was checked against the original hg19 RepeatMasker track, using any overlapping repeat annotation and counting each element once.

| Original repeat overlap | Input elements | Failed lifts | Lift failure rate | Total excluded, including mismatches |
| --- | ---: | ---: | ---: | ---: |
| Inside repeats | 4,630 | 2 | 0.0432% | 25 |
| Outside repeats | 8,236 | 2 | 0.0243% | 26 |

The absolute failure-rate difference is **+0.0189 percentage points**, below the report's descriptive **5-percentage-point** threshold for a large difference. This does not establish absence of evolutionary-class bias: repeat overlap is only a proxy, and sequence-mismatch exclusions also matter. [Repeat-bias analysis and population breakdown](results/tables/liftover_hg19_to_hg38.md)

### Pilot selection and reference verification

A single chromosome keeps the full path from coordinates to comparative evidence small enough to inspect before scaling. **chr8** was specified as the pilot chromosome and recorded in `PILOT_CHROM`; it was not selected after comparing mapping outcomes. Its **390 canonical-motif elements** meet the expected **300–400** range and the operational **150–600** guard. This is an engineering pilot, not a demonstrated representative sample of the genome. Contrast elements are excluded. [Pilot report](results/tables/pilot_setup.md)

All pilot intervals are **200 bp**. A fixed random sample extracted from the written BED gives **10/10 exact oriented matches** to the verified hg38 CSV. The report lists the sampled IDs and coordinates. [Regression evidence](results/tables/pilot_setup.md)

| Required assembly | Sequences observed / UCSC | Bases observed / UCSC | Whole-file publisher MD5 | Result |
| --- | ---: | ---: | --- | --- |
| hg38 | 455 / 455 | 3,209,286,105 / 3,209,286,105 | Match | PASS |
| mm39 | 61 / 61 | 2,728,222,451 / 2,728,222,451 | Match | PASS |

Counts include alternate and unplaced sequences, not only primary chromosomes. Exact digest values and per-reference evidence are in the [mouse genome gate](results/tables/pilot_mouse_lift.md). This table makes no verification claim about references unused by the mouse workflow.

### Mouse minMatch sweep

**minMatch is the fraction of source bases successfully remapped, not nucleotide identity.** It cannot establish sequence conservation or preservation of the REST motif. [UCSC parameter explanation](https://groups.google.com/a/soe.ucsc.edu/g/genome/c/YbyOe7AEXj4)

Only minMatch varied in the pilot sweep; `-multiple` remained disabled. The width band describes the unpadded target element. Exact band endpoints are included, so a width of **160 or 240 bp** is not counted as exceeding **20%** deviation from **200 bp**. [Sweep report](results/tables/pilot_mouse_lift.md)

| minMatch | Lifted / input | Lift rate | Median target width (bp) | Within 160–240 bp | Outside 160–240 bp |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.95 | 63 / 390 | 16.15% | 205 | 56 | 7 |
| **0.5** | **227 / 390** | **58.21%** | **189** | **163** | **64** |
| 0.2 | 247 / 390 | 63.33% | 185 | 165 | 82 |
| 0.1 | 252 / 390 | 64.62% | 184 | 165 | 87 |

**`MINMATCH_CROSS_SPECIES = 0.5` is the mouse-pilot policy.** Moving to 0.2 adds only **2** net in-band mappings while the out-of-band count increases from **64 to 82**. The mapped sets are not strictly nested: this transition gains **21 IDs** and loses **1** as duplicated. Moving to 0.1 adds **5** mappings, all out of band. The selected value captures most width-stable coverage without those additional distorted spans. Width is a diagnostic, not a proof of mapping precision; all successful mappings remain in the audit. This value has not been validated as optimal for other species. [Measured tradeoff and rationale](results/tables/pilot_mouse_lift.md)

The final element lift retains **227/390 (58.21%)** and records **163** failures, whose verbatim categories are:

| Unmapped reason | Elements |
| --- | ---: |
| Deleted in new | 130 |
| Partially deleted in new | 33 |

Both outputs and an audit row for every input are preserved. [Final mouse lift](results/tables/pilot_mouse_lift.md)

### Motif-window coverage and the three-way call

At `minMatch=0.5`, an element can pass even if its motif-bearing region does not map. Scanning such a target interval and treating “no hit” as motif loss would mix alignment gaps with biological absence. An independent lift of the **21 bp** motif window, at the same **0.5** threshold, tracks whether the relevant region is available for scoring. [Motif-window report](results/tables/pilot_mouse_orthologs.md)

| Element track | Motif-window track | Call | Observed elements | Interpretation |
| --- | --- | --- | ---: | --- |
| Lifted | Lifted | MAPPED | 217 | Eligible for sequence extraction after chromosome, strand, and containment checks |
| Lifted | Unmapped | UNDETERMINED | 10 | Element mapped without a mapped motif window; exclude mouse downstream |
| Unmapped | Lifted | UNDETERMINED | 5 | Motif-only mapping; retain its coordinates, but exclude mouse downstream |
| Unmapped | Unmapped | ABSENT | 158 | Candidate alignment loss; not proven biological loss |

The totals are **217 MAPPED + 15 UNDETERMINED + 158 ABSENT = 390**, asserted in code. Motif-only outcomes deliberately remain UNDETERMINED. All **217/217** jointly lifted records pass chromosome/strand agreement and motif-within-element containment. A discordant joint lift is a hard failure, not a silently accepted ortholog. [Call accounting](results/tables/pilot_mouse_orthologs.md)

A passing motif-window lift can still map only part of the window. **MAPPED establishes mapping availability; it does not establish an intact or functional motif.** No motif rescoring has yet been reported.

### Extracted mouse sequences

Only MAPPED elements enter `orthologs.mm39.fa`. Extraction uses the lifted **element** interval, with **15 bp** padding on each side, BED strand, and `twoBitToFa -noMask`. Width flags describe the original unpadded interval. The sidecar keeps every pilot element, including unmapped reasons and partially available coordinates. [Extraction report and full length-frequency table](results/tables/pilot_mouse_orthologs.md)

| Check | Observed |
| --- | ---: |
| FASTA records, equal to MAPPED count | 217 |
| Sidecar rows, equal to pilot input | 390 |
| MAPPED elements within 160–240 bp before padding | 157 |
| MAPPED elements outside 160–240 bp before padding | 60 |
| Extracted sequence lengths, minimum–maximum from the reported distribution | 136–3,285 bp |
| Extracted sequences under 50 bp | 0 |
| Sequences containing N | 0 |
| Elements with padding clipped at a chromosome boundary | 0 |

Long or short target spans are retained and flagged rather than removed to improve apparent coverage. These sequences are ready for a scoring stage; they are not yet evidence of retained or lost REST binding.

### Completion boundary

| Component | Completed work | Remaining work |
| --- | --- | --- |
| Inputs and assemblies | GEO inspection; sampled source-assembly discrimination; exhaustive human conversion and sequence verification | Resolve activity aggregation, missing measurements, repression sign, and repression-call definitions |
| Mouse pilot | Required-reference checks, canonical BED, sequence regression, minMatch sweep, dual-track calls, padded sequence extraction | Score motifs and assess how mapping uncertainty and width distortion affect interpretation |
| Comparative features | Generic source/target lift wrapper; conservation and motif-scoring helper code | Validate production motif inputs and score scale, select conservation tracks, acquire/validate cofactor accessions, run remaining species |
| Evolutionary classification | Helper logic for supplied species states | Connect measured motif evidence and missingness to biological classifications; no production classifications reported |
| Predictive evaluation | Chromosome/sequence-cluster fold builder and evaluator for supplied predictions | Produce real fold assignments, implement training, fit held-out models, and estimate performance |
| Candidate prioritization | Research design and explicit metric interfaces | Define and implement gate-quality metrics, rank candidates, and test experimentally |

Model-training functions in [`src/models.py`](src/models.py) and gate-quality functions in [`src/metrics.py`](src/metrics.py) intentionally raise `NotImplementedError`. The proposed consolidated activity/FIMO validation harness, `src/validate.py` and `make validate`, does not exist. A configuration constant or tested helper is not a completed biological analysis.

## Design decisions and why

| Decision | Reason | Boundary |
| --- | --- | --- |
| Convert once to hg38 | Establish a common coordinate system for conservation, ENCODE features, and species chains; keep source coordinates for audit | Downstream code reads verified hg38 columns, never coordinate text embedded in IDs |
| Use liftOver chains rather than parse MAF blocks | The immediate task is interval projection with explicit failure accounting. Existing chains support a small, inspectable workflow without building a multi-alignment parser | A lifted span is not a base-by-base multi-species alignment or proof of unique biological orthology |
| Keep reference genomes as `.2bit` | `twoBitToFa` extracts requested intervals directly; whole-genome FASTA expansion is unnecessary | Both file integrity and extracted sequences still require validation |
| Read conservation bigWigs remotely | The implemented helper requests genomic windows through pyBigWig and checks HTTP byte-range support; avoids whole-track downloads | Production track URLs remain unset; missing coverage is reported, not imputed as zero |
| Lift element and motif window separately | Whole-element coverage can conceal failure to align the motif region | Partial-window alignment remains possible; UNDETERMINED is excluded for that species |
| Keep width distortions and unmapped records | Attrition and alignment uncertainty must remain visible rather than disappear through filtering | Width stability is only a diagnostic; no motif-retention claim follows from it |
| Gate references by dependency | A step needs valid references it actually reads | Every required reference still faces the full integrity gate |
| Group folds by chromosome and sequence similarity | Prevent local genomic context and clustered related sequences from appearing on both sides of evaluation | Chromosome/cluster connected components are kept intact; too few components cause failure rather than forced splitting |

The fold implementation links chromosomes joined by a CD-HIT-EST cluster transitively, then assigns whole components deterministically. It checks that neither a chromosome nor a cluster crosses folds. Greedy clustering cannot guarantee detection of every related sequence, and this single-chromosome pilot is not the dataset for the planned chromosome-held-out model comparison. No real-data folds or model results have been produced. [Fold implementation](src/folds.py) · [Conservation reader](src/conservation.py)

## Verification

### Integrity and accounting checks

The implemented checks below are distinct from the unfinished activity-validation work. Failures remain visible; thresholds and interval lengths are not changed to obtain a pass.

| Layer | Checks performed | Evidence or implementation |
| --- | --- | --- |
| Downloads and caches | Explicit HTTPS downloads; size limits or expected transfer size; resumable reference transfers; cache URL/digest agreement; provenance appended with source URL, UTC time, size, and local MD5 | [`downloads.py`](src/downloads.py), [`fetch_lift_inputs.py`](src/fetch_lift_inputs.py), [ledger](data/README.md) |
| File format and source identity | Validate actual archive contents rather than trust extensions; restrict oligos to the intended library; unique native IDs; consistent coordinates across cell tables; exact-name or explicit-alias joins only | [Assembly report](results/tables/coordinate_assembly.md), [human conversion](results/tables/liftover_hg19_to_hg38.md) |
| Required genomes | Require nonzero, stable files; successful `twoBitInfo`; exact per-contig equality to UCSC `.chrom.sizes`; sequence counts and total bases; whole-file MD5 equal to publisher manifest; detect file changes during verification | [Genome gate](results/tables/pilot_mouse_lift.md), [`pilot.py`](src/pilot.py) |
| Source coordinates | Check declared source/working build, integer coordinates, valid strand, unique BED-safe IDs, and expected widths; use the original published inputs for human conversion and the verified CSV for pilot selection | [`lift.py`](src/lift.py), [`pilot.py`](src/pilot.py) |
| Human lift | Mapped and unmapped IDs partition the source; retain unchanged failed records and their reasons; enforce the same-species success guard; exact sequence verification in either orientation; no ambiguous-base passes | [Conversion report](results/tables/liftover_hg19_to_hg38.md) |
| Human publication | Verified, mismatched, and failed buckets sum to input overall and by population; only verified IDs exported; source coordinates and exclusions retained; repeat-overlap denominators include failures | [Conversion accounting and bias check](results/tables/liftover_hg19_to_hg38.md) |
| Pilot | Canonical-only selection, count guard, exact widths, valid hg38 sequences; extract the seeded sample from the written BED and require oriented sequence equality | [Pilot regression](results/tables/pilot_setup.md) |
| Mouse sweep | Pilot BED agrees with verified CSV; target intervals lie within the reference; every threshold accounts for every ID; final mapping and failure records equal the selected sweep entry | [Mouse lift](results/tables/pilot_mouse_lift.md) |
| Dual-track join | Motif windows have expected widths and IDs; each track partitions input; replayed element lift equals stored final output; jointly mapped intervals agree in chromosome/strand and containment | [Motif-window report](results/tables/pilot_mouse_orthologs.md) |
| Calls and extraction | Calls exhaustively partition input; sidecar has unique complete IDs; FASTA IDs equal exactly MAPPED IDs; valid uppercase DNA; lengths equal padded coordinate spans; boundary clipping recorded; short sequences and Ns counted | [Extraction report](results/tables/pilot_mouse_orthologs.md) |
| Software contracts | Unit tests cover malformed inputs, external-wrapper outputs, missingness, classification, fold separation, and evaluator behavior; integration tests exercise installed external tools on synthetic fixtures | [`tests/`](tests/), `make test`, `make lint` |

A local MD5 records file identity; it is not automatically publisher verification. Genome gates compare with UCSC manifests. Chain and RepeatMasker hashes are local identity records where publisher checksums were not established. A plausible `.2bit` index or base total alone is insufficient, hence the whole-file checksum check.

The standalone genome-download shell script writes its ledger block after its transfer loop. An interrupted invocation may therefore leave downloads without completed ledger entries; its presence on disk does not substitute for the required-reference gate. The Python download helpers append provenance within each successful transfer. [Download script](scripts/fetch_genomes.sh) · [Download functions](src/fetch_lift_inputs.py)

### Recorded repeatability

Separate paired-run checks compared generated artifacts using SHA-256. These manifests record actual comparisons; they are not automatically regenerated by the Make targets.

| Repeated command | Artifacts compared | Recorded result | Manifest |
| --- | ---: | --- | --- |
| `make lift` | 10 | Byte-identical | [Human conversion](results/tables/liftover_hg19_to_hg38_reproducibility.json) |
| `make pilot-mouse` | 33 | Byte-identical; printed output identical | [Mouse lift](results/tables/pilot_mouse_reproducibility.json) |
| `make pilot-orthologs` | 47 | Byte-identical; printed output identical | [Mouse extraction](results/tables/pilot_mouse_orthologs_reproducibility.json) |

Stable ordering, a fixed pilot sample, and timestamp-free conversion/pilot reports support replay with identical inputs and tools. This does not claim identical bytes across software versions or filesystem layouts. The standalone assembly-evidence report includes its generation timestamp; the byte-identical claims above concern the listed workflows. Input hashes and executed commands are retained in the reports and logs.

## Reproducing

### Environment

Run commands from the repository root:

```bash
conda env create -f environment.yml
conda activate re1-gates
make test
make lint
```

`make env` is equivalent to the environment-creation command. [environment.yml](environment.yml) declares Python dependencies and genomics packages. The completed workflow needs `liftOver`, `twoBitInfo`, and `twoBitToFa`; explicit acquisition uses system `curl`. FIMO, bedtools, and CD-HIT support other helpers and integration tests. Missing executables can cause integration tests to skip, so inspect the test summary rather than interpreting a skip as successful tool verification.

Use `make TARGET PYTHON=/path/to/python` when selecting an interpreter explicitly. Do not run Python with assertions disabled: accounting assertions are part of the validation contract.

### Required local inputs

| Input | Expected location | Purpose |
| --- | --- | --- |
| GEO archive | `data/raw/geo/GSE196171/GSE196171_RAW.tar` | Native intervals and construct measurements |
| Publisher Table S2 | `data/raw/assembly_check/Table_S2_publisher.zip` | Published oligo sequences and explicit aliases |
| Cached reference queries | `data/raw/assembly_check/ucsc/` | Offline source-assembly evidence |
| Human assembly chain | `data/raw/hg19ToHg38.over.chain.gz` | Original-to-working-build conversion |
| RepeatMasker annotation | `data/raw/repeats/hg19.rmsk.txt.gz` | Conversion attrition audit |
| Required genomes | `data/raw/genomes/hg38.2bit`, `data/raw/genomes/mm39.2bit` | Human verification and mouse extraction |
| UCSC reference metadata | `data/raw/genomes/verification/` | Per-assembly `.chrom.sizes` and `md5sum.txt` caches |
| Mouse chain | `data/raw/chains/hg38ToMm39.over.chain.gz` | Element and motif-window projections |
| Provenance ledger | `data/README.md` | Source URLs and cache identity |

A fresh checkout contains code and documentation, not the downloaded references or generated reports. Use the [data catalog](data/README.md) for exact source URLs and recorded identities. The failed `mmc3.zip` response retained in the assembly-check directory is an HTML challenge, not the usable supplement; use `Table_S2_publisher.zip`. [Fetch failure and recovery evidence](results/tables/coordinate_assembly.md)

Acquisition commands are **explicit network operations**, separate from offline reproduction:

```bash
# Requires the GEO archive already at its documented path.
# Fetches the publisher oligo table and small reference-query evidence.
python -m src.assembly_check --fetch

# Fetches the human chain and repeat annotation, without a whole genome.
python -m src.fetch_lift_inputs

# Optional explicit acquisition of hg38 as well.
python -m src.fetch_lift_inputs --genome
```

To acquire the GEO archive or mouse chain through the existing resumable helper, run the following only when those inputs are needed. It checks an existing final file against provenance rather than silently replacing it:

```python
from src import config
from src.assembly_check import GEO_ARCHIVE
from src.fetch_lift_inputs import fetch, remote_size

sources = [
    (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/GSE196171_RAW.tar",
        GEO_ARCHIVE,
    ),
    (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm39.over.chain.gz",
        config.CHAIN_PATHS["mm39"],
    ),
]
for url, path in sources:
    size = remote_size(url)
    print(path.name, size, "bytes")
    fetch(url, path, size)
```

Acquire the verification metadata explicitly for the required pair:

```python
from src.pilot import fetch_metadata

errors = fetch_metadata(("hg38", "mm39"))
if errors:
    raise RuntimeError("; ".join(errors))
```

This metadata helper can also acquire the platform-appropriate UCSC `twoBitInfo` utility if absent. It does not download genomes. The narrower `python -m src.pilot --fetch-metadata` fetches hg38 metadata and then runs the pilot.

Target-genome acquisition is handled manually by [`scripts/fetch_genomes.sh`](scripts/fetch_genomes.sh), which supports `--check-only`, prompts before downloads, and resumes transfers. It covers the target species rather than hg38. Do not launch a competing transfer while that script owns a download. There is no `make data-genomes` target; no analysis target repairs or downloads references.

### Offline replay

With inputs and executables available:

```bash
make assembly-check
make lift
make pilot-orthologs
```

For intermediate inspection, `make pilot` and `make pilot-mouse` are independently available. The last command above rebuilds both before motif-window mapping and extraction. Reports print expected/observed values or explicit errors. Do not bypass a failed gate or substitute an unverified file to obtain an output.

Re-running with identical inputs and tools should reproduce the conversion and pilot artifacts. Compare their hashes and reports as in the [recorded repeatability checks](#recorded-repeatability); the standalone assembly report's timestamp is expected to change. `make test` checks software behavior; real-data assertions in the data targets are separate evidence.

### Data and output handling

Downloaded and generated data remain local and must not be committed. The active [`.gitignore`](.gitignore) excludes data artifacts and most generated results, with documentation/placeholders and fold-output exceptions. Local evidence links can consequently be absent in a fresh checkout. Existing reports document observed runs; regenerating them requires the corresponding inputs.

Before committing source or documentation changes, inspect:

```bash
git status --short
git diff --cached --name-only -- data/
```

Do not stage downloaded/generated data or change historical ledger rows. Download provenance and generated-output identity are separate records: the ledger records acquisition; reports, audits, and repeatability manifests identify derived artifacts.

## Repository layout

| Location | Contents and role |
| --- | --- |
| [`README.md`](README.md) | Technical record, measured results, reproduction instructions, and interpretation boundaries |
| [`README (2).md`](README%20%282%29.md) | Research rationale and proposed study design; use this README for implementation status |
| [`data-README.md`](data-README.md) | Acquisition and data-preparation companion |
| [`data/README.md`](data/README.md) | Local data catalog and append-only download provenance |
| [`src/config.py`](src/config.py) | Assembly identities, paths, pilot chromosome, motif convention, and chosen mouse minMatch; some future-analysis settings remain unset |
| [`src/assembly_check.py`](src/assembly_check.py), [`src/lift.py`](src/lift.py) | Source-assembly evidence, generic liftOver interface, verified human conversion, and repeat-bias audit |
| [`src/pilot.py`](src/pilot.py), [`src/pilot_mouse.py`](src/pilot_mouse.py), [`src/pilot_orthologs.py`](src/pilot_orthologs.py) | Dependency-scoped genome gates, canonical pilot, mouse sweep, motif-window calls, and sequence extraction |
| [`src/downloads.py`](src/downloads.py), [`src/fetch_lift_inputs.py`](src/fetch_lift_inputs.py), [`scripts/`](scripts/) | Explicit acquisition and provenance; manual target-genome download script |
| [`src/_utils.py`](src/_utils.py), [`src/io_mpra.py`](src/io_mpra.py) | Table, coordinate, FASTA, and input-validation helpers; explicit activity-column mappings remain unset |
| [`src/motifs.py`](src/motifs.py), [`src/conservation.py`](src/conservation.py), [`src/orthologs.py`](src/orthologs.py), [`src/classify.py`](src/classify.py) | Scoring, remote-track, orthology, and classification helpers; production comparative features are unfinished |
| [`src/folds.py`](src/folds.py), [`src/models.py`](src/models.py), [`src/metrics.py`](src/metrics.py) | Fold grouping and evaluator; explicit training and gate-metric stubs |
| [`tests/`](tests/) | Unit tests and external-tool integration fixtures |
| [`notebooks/`](notebooks/) | Unexecuted exploration scaffold; not a source of completed scientific results |
| `data/raw/` | Source archives, sequence-query caches, reference genomes, chains, annotations, and verification metadata |
| `data/processed/` | Verified human CSV, pilot BED, lift outputs, per-element audits, and extracted mouse sequences |
| `results/tables/` | Measured reports, assembly comparison rows, and repeatability manifests |
| `results/figures/`, `results/folds/` | Destinations for future figures and real-data fold assignments |
| [`Makefile`](Makefile), [`environment.yml`](environment.yml) | Supported commands and dependency specification |

## Data sources

| Source | Role | Status and local evidence |
| --- | --- | --- |
| [GEO GSE196171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171) and [supplementary directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/) | MPRAduo activity tables and native interval metadata | Inspected; assembly-resolved coordinate inputs used in human conversion; activity preparation unfinished |
| [Publisher Table S2 archive](https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22001926-mmc3.zip) | Library oligo sequences | Used in [assembly discrimination](results/tables/coordinate_assembly.md) and [full human verification](results/tables/liftover_hg19_to_hg38.md) |
| [UCSC sequence API](https://genome.ucsc.edu/goldenPath/help/api.html) | Small reference queries at published coordinates | Cached observations and query hashes in [comparison CSV](results/tables/coordinate_assembly_comparisons.csv) |
| [UCSC hg38 reference](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) and [mm39 reference](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/) | `.2bit` DNA, chromosome sizes, publisher checksums | Required references pass the [mouse gate](results/tables/pilot_mouse_lift.md) |
| [UCSC hg19 chains](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) and [hg38 chains](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) | Human assembly conversion and comparative interval projection | hg19→hg38 and hg38→mm39 used; exact input digests in their reports |
| [UCSC hg19 RepeatMasker](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz) | Repeat-overlap diagnostic for human conversion | Used in [attrition comparison](results/tables/liftover_hg19_to_hg38.md) |
| [JASPAR](https://jaspar.elixir.no/) | Intended REST motif input | Production matrix/scoring validation remains outstanding |
| [ENCODE](https://www.encodeproject.org/) and publication supplements | Intended cofactor/accessibility evidence and accession reconciliation | No validated production feature set reported |
| UCSC phyloP/phastCons tracks | Intended conservation features | Remote-reading helper exists; specific track URLs remain unset |
| [Author duoREST repository](https://github.com/tewhey-lab/duoREST) | Reference implementation for understanding the published analysis | External reference, not an executed dependency of this workflow |

## References

- **Mouri et al.** *Whole-genome functional characterization of RE1 silencers using a modified massively parallel reporter assay.* Primary assay and sequence source. [Article and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903721/) · [DOI](https://doi.org/10.1016/j.xgen.2022.100234)
- **UCSC Genome Browser.** Reference assemblies, chains, sequence API, and liftOver parameter semantics. [Downloads](https://hgdownload.soe.ucsc.edu/downloads.html) · [API](https://genome.ucsc.edu/goldenPath/help/api.html) · [minMatch explanation](https://groups.google.com/a/soe.ucsc.edu/g/genome/c/YbyOe7AEXj4)
- **MEME Suite FIMO.** Motif-scanning software used by the implemented wrapper; production scoring is pending. [Documentation](https://meme-suite.org/meme/doc/fimo.html)
- **pyBigWig.** Indexed remote conservation-track access. [Project documentation](https://github.com/deeptools/pyBigWig)
- **CD-HIT.** Sequence-clustering software used by the fold-construction helper. [Project documentation](https://github.com/weizhongli/cdhit)

Project measurements in this README are linked to files in `results/tables/`. The conversion-attrition percentage is explicitly calculated from those reported counts; it is not an additional measurement. Identifiers such as assembly names and accessions identify inputs, not experimental results.

## Scope and known limitations

| Limitation | Consequence for interpretation |
| --- | --- |
| Episomal reporter versus integrated circuit | Reporter repression does not establish behavior after genomic integration, chromatin remodeling, or coupling to other circuit components |
| Cell lines versus neurons | The assayed cell lines do not demonstrate performance in mature neurons, organoids, or neural computing systems |
| Candidate list versus demonstrated gate | A computational ranking would nominate sequences for testing; switching behavior, dynamic range, leakiness, and circuit robustness remain unmeasured here |
| ABSENT is alignment loss | Failure of both lift tracks can reflect chain coverage, assembly differences, or evolutionary change; it does not prove a biological deletion or human-specific origin |
| MAPPED is not motif retention | Partial motif-window alignment and sequence changes remain possible; a successful lift does not establish REST binding |
| Target-span distortion | Long or shortened spans can complicate rescoring; width flags remain available for sensitivity analyses |
| Pilot selection | A canonical-element chromosome pilot establishes an executable path, not generalization to other chromosomes, contrast elements, or species |
| Assay outcome not finalized | Enhancer-context handling, missingness, repression sign, and positive-label definitions must be settled before fitting models; left-tail repression has not yet been validated |
| Repeat overlap is a proxy | The conversion diagnostic does not directly label primate-specific elements or rule out differential evolutionary-class attrition |
| Conservation is observational | Predictive improvement, if found, would not establish that evolutionary retention causes repression |
| Evaluation remains unrun | Grouped folds and metric functions cannot establish out-of-sample performance until real inputs, models, and held-out predictions exist |

These boundaries separate verified inputs and mapping evidence from proposed comparative interpretation and eventual experimental testing. The [research-design companion](README%20%282%29.md) develops the hypothesis and proposed analyses.

## Team

**Gabrielle DeCuir and Vrinda Sehgal** · Johns Hopkins Applied Comparative Genomics.

The project builds on the published MPRAduo study and public reference resources. Repository validation results are this project's work; source-study measurements remain attributed to their authors.
