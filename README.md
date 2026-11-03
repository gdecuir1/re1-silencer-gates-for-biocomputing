# re1-silencer-gates-for-biocomputing
Mining human RE1 silencer elements as candidate genetic NOT gates for neural biocomputing. Comparative genomics, MPRA data, and ML to rank REST-binding sites by predicted gate quality. JHU Applied Comparative Genomics, Fall 2026.

## Scaffold status

Data loading/validation, FIMO, remote conservation queries, liftOver, sequence
extraction, state classification, fold generation, and held-out metric evaluation
are implemented. Model training and gate-quality metrics deliberately raise
`NotImplementedError`. No research data, computed folds, or scientific results are
provided. Tests use explicitly synthetic fixtures in memory or temporary folders.

## Environment and checks

Use Python 3.11 and the dependencies in `environment.yml`:

```bash
make env
conda activate re1-gates
make test
make lint
python -c "import src.config"
```

`make env` explicitly installs software, not datasets. If the environment already
exists, use `conda env update -f environment.yml` instead. `pytest` and `ruff` are
the added development dependencies. For an alternate interpreter, use
`make test lint PYTHON=python`. `make clean` removes tool/Python caches only.
External-tool tests run when FIMO, CD-HIT-EST, liftOver, and bedtools are on PATH;
otherwise those checks skip. Other tests require no executables or network access.

## Configure real inputs

All default paths derive from the location of `src/config.py`, independently of
the working directory. Nothing downloads data automatically.

1. Put the verified MPRA file at `MPRA_PATH` (default `data/raw/mpra.tsv`), or pass
   its path to `load_mpra`. Record provenance in `data/README.md`.
2. Fill `MPRA_COLUMNS` with canonical-to-source headers for `element_id` and
   `sequence`. For wide tables, fill `MPRA_REPRESSION_COLUMNS` with one source
   repression header for each of the four configured cell types. For long tables,
   set `MPRA_LAYOUT = "long"` and also map `cell_type` and `repression` in
   `MPRA_COLUMNS`. File delimiter is `MPRA_SEPARATOR`.
3. Map optional `chrom`, `start`, `end`, and `strand` together if available, or
   join independently verified coordinates by element ID before genomic queries.
   Coordinates must already be zero-based, half-open hg38 intervals. Sequence
   orientation must be explicitly `+` or `-`; the loader never guesses it.
4. Supply the existing JASPAR MA0138.2 **MEME-format** file at `MOTIF_PATH`, and
   existing chain/genome files at the configured paths or explicit function paths.
5. Fill `BIGWIG_URLS` with verified hg38 phyloP and phastCons HTTP(S) URLs and
   record the chosen track versions in `data/README.md`. No track is assumed.

The loader requires exactly 8,436 unique 200 bp elements and one finite repression
measurement per element/cell pair. It rejects missing cells, duplicates, inconsistent
metadata, invalid DNA, and malformed coordinates. It does not select a source-table
subset, average replicates, change signs, or derive repression from raw activity.
Those scientific decisions require verification of the source schema and assay scale.

## Coordinate and state conventions

`MOTIF_SPAN = (91, 111)` is **zero-based inclusive**, so the motif spans 21 bases.
For forward-oriented elements the bigWig interval is `[start + 91, start + 112)`;
reverse-oriented elements use `[end - 112, end - 91)`. The orientation is that of
the supplied element sequence, not the strand of its best FIMO hit.

FIMO runs on both strands as specified by the supplied matrix; provide a DNA matrix
with `strands: + -`. Its text output is streamed to a temporary file and reduced in
chunks with p-value threshold 1, avoiding default p-value filtering and stored-hit
limits. Best hits maximize score, with deterministic ties; FIMO coordinates remain
one-based inclusive. `retains_motif` uses `score >= 20.86`. Unscorable sequences retain
missing scores and missing motif-retention calls. A matrix's background frequencies
and scoring settings must be checked against the published cutoff before scientific
interpretation. See the [FIMO documentation](https://meme-suite.org/meme/doc/fimo.html).

Conservation uses pyBigWig/libBigWig remote reads only after a one-byte HTTP range
probe succeeds. Means are calculated from the actual 21 per-base values, with covered
base counts alongside them. Missing values remain missing; chromosome mismatches and
out-of-bounds windows raise. See [pyBigWig documentation](https://github.com/deeptools/pyBigWig).

`lift_elements(elements, chain, min_match=...)` requires an explicit minimum mapped
fraction and sends the **full 200 bp** source interval to liftOver. It returns mapped
and unmapped BED6 tables; failed records include original coordinates, raw BED lines,
and reasons. Target intervals may differ in length. `extract_sequences` uses stranded
bedtools extraction and checks for skipped/truncated records. Chain provenance,
cross-species mapping quality, and target motif correspondence still require analysis;
a best hit anywhere in the lifted interval is not itself proof of motif orthology.

For classification, supply one explicit state per target species:

| State | Meaning |
| --- | --- |
| `retained` | Lifted and motif retention established |
| `lifted_no_motif` | Lifted and scored, with no retained motif |
| `failed_lift` | No mapped interval; not biological absence |
| `lifted_unscorable` | Mapped, but motif presence cannot be evaluated |

The deepest observed retention wins: opossum supports `therian_retention`;
mouse/dog/cow support `eutherian_retention`; macaque supports `primate_retention`.
The explicit rule table is in config. These are descriptive evidence classes,
not inferred divergence dates. With no retained motif, incomplete evidence is
`unresolved`; otherwise the label is `no_observed_outgroup_motif`. No call asserts
human specificity from a failed lift. The matrix output preserves all original states
and adds an evidence-completeness flag.

## Leakage guard and evaluation

Pass one row per element (`element_id`, `chrom`, `sequence`) to `build_folds`.
It runs nucleotide CD-HIT (`cd-hit-est`) on both strands, then joins chromosomes
connected by any shared cluster, including transitive connections. Entire connected
components go to the same fold. Both chromosomes and CD-HIT clusters are disjoint
between train and test. Too few components raises instead of weakening the grouping.

The configurable defaults (95% global identity, 95% coverage of each sequence,
five folds) are operational choices, not verified biological cutoffs. Assignments
are deterministic and save cluster IDs, component IDs, and thresholds to
`results/folds/assignments.csv`. Retain tool versions and input provenance with any
real run. [CD-HIT uses greedy representative clustering](https://github.com/weizhongli/cdhit);
the guarantee is that its clusters never split across folds, not that every possible
pair above a similarity threshold has been exhaustively identified. An all-pairs
similarity audit would be a separate validation step.

`evaluate_folds` accepts held-out predictions and explicit binary positive labels;
it returns Spearman and AUROC per fold and cell type. Higher predictions must mean
stronger repression/more likely positive. Undefined statistics remain NaN with a
reason. The evaluator cannot establish that predictions were actually held out;
future model implementations must use the saved folds and training-only preprocessing.

## Repository contents

```text
src/                  Configuration, mechanical processing, and explicit model stubs
tests/                Data-free tests and optional installed-tool integration checks
notebooks/01_explore.ipynb  Load verified MPRA data and plot repression by cell type
data/README.md        Source provenance inventory
data/raw/             Local source files (ignored)
data/processed/       Derived inputs (ignored)
results/figures/       Generated figures (ignored)
results/tables/        Generated tables (ignored)
results/folds/         Reproducible fold definitions (allowed in Git)
```

The original repository contained `gitignore`, which Git does not apply. The active
`.gitignore` now excludes data and generated results, allowing only README files,
`.gitkeep` placeholders, and fold definitions. Large genomic formats remain ignored.
The original notes (`README (2).md`, `data-README.md`, and `gitignore`) are retained.
