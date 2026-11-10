# Data Catalog and Provenance Ledger

## Week 2 motif-window coverage and mouse sequences — 2026-11-08

The independent 21 bp motif-window track uses the same minMatch 0.5 as the element track. The sidecar accounts for **217 MAPPED + 15 UNDETERMINED + 158 ABSENT = 390**. The 15 UNDETERMINED cases comprise ten element-only and five motif-only lifts; the latter policy was explicitly approved. All 15 are excluded for mouse downstream and absent from the FASTA. ABSENT denotes an alignment-based candidate loss, not a confirmed biological absence.

The FASTA contains **217** oriented mm39 element sequences with 15 bp padding on each side. **60** MAPPED elements have original widths outside 160–240 bp. Sequence lengths are **136–3,285 bp**, with **zero** under 50 bp, **zero** containing N, and **zero** requiring chromosome-boundary padding clipping. No additional downloads were needed.

See the [full report](../results/tables/pilot_mouse_orthologs.md), [data contract](../data-README.md#motif-window-coverage-and-ortholog-fasta), and [reproduction evidence](../results/tables/pilot_mouse_orthologs_reproducibility.json). Two real `make pilot-orthologs` runs produced 47 byte-identical artifacts and identical printed output.

New derived files (SHA-256; source identities are recorded in the report and inputs.json):

| Derived file | SHA-256 |
| --- | --- |
| `data/processed/pilot_mouse/orthologs.mm39.fa` | `5f8092ad6a5a38a91bf918a61ee7149682530b7ffb34683eb6ba368b01b6d892` |
| `data/processed/pilot_mouse/orthologs.mm39.csv` | `85117ff6a52bae9c6e4e40b9a3b176034b2e83d035a6bc02551650d953a10180` |
| `data/processed/pilot_mouse/motif_window/input.hg38.bed` | `74639edaa17f9ee18b23b92f0fc76360a1a80424eb689d68321a4d5dbb684c42` |
| `data/processed/pilot_mouse/motif_window/liftOver.log` | `3ba5adda83a5eb97ffdb4a9cdd71188c9cb322a664889816946c335adf409c0b` |
| `data/processed/pilot_mouse/motif_window/lifted.mm39.bed` | `b4369aac59b834755553d13a5431341dbc0e5637abbdcd81edfa991908023efa` |
| `data/processed/pilot_mouse/motif_window/unmapped.hg38.bed` | `d9ee73956a21773767969436850c54ad12c68bc89aff3d2855c8cb50c3930781` |
| `data/processed/pilot_mouse/ortholog_extraction/element_track/liftOver.log` | `079308041fdb4dee3a0fffeac6eadbdfb159cc4893665f04b056c8dc7f1b1689` |
| `data/processed/pilot_mouse/ortholog_extraction/element_track/lifted.mm39.bed` | `bf898c197ead2ad7c66c6b0a4bb939c679b27cc5302f1449cc324230894d891d` |
| `data/processed/pilot_mouse/ortholog_extraction/element_track/unmapped.hg38.bed` | `a7f7097fb1e8a227269cd267665e050e54ad4297edf4b2f337ea0f31ed7bf612` |
| `data/processed/pilot_mouse/ortholog_extraction/expected_pilot.hg38.bed` | `0b52b57f2676d476a8a4a131a66ecc2cff187dbb2f27a8773b37b474b16e5930` |
| `data/processed/pilot_mouse/ortholog_extraction/inputs.json` | `3d79eb85f22b2f8dd74d7e8f0caa65e49fe04adafae912d00e3683e566d7573f` |
| `data/processed/pilot_mouse/ortholog_extraction/padded_elements.mm39.bed` | `8e53b1e269f2db6b8222251361bd5252abac57fe13cb8540fdc952c6b0330bcb` |
| `data/processed/pilot_mouse/ortholog_extraction/twoBitToFa.log` | `d3c8a825e386e15b5af05748d96690cf8a0356eea4a73057d0b1882e9172e111` |

Existing provenance and derivation snapshots below are preserved unchanged.


## Week 2 chr8 pilot and mouse lift — 2026-11-08

The pilot is now published: **390 canonical chr8 elements**, all 200 bp, with **10/10 exact hg38 regressions**. Genome dependencies are scoped by consumer: hg38 for the pilot; hg38/mm39 for the mouse lift; hg38 plus all five targets for the future multi-species step. The integrity verifier itself is unchanged. This work did not inspect or change bosTau9 or monDom5 downloads.

The measured mouse sweep selected **minMatch 0.5**: **227 mapped and 163 unmapped**, with 64 mapped spans flagged for >20% length deviation. All outcomes remain in the audit. The new hg38-to-mm39 chain's source URL, UTC completion timestamp, byte size, and local MD5 were appended automatically to the download ledger below. No old provenance row was edited.

See the [pilot report](../results/tables/pilot_setup.md), [mouse sweep and failure reasons](../results/tables/pilot_mouse_lift.md), and [two-run reproduction evidence](../results/tables/pilot_mouse_reproducibility.json). All 33 artifacts and printed output were identical across two real-data runs.

Derived-file inventory for this run (SHA-256, not download MD5):

| Derived file | SHA-256 |
| --- | --- |
| `data/processed/pilot_chr8.bed` | `0b52b57f2676d476a8a4a131a66ecc2cff187dbb2f27a8773b37b474b16e5930` |
| `data/processed/pilot_mouse/expected_pilot.hg38.bed` | `0b52b57f2676d476a8a4a131a66ecc2cff187dbb2f27a8773b37b474b16e5930` |
| `data/processed/pilot_mouse/final/element_audit.csv` | `fa26d21cda12f1ece4f2a93bc1d56bc80e1d5818e04b28552bf24e133b11b826` |
| `data/processed/pilot_mouse/final/liftOver.log` | `dc92bd0e2672c4302ab84d71ff3ea3cd0f44fdc6758069e88d0534cf5b037a58` |
| `data/processed/pilot_mouse/final/lifted.mm39.bed` | `bf898c197ead2ad7c66c6b0a4bb939c679b27cc5302f1449cc324230894d891d` |
| `data/processed/pilot_mouse/final/unmapped.hg38.bed` | `a7f7097fb1e8a227269cd267665e050e54ad4297edf4b2f337ea0f31ed7bf612` |
| `data/processed/pilot_mouse/minmatch_0.1/element_audit.csv` | `bebeea11f033f33143a22280dbc9c2bfe7fb6b5b68d5c578fc89a9a86afec3fa` |
| `data/processed/pilot_mouse/minmatch_0.1/liftOver.log` | `570e970856baaca3f4b875f33c6f98fd604dfba3a71f2ae3e72c78c1b7035ef5` |
| `data/processed/pilot_mouse/minmatch_0.1/lifted.mm39.bed` | `f4e25d74c471ae7e4facc83587d97cd051d25ba7b3289f1bc5e209593b14517f` |
| `data/processed/pilot_mouse/minmatch_0.1/unmapped.hg38.bed` | `8b1a151284ee67bc730f030c7c0eda537d47d9789993a70d64fdb92d2adca0d5` |
| `data/processed/pilot_mouse/minmatch_0.2/element_audit.csv` | `28be399d6438b49c83ed3124e8dc31f85a0d18f771f2ce5b62a5605340340e96` |
| `data/processed/pilot_mouse/minmatch_0.2/liftOver.log` | `2497fc11ccafde916113b0c53425ae3a5154be20dc9efe88f83f91739b7c2e12` |
| `data/processed/pilot_mouse/minmatch_0.2/lifted.mm39.bed` | `da8a34e4e2a379a48fc787566dda61be363ae87714b955874d23c2aebc73fbc3` |
| `data/processed/pilot_mouse/minmatch_0.2/unmapped.hg38.bed` | `648f6ce258523b007863b51949d36235829629fab0810b79449c56962dccf3aa` |
| `data/processed/pilot_mouse/minmatch_0.5/element_audit.csv` | `fa26d21cda12f1ece4f2a93bc1d56bc80e1d5818e04b28552bf24e133b11b826` |
| `data/processed/pilot_mouse/minmatch_0.5/liftOver.log` | `14c9c48efa8e337da68d059629e28a45ce04b4024ece8b275f27129cfd622139` |
| `data/processed/pilot_mouse/minmatch_0.5/lifted.mm39.bed` | `bf898c197ead2ad7c66c6b0a4bb939c679b27cc5302f1449cc324230894d891d` |
| `data/processed/pilot_mouse/minmatch_0.5/unmapped.hg38.bed` | `a7f7097fb1e8a227269cd267665e050e54ad4297edf4b2f337ea0f31ed7bf612` |
| `data/processed/pilot_mouse/minmatch_0.95/element_audit.csv` | `112125dfabc1b5376294f0a23adfa6fe6899203aee127af89507872c651ab385` |
| `data/processed/pilot_mouse/minmatch_0.95/liftOver.log` | `246bc3c1a8ab4ef1546f20014e2d687ac2e19d441e7bce44deb866be6ce545e2` |
| `data/processed/pilot_mouse/minmatch_0.95/lifted.mm39.bed` | `bad9a9ff32e146c3569a2861952ff96fb85ab2a9e83e2f391ed45f483bb450e1` |
| `data/processed/pilot_mouse/minmatch_0.95/unmapped.hg38.bed` | `77530c8d30cc55b0d5c5dc320f1c880cbaac0b21abdf7719a72713a9627b6f14` |
| `data/processed/pilot_mouse/sweep.csv` | `890fd944414ad9cbae94d177bebead874e2cadba3ffe3c7c930be24b720b017e` |
| `data/processed/pilot_setup/hg38.twoBitInfo.log` | `912f6471944218f892fbf282321c1f77fd043818aa87ddb74ad1a48827691396` |
| `data/processed/pilot_setup/hg38.twoBitInfo.txt` | `d06fcd768c9451398eb8e37d486c94e464efc5aca4991af4b05a423dcce3cd6f` |
| `data/processed/pilot_setup/mm39.twoBitInfo.log` | `91db5b9198a6e978b6aad86838ad7707a6337569627db90e6cf16aa8bd065493` |
| `data/processed/pilot_setup/mm39.twoBitInfo.txt` | `90be49c9da524f66246749e3d8fe41f9a124d6cc5989bdd6640babcab76de894` |
| `data/processed/pilot_setup/regression.csv` | `3d48d39e18d7ff61c44d4218e9f4db340b6843556742b61979234db8355ffc13` |
| `data/processed/pilot_setup/sample.bed` | `6afefca381cc62b4149c8aff59a9a08da3bd5a01677af4c11da50e1e307d51b2` |
| `data/processed/pilot_setup/sample.fa` | `a96eb8b86f2b6403a340bb5485777fa4e95912133a3949623cfb90a6b4e3f45c` |
| `data/processed/pilot_setup/twoBitToFa.log` | `a5ebb55ca8bc95393f14b5fc9af9ff46fa980bc7de31579a5dbb47a54475c081` |

These are derivations, not downloads. Existing rows below remain dated historical snapshots.


## Coordinate conversion update — 2026-11-08

This update supersedes older conversion/name-reconciliation status statements in the **2026-11-03 catalog snapshot below**. All existing catalog and provenance rows are retained verbatim. New downloads append at the end of the existing ledger.

**12,815 verified, 47 lifted with sequence mismatches, and 4 failed lifts out of 12,866 inputs.** All 12,866 library sequences were resolved by exact names or explicit publisher aliases. No adapter trimming was needed. Downstream coordinates come from `data/processed/elements_hg38.csv`, with source hg19 coordinates retained only for traceability. See the [operating guide](../data-README.md#verified-hg38-handoff) and [complete report](../results/tables/liftover_hg19_to_hg38.md).

New reference acquisitions: `data/raw/hg19ToHg38.over.chain.gz`, `data/raw/repeats/hg19.rmsk.sql`, `data/raw/repeats/hg19.rmsk.txt.gz`, `data/raw/genomes/hg38.md5sum.txt`, and `data/raw/genomes/hg38.2bit`. Their exact source URLs, timestamps, sizes, and local MD5s are appended below. hg38's local MD5 also matches UCSC's published digest. The hg19 database `md5sum.txt` request returned HTTP 404; no upstream checksum verification is claimed for RepeatMasker. Initial sandboxed requests failed DNS resolution and were retried with approved network access.

Generated data inventory (not downloads; derived by `make lift` from the recorded source inputs):

| Generated file | MD5 after verified replay |
| --- | --- |
| `data/processed/elements_hg38.csv` | `ea9e93a9c9a831cac9c134ee9cfb3474` |
| `data/processed/liftover_hg19_to_hg38/element_audit.csv` | `637fb18a5f003e3f01d8f70e59e74122` |
| `data/processed/liftover_hg19_to_hg38/input.hg19.bed` | `e0406ca2118970454dc8fdebf73dd0a1` |
| `data/processed/liftover_hg19_to_hg38/input_elements.csv` | `52a8e5087eef6c0cb304e43c12ded227` |
| `data/processed/liftover_hg19_to_hg38/liftOver.log` | `c3807eb996f640a33f5568e4b0689bb2` |
| `data/processed/liftover_hg19_to_hg38/lifted.hg38.bed` | `83e5f99992e94aed4737da435933f70e` |
| `data/processed/liftover_hg19_to_hg38/lifted.hg38.fa` | `e11dc983f0021c15183a56d2a9152a4a` |
| `data/processed/liftover_hg19_to_hg38/twoBitToFa.log` | `443ef84e1ec966fcf1dafe64023c0463` |
| `data/processed/liftover_hg19_to_hg38/unmapped.hg19.bed` | `41a56f03508f55db8a842e10aa168f0e` |

The full per-element audit preserves failures and mismatches. Download ledger records are immutable; derived-file hashes describe this conversion snapshot. A future change to input data or code requires a new run and a new dated inventory, not rewriting these rows.

---


**Authoritative record of local source files, acquisition evidence, and download provenance.**

[Main operating guide](../README.md) · [Research rationale](../README%20%282%29.md) · [Week 1 workflow](../data-README.md)

Catalog status: **2026-11-03**. Data payloads are local and ignored by Git; this README is versioned documentation. Existing download records are historical evidence and must not be rewritten to match a new file or interpretation.

> **Assembly:** native GEO `motif` and `noMotif` intervals are supported as **hg19**. The project's hg38 target is a separate coordinate system. Source intervals have not been converted. The resolved assembly does not imply that every element, variant, or library has been sequence-validated.

## Contents

1. [Catalog scope and status](#catalog-scope-and-status)
2. [Storage map](#storage-map)
3. [GEO supplementary inventory](#geo-supplementary-inventory)
4. [Table schemas and missing values](#table-schemas-and-missing-values)
5. [Native population and identifier semantics](#native-population-and-identifier-semantics)
6. [Published sequence supplement](#published-sequence-supplement)
7. [Assembly evidence cache](#assembly-evidence-cache)
8. [Acquisition failures and content validation](#acquisition-failures-and-content-validation)
9. [Pending inputs and configured locations](#pending-inputs-and-configured-locations)
10. [Provenance contract](#provenance-contract)
11. [Historical planning inventory](#historical-planning-inventory)
12. [Phase A download ledger](#phase-a-discovery-downloads)
13. [Assembly-investigation download ledger](#coordinate-assembly-investigation-downloads)

## Catalog scope and status

The original discovery and assembly follow-up acquired **101 files/responses totaling 33,071,763 bytes**. This total excludes README files and `.gitkeep` placeholders. The acquisition record includes one failed HTML response saved under a `.zip` filename; file presence and a checksum are not proof that its content is valid scientific input.

| Evidence category | State | What it establishes |
| --- | --- | --- |
| GEO directory inventory | Inspected and saved | Three supplementary-directory entries and their observed sizes |
| Six tables inside the GEO archive | Read directly from gzip members | Row counts, field names, inferred types, native populations, and interval metadata |
| Noncanonical K562 results table | Inspected separately | A distinct processed table with entirely missing coordinate columns |
| Table S2 publisher ZIP | Downloaded and opened | Oligo names, sequence strings, and library labels |
| UCSC sequence responses | Downloaded as bounded windows | Reference DNA for the sampled assembly comparisons |
| Prepared MPRA modeling input | Not produced | Population/context/sign decisions remain unresolved |
| Whole reference genomes | Not downloaded | Full reference-based downstream processing remains pending |

A fresh checkout contains the documentation, not these ignored payloads. Acquisition instructions are in the [Week 1 guide](../data-README.md#acquisition-and-replay-procedures).

## Storage map

Paths below are relative to the repository root. Links from this README are relative to `data/`.

| Repository path | Contents | Treatment |
| --- | --- | --- |
| `data/raw/geo/GSE196171/` | Archive, processed table, inventory, metadata, directory snapshot | Retain original downloaded bytes |
| `data/raw/assembly_check/Table_S2_publisher.zip` | Valid published sequence supplement | Retain archive; inspect its member directly |
| `data/raw/assembly_check/mmc3.zip` | Failed PMC HTML response | Preserve failure evidence; never parse as the sequence supplement |
| `data/raw/assembly_check/ucsc/` | 94 assembly-specific JSON responses | Verify URL/assembly/coordinates/hash before reuse |
| `data/processed/` | Future prepared inputs | Currently contains only `.gitkeep` |
| `results/tables/` | Generated discovery and assembly reports | Derived artifacts, not original downloads |
| `results/figures/` | Future validation/scientific figures | No figures generated yet |
| `results/folds/` | Future reproducible fold definitions | No real-data assignments generated yet |

The six gzip tables inside `GSE196171_RAW.tar` were not separately extracted to disk. They have one archive acquisition record, plus the member inventory below. Do not count archive members as six additional downloaded files.

## GEO supplementary inventory

Source: [NCBI supplementary directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/).

| Directory entry | Observed bytes | Local file | Meaning |
| --- | ---: | --- | --- |
| `GSE196171_RAW.tar` | 31,252,480 | [Archive](raw/geo/GSE196171/GSE196171_RAW.tar) | Container for six compressed tables; its RAW filename does not mean all fields are unprocessed counts |
| `GSE196171_nonCan_K562_results.run1.txt.gz` | 308,452 | [Noncanonical results](raw/geo/GSE196171/GSE196171_nonCan_K562_results.run1.txt.gz) | Separate K562 experiment |
| `filelist.txt` | 596 | [NCBI inventory](raw/geo/GSE196171/filelist.txt) | Archive/member metadata |

Two additional discovery snapshots were saved: [series metadata](raw/geo/GSE196171/series_metadata.soft.txt), 3,685 bytes, and [directory HTML](raw/geo/GSE196171/supplementary_directory.html), 832 bytes. These are metadata responses, not additional supplementary-directory datasets.

### Archive members

The byte counts below are sizes of the compressed gzip members inside the TAR, not decompressed text sizes. Table row counts exclude headers.

| Member | Compressed bytes | Rows | Columns | Experiment |
| --- | ---: | ---: | ---: | --- |
| `GSM5861606_MPRAduo_counts_benchmark_run1.txt.gz` | 1,282,800 | 27,169 | 20 | Benchmark run 1 |
| `GSM5861607_MPRAduo_counts_benchmark_run2.txt.gz` | 1,511,379 | 27,169 | 20 | Benchmark run 2 |
| `GSM5861611_MPRAduo_counts_WGREST_GM12878.txt.gz` | 7,144,960 | 139,325 | 20 | WGREST / GM12878 |
| `GSM5861612_MPRAduo_counts_WGREST_HepG2.txt.gz` | 7,125,969 | 139,325 | 20 | WGREST / HepG2 |
| `GSM5861613_MPRAduo_counts_WGREST_K562.txt.gz` | 7,086,099 | 139,325 | 20 | WGREST / K562 |
| `GSM5861614_MPRAduo_counts_WGREST_SKNSH.txt.gz` | 7,094,836 | 139,325 | 20 | WGREST / source label `SK.N.SH` |

The separate noncanonical table has **4,975 rows and 20 columns**. It is not a substitute for the canonical 8,436-element population.

## Table schemas and missing values

The tables below enumerate every original column, grouped by inferred dtype. The source files are read with `pandas.read_csv(..., sep="\t", low_memory=False)` and default NA recognition. Dtypes describe parser output, not declared database types or biological measurements. In particular, `float64` can represent integer-valued counts with missing entries, and an entirely missing column also infers a floating dtype.

The original discovery used pandas 3.1.0; the schema tables below were regenerated from the local files during this documentation update. The detailed [discovery report](../results/tables/week1_discovery.md) also records each column's nonmissing count and original order.

### GSM5861606_MPRAduo_counts_benchmark_run1.txt.gz

| Rows | Observed dtype | Columns with this dtype |
| ---: | --- | --- |
| 27,169 | `str` | `id`, `library`, `enhancer`, `silencer`, `chromosome`, `start`, `end`, `project` |
| 27,169 | `int64` | `cDNA_01`, `cDNA_02`, `cDNA_03`, `cDNA_04`, `plmd_01`, `plmd_02`, `plmd_03`, `plmd_04`, `PlasmidsBCsum` |
| 27,169 | `float64` | `plmean`, `log2FoldChange`, `lfcSE` |

### GSM5861607_MPRAduo_counts_benchmark_run2.txt.gz

| Rows | Observed dtype | Columns with this dtype |
| ---: | --- | --- |
| 27,169 | `str` | `id`, `library`, `enhancer`, `silencer`, `chromosome`, `start`, `end`, `project` |
| 27,169 | `int64` | `cDNA_01`, `cDNA_02`, `cDNA_03`, `cDNA_04`, `plmd_01`, `plmd_02`, `plmd_03`, `plmd_04`, `PlasmidsBCsum` |
| 27,169 | `float64` | `plmean`, `log2FoldChange`, `lfcSE` |

### GSM5861611_MPRAduo_counts_WGREST_GM12878.txt.gz

| Rows | Observed dtype | Columns with this dtype |
| ---: | --- | --- |
| 139,325 | `str` | `id`, `celltype`, `enhancer`, `silencer`, `chromosome`, `project` |
| 139,325 | `float64` | `start`, `end`, `GM12878_1`, `GM12878_2`, `GM12878_3`, `GM12878_4`, `plasmid_1`, `plasmid_2`, `plasmid_3`, `plasmid_4`, `PlasmidsBCsum`, `plmean`, `log2FoldChange`, `lfcSE` |

### GSM5861612_MPRAduo_counts_WGREST_HepG2.txt.gz

| Rows | Observed dtype | Columns with this dtype |
| ---: | --- | --- |
| 139,325 | `str` | `id`, `celltype`, `enhancer`, `silencer`, `chromosome`, `project` |
| 139,325 | `float64` | `start`, `end`, `HepG2_1`, `HepG2_2`, `HepG2_3`, `HepG2_4`, `plasmid_1`, `plasmid_2`, `plasmid_3`, `plasmid_4`, `PlasmidsBCsum`, `plmean`, `log2FoldChange`, `lfcSE` |

### GSM5861613_MPRAduo_counts_WGREST_K562.txt.gz

| Rows | Observed dtype | Columns with this dtype |
| ---: | --- | --- |
| 139,325 | `str` | `id`, `celltype`, `enhancer`, `silencer`, `chromosome`, `project` |
| 139,325 | `float64` | `start`, `end`, `K562_1`, `K562_2`, `K562_3`, `K562_4`, `plasmid_1`, `plasmid_2`, `plasmid_3`, `plasmid_4`, `PlasmidsBCsum`, `plmean`, `log2FoldChange`, `lfcSE` |

### GSM5861614_MPRAduo_counts_WGREST_SKNSH.txt.gz

| Rows | Observed dtype | Columns with this dtype |
| ---: | --- | --- |
| 139,325 | `str` | `id`, `celltype`, `enhancer`, `silencer`, `chromosome`, `project` |
| 139,325 | `float64` | `start`, `end`, `SK.N.SH_1`, `SK.N.SH_2`, `SK.N.SH_3`, `SK.N.SH_4`, `plasmid_1`, `plasmid_2`, `plasmid_3`, `plasmid_4`, `PlasmidsBCsum`, `plmean`, `log2FoldChange`, `lfcSE` |

### GSE196171_nonCan_K562_results.run1.txt.gz

| Rows | Observed dtype | Columns with this dtype |
| ---: | --- | --- |
| 4,975 | `str` | `ID`, `project` |
| 4,975 | `float64` | `SNP`, `chr`, `pos`, `ref_allele`, `alt_allele`, `allele`, `window`, `strand`, `haplotype`, `start`, `stop`, `ctrl_mean`, `exp_mean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj` |

### Important interpretation details

| Field or pattern | Observed meaning / limitation |
| --- | --- |
| WGREST `id` | Unique construct identifier, combining enhancer and silencer context |
| WGREST `silencer` | Candidate element identifier; repeated across five enhancer contexts |
| `project` | Population/category label; exact native labels are `motif` and `noMotif` |
| `celltype` | Source cell label; SK-N-SH is spelled `SK.N.SH` in the file |
| `chromosome`, `start`, `end` | Native interval metadata; source chromosome labels omit the `chr` prefix |
| `log2FoldChange`, `lfcSE` | Processed activity/statistical fields; some values are missing |
| WGREST count columns | Cell-specific replicate columns and plasmid counts; not substitutes for a verified outcome definition |
| Benchmark `start`/`end` | Infer as strings because nonnumeric entries exist; do not coerce silently |
| Noncanonical `chr`, `start`, `stop` | Entirely missing despite the presence of those headers |

The WGREST files contain neither sequence nor strand columns, nor `padj` or a validated binary repression-call column. The noncanonical file does contain `pvalue` and `padj`, but belongs to a different experiment and has no usable interval columns. Those statistics must not be transferred to native elements by assumption.

## Native population and identifier semantics

The following counts apply to each of the four WGREST tables:

| Exact source category | Distinct `silencer` IDs | Construct rows | Complete coordinate rows | Nonmissing `log2FoldChange` rows |
| --- | ---: | ---: | ---: | ---: |
| `motif` | 8,436 | 42,180 | 42,180 | 37,256 |
| `noMotif` | 4,430 | 22,150 | 22,150 | 19,502 |

Every native element has five enhancer contexts: En02, En09, En11, En19, and En21. All complete native intervals have `end - start = 200`. The four cell tables agree on the 12,866 distinct native interval definitions.

These counts do not establish a fully observed modeling matrix. A present element can have a missing measurement in one or more contexts. The shared missing-value totals shown above do not, by themselves, prove that the same measurements are missing in every cell type. Usable K562/SK-N-SH overlap must be computed for an explicit eligibility/context definition.

Do not use a substring match on `motif` to define the canonical native set: categories such as `motif,maf1` and `motif,rare` belong to other source populations. Do not treat construct rows as independent native-element counts.

## Published sequence supplement

| Property | Observed value |
| --- | --- |
| Source | Publisher's Table S2, associated with DOI 10.1016/j.xgen.2022.100234 |
| Local archive | [Table_S2_publisher.zip](raw/assembly_check/Table_S2_publisher.zip) |
| Archive size | 1,470,375 bytes |
| Inspected member | `xgen_100234_mmc03.tsv` |
| Member size | 6,475,283 bytes of uncompressed TSV |
| Columns | `oligo_name`, `sequence`, `library` |
| Total rows | 26,682 |
| `WG_REST` rows | 24,000 |
| `benchmark` rows | 1,687 |
| `non_canonical` rows | 995 |

Oligo names recur across libraries. The assembly audit filters `library == "WG_REST"` before joining names; that subset has unique oligo names. A join across all libraries fails the intended one-to-one relationship and is not an appropriate workaround.

| Native set | GEO elements | Direct names matched to WG_REST oligos | Still unmatched |
| --- | ---: | ---: | ---: |
| Canonical | 8,436 | 5,362 | 3,074 |
| No-motif | 4,430 | 4,430 | 0 |
| Combined | 12,866 | 9,792 | 3,074 |

All 9,792 directly matched sequences are 200 bp. The unmatched canonical names are retained as unresolved identifiers. No sequence was fabricated, no alias was guessed, and no expected population count was lowered to remove them.

## Assembly evidence cache

The audit selects the first directly matched native interval by ascending position within each available project/chromosome group, before querying the references. This yields 24 canonical and 23 no-motif samples, spanning the available chromosome groups.

| Candidate assembly | Exact full-sequence matches | Diagnostic start −1 matches | Diagnostic start +1 matches |
| --- | ---: | ---: | ---: |
| hg19 | 47/47 | 0/47 | 0/47 |
| hg38 | 10/47 | 1/47 | 0/47 |

Every hg19 match is forward-oriented. The result supports hg19 zero-based, half-open native intervals. Ten loci happen to match both builds at the same numbers and do not discriminate the assemblies. The one shifted hg38 match is a diagnostic observation, not a rule to shift every source interval.

Each file under `raw/assembly_check/ucsc/` records a request of the form:

```text
https://api.genome.ucsc.edu/getData/sequence?genome=<assembly>;chrom=<chromosome>;start=<source_start_minus_one>;end=<source_end_plus_one>
```

This is a URL pattern, not a fabricated observed response. Actual URLs are recorded individually in the ledger and comparison CSV. The audit checks the response assembly, chromosome, start, end, and 202-base sequence length against each request before comparing the inner 200 bases. It compares both DNA orientations and diagnostic adjacent windows.

The cache comprises 94 small JSON responses rather than downloaded reference genomes. For the underlying interface, see the [UCSC API coordinate documentation](https://genome.ucsc.edu/goldenPath/help/api.html).

**Scope of the conclusion:** source-assembly identification for the sampled native populations. It does not validate all 12,866 genomic sequences, resolve the 3,074 unmatched canonical IDs, establish variant-coordinate builds, verify motif placement, or perform hg19 → hg38 conversion.

## Acquisition failures and content validation

| Event | Observed response | Resolution / retained evidence |
| --- | --- | --- |
| Initial sandbox access to NCBI/UCSC | Host-resolution errors | Authorized network retries succeeded; errors were reported |
| GEO HTML access through the browsing tool | Browser/reCAPTCHA response | Direct text/directory retrieval supplied the discovery metadata |
| PMC Table S2 endpoint | HTTP response containing 1,817 bytes of HTML, saved as `mmc3.zip` | Retained as failure evidence; not used as sequence input |
| Publisher article page through the browsing tool | HTTP 403 | The separate publisher supplement URL succeeded |
| Publisher Table S2 download | Valid 1,470,375-byte ZIP | Opened, schema inspected, and used in the assembly audit |

The failed PMC response is deliberately distinguishable from the valid archive. Do not overwrite it with a successful download while keeping its original hash record. A future replacement must have its own path and provenance history.

## Pending inputs and configured locations

**Configured paths express expectations, not acquired files.** The following inputs are not present as completed production inputs in the current catalog.

| Input | Intended location or source | Status / prerequisite |
| --- | --- | --- |
| Prepared MPRA table | `data/raw/mpra.tsv` is the current default `MPRA_PATH`; explicit paths can be supplied | Not produced; do not confuse this conventional path with the original GEO tables |
| Full matched element FASTA | `data/processed/elements.fasta` | Not produced |
| REST MEME matrix | `data/raw/MA0138.2.meme` | Not downloaded |
| Native hg19 → hg38 chain | Must be explicitly selected and recorded | Not configured/downloaded by the current target-chain mapping |
| hg38 → rheMac10 | `data/raw/chains/hg38ToRheMac10.over.chain.gz` | Expected path only |
| hg38 → mm39 | `data/raw/chains/hg38ToMm39.over.chain.gz` | Expected path only |
| hg38 → canFam6 | `data/raw/chains/hg38ToCanFam6.over.chain.gz` | Expected path only |
| hg38 → bosTau9 | `data/raw/chains/hg38ToBosTau9.over.chain.gz` | Expected path only |
| hg38 → monDom5 | `data/raw/chains/hg38ToMonDom5.over.chain.gz` | Expected path only |
| Target FASTAs | `data/raw/genomes/{assembly}.fa` | Not downloaded |
| Human hg38 `.2bit` | To be specified in the manual genome acquisition workflow | Not downloaded |
| phyloP / phastCons | `BIGWIG_URLS` currently contains unset entries | Track identity and hg38 compatibility must be verified; use remote byte ranges |
| Table S11 and other required supplements | Primary paper supplement endpoints | Not acquired during the assembly task |
| ENCODE occupancy data | Accessions must be verified from the relevant source tables | No completed occupancy feature input |
| `duoREST` reference checkout | Proposed gitignored `vendor/` directory | Not cloned; only selected remote source files were read |
| Companion GSE212253 data | Original project plan | Not downloaded in the current workflow |

The [Week 1 guide](../data-README.md) defines the remaining acquisition sequence and validation requirements. No genome-download target or six-check validation target is implemented yet.

## Provenance contract

### Required fields and meanings

| Field | Meaning | Important qualification |
| --- | --- | --- |
| File | Repository-relative path to the acquired bytes | Archive members are identified through their parent unless separately materialized |
| Source URL | Exact requested source, including query parameters | A similarly named file at another URL is a different source identity |
| Downloaded UTC | Completion timestamp in UTC | Not the publisher's creation date or a later cache-read time |
| Size | Actual local file length in bytes | Not an approximate human-readable directory size |
| MD5 | Locally computed identity digest | Not evidence of a publisher-checksum match unless explicitly compared |

No upstream MD5 comparison is claimed for the current downloads. Hashes make later file changes detectable against the recorded ledger; they do not independently establish biological validity or publisher authenticity.

### Small downloader behavior

[`download_small_file`](../src/downloads.py) performs no network operation at import time. Its default is cache-only (`fetch=False`), with a 10,000,000-byte cap. An explicit cap may be supplied up to 100,000,000 bytes. It requires HTTPS and a destination below the repository's `data/` directory.

For an existing file, it requires a ledger row matching the path, source URL, and current local MD5. A changed or unrecorded cache raises an error. For a new file, it uses a temporary `.part` path, checks transfer success and size, records a local digest, moves the completed file to its destination, and appends a row internally. It does not provide resumable whole-genome transfers or validate every possible file format; callers must validate the content they requested.

### Maintaining this document

- Preserve every existing acquisition row. Corrections to interpretation belong in explanatory prose or a new record, not altered historical hashes or timestamps.
- Keep original source bytes separate from derived data.
- Record every new actual file/response acquired by the supported downloader.
- Keep a usable ledger table at the end: the downloader appends rows to the end of this file.
- Do not append prose beneath the final open ledger table without first arranging a suitable table for future rows.
- Keep local raw payloads ignored. This README and other documentation can be versioned without adding the data files.
- Use `git status --short`, `git diff --cached --name-only`, and `git check-ignore` before a handoff; do not force-add ignored datasets.

## Historical planning inventory

The original planning tables are retained below so their existing rows remain unchanged. Blank cells are historical placeholders, **not evidence of completed downloads**. The current catalog and observed download ledgers take precedence for acquisition status.

<details>
<summary>Original source, remote-track, and target-genome planning notes</summary>

No data files are committed to this repository. This file is the record of what was
downloaded, from where, and when. Fill a row the moment you download something; a file
whose provenance you cannot reconstruct in December is a file you cannot put in the paper.

### Inputs

| File | Source | URL | Downloaded | Notes |
| --- | --- | --- | --- | --- |
| | GEO GSE196171 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171 | | MPRA per-element activity |
| | GEO GSE212253 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212253 | | Companion BRB-seq |
| | Cell Genomics supplement | https://doi.org/10.1016/j.xgen.2022.100234 | | Table S11 = ENCODE accessions |
| | UCSC | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/ | | hg38 genome |
| | UCSC | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ | | Chain files per species |
| | JASPAR | https://jaspar.elixir.no/matrix/MA0138.2 | | REST profile, MEME format |
| | ENCODE | https://www.encodeproject.org | | REST, RCOR1, SIN3A, EHMT2 in K562 |

### Queried remotely, not downloaded

phyloP and phastCons bigWigs are read over HTTP by byte range rather than downloaded. The
genome-wide files run to several gigabytes each and only the 8,436 motif windows are needed.

| Track | URL | Notes |
| --- | --- | --- |
| phyloP | | Confirm exact file name from the UCSC directory listing |
| phastCons | | |

### Target genomes

| Species | Assembly | Downloaded | Size |
| --- | --- | --- | --- |
| Rhesus macaque | | | |
| Mouse | | | |
| Dog | | | |
| Cow | | | |
| Opossum | | | |

</details>

## Phase A discovery downloads

MD5 values below are computed locally. No upstream MD5 was supplied; these values have not been compared against publisher checksums.

| File | Source URL | Downloaded (UTC) | Size (bytes) | MD5 (local) |
| --- | --- | --- | ---: | --- |
| `data/raw/geo/GSE196171/GSE196171_nonCan_K562_results.run1.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/GSE196171_nonCan_K562_results.run1.txt.gz | 2026-11-03T15:40:30+00:00 | 308452 | `2783c20b489808ce39284646caebc922` |
| `data/raw/geo/GSE196171/filelist.txt` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/filelist.txt | 2026-11-03T15:40:33+00:00 | 596 | `fe8819abe5f7f49ff9848936fe2ed4a3` |
| `data/raw/geo/GSE196171/series_metadata.soft.txt` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171&targ=self&form=text&view=full | 2026-11-03T15:40:35+00:00 | 3685 | `822bab1d23c915e706f9054b1c45b142` |
| `data/raw/geo/GSE196171/supplementary_directory.html` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/ | 2026-11-03T15:40:36+00:00 | 832 | `b0c32e31b6cf636158ea44467ac8202f` |
| `data/raw/geo/GSE196171/GSE196171_RAW.tar` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/GSE196171_RAW.tar | 2026-11-03T15:41:04+00:00 | 31252480 | `6fb66e9dac40962ffc5b8bc955dd96cc` |

## Coordinate assembly investigation downloads

MD5 values are local file identities, not publisher-verified checksums. The first
PMC response below (`mmc3.zip`, 1,817 bytes) is an HTML browser challenge, **not a
valid ZIP**. `Table_S2_publisher.zip` is the successfully opened sequence supplement.
The UCSC JSON files contain only 202 bp query windows, with response coordinates
and assembly checked against the request. All original provenance rows are retained.

| File | Source URL | Downloaded (UTC) | Size (bytes) | MD5 (local) |
| --- | --- | --- | ---: | --- |
| `data/raw/assembly_check/mmc3.zip` | https://pmc.ncbi.nlm.nih.gov/articles/instance/9903721/bin/mmc3.zip | 2026-11-03T15:49:49+00:00 | 1817 | `a74c72d18c52a871d1aeee46d2e31093` |
| `data/raw/assembly_check/Table_S2_publisher.zip` | https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22001926-mmc3.zip | 2026-11-03T15:50:37+00:00 | 1470375 | `74baa0e59cef1f6af284620edd9af9b0` |
| `data/raw/assembly_check/ucsc/hg19_chr1_1235903_1236105.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr1;start=1235903;end=1236105 | 2026-11-03T15:54:13+00:00 | 357 | `f4f351a41094cda0bd3285428cbae351` |
| `data/raw/assembly_check/ucsc/hg38_chr1_1235903_1236105.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=1235903;end=1236105 | 2026-11-03T15:54:14+00:00 | 357 | `95cccebfdc0c3325d4d1564bd4616ed7` |
| `data/raw/assembly_check/ucsc/hg19_chr10_329877_330079.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr10;start=329877;end=330079 | 2026-11-03T15:54:16+00:00 | 356 | `1341e5ca73a385ff91d80f14ca523815` |
| `data/raw/assembly_check/ucsc/hg38_chr10_329877_330079.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=329877;end=330079 | 2026-11-03T15:54:17+00:00 | 356 | `b07eff3ae78cd20183690b765e2a982b` |
| `data/raw/assembly_check/ucsc/hg19_chr11_280224_280426.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr11;start=280224;end=280426 | 2026-11-03T15:54:18+00:00 | 356 | `737a138954096f00ffa5c560d325235b` |
| `data/raw/assembly_check/ucsc/hg38_chr11_280224_280426.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr11;start=280224;end=280426 | 2026-11-03T15:54:20+00:00 | 356 | `4f6c9f546df790d1b876566d72276c77` |
| `data/raw/assembly_check/ucsc/hg19_chr12_176237_176439.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr12;start=176237;end=176439 | 2026-11-03T15:54:21+00:00 | 356 | `22fc0c1aefb4f4ed0ae3bbc336b1de9a` |
| `data/raw/assembly_check/ucsc/hg38_chr12_176237_176439.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr12;start=176237;end=176439 | 2026-11-03T15:54:23+00:00 | 356 | `9872890aa5aebe25c41d85e8b939f394` |
| `data/raw/assembly_check/ucsc/hg19_chr13_21635562_21635764.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr13;start=21635562;end=21635764 | 2026-11-03T15:54:24+00:00 | 360 | `8369a41d60120dd051b9851f8317c4be` |
| `data/raw/assembly_check/ucsc/hg38_chr13_21635562_21635764.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr13;start=21635562;end=21635764 | 2026-11-03T15:54:26+00:00 | 360 | `32e54b0662d2d46f602e86d4d753e56d` |
| `data/raw/assembly_check/ucsc/hg19_chr14_20937583_20937785.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr14;start=20937583;end=20937785 | 2026-11-03T15:54:27+00:00 | 360 | `6fc754fdfc9dba5d6a9ffc62adde5170` |
| `data/raw/assembly_check/ucsc/hg38_chr14_20937583_20937785.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr14;start=20937583;end=20937785 | 2026-11-03T15:54:28+00:00 | 360 | `be0c732c41fc08455729ffb6d9aaa702` |
| `data/raw/assembly_check/ucsc/hg19_chr15_23172924_23173126.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr15;start=23172924;end=23173126 | 2026-11-03T15:54:30+00:00 | 360 | `209228b92acdf7c33445ff93f0921a2e` |
| `data/raw/assembly_check/ucsc/hg38_chr15_23172924_23173126.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr15;start=23172924;end=23173126 | 2026-11-03T15:54:31+00:00 | 360 | `fc9dfb0c6b8cb7870c788ef358142a6f` |
| `data/raw/assembly_check/ucsc/hg19_chr16_750545_750747.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr16;start=750545;end=750747 | 2026-11-03T15:54:33+00:00 | 356 | `956a9f7e4f4659a58a499418fac873b9` |
| `data/raw/assembly_check/ucsc/hg38_chr16_750545_750747.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr16;start=750545;end=750747 | 2026-11-03T15:54:34+00:00 | 356 | `498a1f39918a5b1ca800adae92cb8fbb` |
| `data/raw/assembly_check/ucsc/hg19_chr17_43638_43840.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr17;start=43638;end=43840 | 2026-11-03T15:54:35+00:00 | 354 | `27fc96a83cd84e7853f3b0b7a6534031` |
| `data/raw/assembly_check/ucsc/hg38_chr17_43638_43840.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr17;start=43638;end=43840 | 2026-11-03T15:54:37+00:00 | 354 | `6bb4e154ccfa1bcf3f21f64050efd440` |
| `data/raw/assembly_check/ucsc/hg19_chr18_158525_158727.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr18;start=158525;end=158727 | 2026-11-03T15:54:38+00:00 | 356 | `6620355e431f0be99dcc19ccf8e67b04` |
| `data/raw/assembly_check/ucsc/hg38_chr18_158525_158727.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr18;start=158525;end=158727 | 2026-11-03T15:54:39+00:00 | 356 | `ff199d5f5fed47b3bc120891d3eca688` |
| `data/raw/assembly_check/ucsc/hg19_chr19_482602_482804.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr19;start=482602;end=482804 | 2026-11-03T15:54:41+00:00 | 356 | `08ffd9a55692c1316c31c319fb7764da` |
| `data/raw/assembly_check/ucsc/hg38_chr19_482602_482804.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr19;start=482602;end=482804 | 2026-11-03T15:54:42+00:00 | 356 | `018ab86467ca2b1b6f66dc44fb7831d1` |
| `data/raw/assembly_check/ucsc/hg19_chr2_423528_423730.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr2;start=423528;end=423730 | 2026-11-03T15:54:44+00:00 | 355 | `3a4b0987aafc234842876769856d2cc0` |
| `data/raw/assembly_check/ucsc/hg38_chr2_423528_423730.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr2;start=423528;end=423730 | 2026-11-03T15:54:45+00:00 | 355 | `4b255bc186d8522b1c5b4aa97afca48a` |
| `data/raw/assembly_check/ucsc/hg19_chr20_352447_352649.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr20;start=352447;end=352649 | 2026-11-03T15:54:47+00:00 | 356 | `595a2ad43d4e85bd179bf32911ec43d0` |
| `data/raw/assembly_check/ucsc/hg38_chr20_352447_352649.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr20;start=352447;end=352649 | 2026-11-03T15:54:49+00:00 | 356 | `4c59961f1300d9583d12a06833497445` |
| `data/raw/assembly_check/ucsc/hg19_chr21_15219147_15219349.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr21;start=15219147;end=15219349 | 2026-11-03T15:54:50+00:00 | 360 | `e3abd04283ceeb8c4a7ce20b81d04758` |
| `data/raw/assembly_check/ucsc/hg38_chr21_15219147_15219349.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr21;start=15219147;end=15219349 | 2026-11-03T15:54:51+00:00 | 360 | `8c0356f5ad033b6507098c611fb90367` |
| `data/raw/assembly_check/ucsc/hg19_chr22_17601303_17601505.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr22;start=17601303;end=17601505 | 2026-11-03T15:54:53+00:00 | 360 | `3bc32d82718b395567806cf44346e5ee` |
| `data/raw/assembly_check/ucsc/hg38_chr22_17601303_17601505.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr22;start=17601303;end=17601505 | 2026-11-03T15:54:54+00:00 | 360 | `cf1313ce8ead350569b7e9286d9bb19b` |
| `data/raw/assembly_check/ucsc/hg19_chr3_174151_174353.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr3;start=174151;end=174353 | 2026-11-03T15:54:56+00:00 | 355 | `97f76b5c06cd39736a136e701597b1f4` |
| `data/raw/assembly_check/ucsc/hg38_chr3_174151_174353.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr3;start=174151;end=174353 | 2026-11-03T15:54:57+00:00 | 355 | `90ad635b3bf6acbeef567ce6b43665d5` |
| `data/raw/assembly_check/ucsc/hg19_chr4_569086_569288.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr4;start=569086;end=569288 | 2026-11-03T15:54:59+00:00 | 355 | `bfb58988fa6b3941dadff80545faacaf` |
| `data/raw/assembly_check/ucsc/hg38_chr4_569086_569288.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr4;start=569086;end=569288 | 2026-11-03T15:55:00+00:00 | 355 | `1631399e36fb039dfdb03323654673eb` |
| `data/raw/assembly_check/ucsc/hg19_chr5_57199_57401.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr5;start=57199;end=57401 | 2026-11-03T15:55:01+00:00 | 353 | `2c5e6d1feb0316065d4ef3b45de28a3c` |
| `data/raw/assembly_check/ucsc/hg38_chr5_57199_57401.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr5;start=57199;end=57401 | 2026-11-03T15:55:03+00:00 | 353 | `c5b07c40404c4bd1ff9ed742cff1e475` |
| `data/raw/assembly_check/ucsc/hg19_chr6_1268734_1268936.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr6;start=1268734;end=1268936 | 2026-11-03T15:55:04+00:00 | 357 | `d9b9fab9e576497f64eb25106002f7b2` |
| `data/raw/assembly_check/ucsc/hg38_chr6_1268734_1268936.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr6;start=1268734;end=1268936 | 2026-11-03T15:55:06+00:00 | 357 | `7e6f174489c0ac543da5b372b01af443` |
| `data/raw/assembly_check/ucsc/hg19_chr7_577756_577958.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr7;start=577756;end=577958 | 2026-11-03T15:55:07+00:00 | 355 | `47abd8c4be02b7a8a5f5a9cbe226ec90` |
| `data/raw/assembly_check/ucsc/hg38_chr7_577756_577958.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr7;start=577756;end=577958 | 2026-11-03T15:55:08+00:00 | 355 | `d4fb5adfa1cf2ad6b0a86f44c31cfc80` |
| `data/raw/assembly_check/ucsc/hg19_chr8_529818_530020.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr8;start=529818;end=530020 | 2026-11-03T15:55:10+00:00 | 355 | `96e15715d5a4c0879d9ed44931872d9c` |
| `data/raw/assembly_check/ucsc/hg38_chr8_529818_530020.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr8;start=529818;end=530020 | 2026-11-03T15:55:11+00:00 | 355 | `cc609d4be96ffd86823a9f05f022c6b1` |
| `data/raw/assembly_check/ucsc/hg19_chr9_228804_229006.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr9;start=228804;end=229006 | 2026-11-03T15:55:12+00:00 | 355 | `fe81fd326e9997af4be01c5a877c4912` |
| `data/raw/assembly_check/ucsc/hg38_chr9_228804_229006.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr9;start=228804;end=229006 | 2026-11-03T15:55:14+00:00 | 355 | `455f71d74b443ff24fb4b08633fbf4ea` |
| `data/raw/assembly_check/ucsc/hg19_chrX_2527390_2527592.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chrX;start=2527390;end=2527592 | 2026-11-03T15:55:15+00:00 | 357 | `38b0941fc3559e013945b933eab5061a` |
| `data/raw/assembly_check/ucsc/hg38_chrX_2527390_2527592.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrX;start=2527390;end=2527592 | 2026-11-03T15:55:17+00:00 | 357 | `92218178982e908af6c00cdc9724fd51` |
| `data/raw/assembly_check/ucsc/hg19_chrY_16343681_16343883.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chrY;start=16343681;end=16343883 | 2026-11-03T15:55:18+00:00 | 359 | `ff23a0d1b264eaccb563331c73ddf96a` |
| `data/raw/assembly_check/ucsc/hg38_chrY_16343681_16343883.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrY;start=16343681;end=16343883 | 2026-11-03T15:55:19+00:00 | 359 | `80c0c78ee473ba310202ac2a0b6607ff` |
| `data/raw/assembly_check/ucsc/hg19_chr1_956517_956719.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr1;start=956517;end=956719 | 2026-11-03T15:55:21+00:00 | 355 | `bf9d1285c82531fc2b3eecb9d98be5a9` |
| `data/raw/assembly_check/ucsc/hg38_chr1_956517_956719.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=956517;end=956719 | 2026-11-03T15:55:22+00:00 | 355 | `6d18ae19e4072fe9766242e8b033ad43` |
| `data/raw/assembly_check/ucsc/hg19_chr10_120009_120211.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr10;start=120009;end=120211 | 2026-11-03T15:55:24+00:00 | 356 | `f63da98a930ca415b94a4cb54770adb4` |
| `data/raw/assembly_check/ucsc/hg38_chr10_120009_120211.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=120009;end=120211 | 2026-11-03T15:55:25+00:00 | 356 | `a69ee43d8110da21c82f8c59260e7026` |
| `data/raw/assembly_check/ucsc/hg19_chr11_748184_748386.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr11;start=748184;end=748386 | 2026-11-03T15:55:26+00:00 | 356 | `f44c7dd41073a9c02c00fd023541707d` |
| `data/raw/assembly_check/ucsc/hg38_chr11_748184_748386.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr11;start=748184;end=748386 | 2026-11-03T15:55:28+00:00 | 356 | `c3c83776ce40396e1fb29143cf4530e8` |
| `data/raw/assembly_check/ucsc/hg19_chr12_84663_84865.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr12;start=84663;end=84865 | 2026-11-03T15:55:29+00:00 | 354 | `787dc4926bb0c4056f400b2ccae30b37` |
| `data/raw/assembly_check/ucsc/hg38_chr12_84663_84865.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr12;start=84663;end=84865 | 2026-11-03T15:55:31+00:00 | 354 | `a79589b8ebd1af85a380fbe13550abd7` |
| `data/raw/assembly_check/ucsc/hg19_chr13_21127571_21127773.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr13;start=21127571;end=21127773 | 2026-11-03T15:55:32+00:00 | 360 | `117e73ca21a8ef014bc65cf0c9be27ca` |
| `data/raw/assembly_check/ucsc/hg38_chr13_21127571_21127773.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr13;start=21127571;end=21127773 | 2026-11-03T15:55:33+00:00 | 360 | `924fbcce1cc5e6ea86184c3a602805b6` |
| `data/raw/assembly_check/ucsc/hg19_chr14_21100302_21100504.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr14;start=21100302;end=21100504 | 2026-11-03T15:55:35+00:00 | 360 | `7ed84683f7052dab1b42e5e25719022d` |
| `data/raw/assembly_check/ucsc/hg38_chr14_21100302_21100504.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr14;start=21100302;end=21100504 | 2026-11-03T15:55:37+00:00 | 360 | `5a1a60a4c45003fadbf98a9ce16953a1` |
| `data/raw/assembly_check/ucsc/hg19_chr15_25123368_25123570.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr15;start=25123368;end=25123570 | 2026-11-03T15:55:39+00:00 | 360 | `3942bb2f414c3cdc0f55a4c335044325` |
| `data/raw/assembly_check/ucsc/hg38_chr15_25123368_25123570.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr15;start=25123368;end=25123570 | 2026-11-03T15:55:40+00:00 | 360 | `95ea3189686a5f15be3d0661cdbc2538` |
| `data/raw/assembly_check/ucsc/hg19_chr16_962794_962996.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr16;start=962794;end=962996 | 2026-11-03T15:55:42+00:00 | 356 | `f0e0818769c66cb667c22423af5b4653` |
| `data/raw/assembly_check/ucsc/hg38_chr16_962794_962996.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr16;start=962794;end=962996 | 2026-11-03T15:55:43+00:00 | 356 | `a6efc4aa603682b287b3f2e271928abb` |
| `data/raw/assembly_check/ucsc/hg19_chr17_465146_465348.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr17;start=465146;end=465348 | 2026-11-03T15:55:44+00:00 | 356 | `10e1806359450704e3b526e41085ec1a` |
| `data/raw/assembly_check/ucsc/hg38_chr17_465146_465348.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr17;start=465146;end=465348 | 2026-11-03T15:55:46+00:00 | 356 | `5f979fa7cf3c673bdbe7d68a5c7cd092` |
| `data/raw/assembly_check/ucsc/hg19_chr18_961862_962064.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr18;start=961862;end=962064 | 2026-11-03T15:55:47+00:00 | 356 | `4de65a2bd387077bd114dbb1a80aed97` |
| `data/raw/assembly_check/ucsc/hg38_chr18_961862_962064.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr18;start=961862;end=962064 | 2026-11-03T15:55:48+00:00 | 356 | `907007bbffec9d78a41b7725345dbd9b` |
| `data/raw/assembly_check/ucsc/hg19_chr19_277633_277835.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr19;start=277633;end=277835 | 2026-11-03T15:55:50+00:00 | 356 | `9909be0aa5a0a4a47bce10e9c20abe35` |
| `data/raw/assembly_check/ucsc/hg38_chr19_277633_277835.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr19;start=277633;end=277835 | 2026-11-03T15:55:51+00:00 | 356 | `fa4bbe4b3da9bcd9096e4fdaf0556afb` |
| `data/raw/assembly_check/ucsc/hg19_chr2_529491_529693.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr2;start=529491;end=529693 | 2026-11-03T15:55:53+00:00 | 355 | `b5d88e6c4435f8c80e9a52e07cc3c178` |
| `data/raw/assembly_check/ucsc/hg38_chr2_529491_529693.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr2;start=529491;end=529693 | 2026-11-03T15:55:54+00:00 | 355 | `a4b6a1fe9d33f7e96d3f9cf9929aaea8` |
| `data/raw/assembly_check/ucsc/hg19_chr20_279154_279356.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr20;start=279154;end=279356 | 2026-11-03T15:55:56+00:00 | 356 | `c914c483d376c7e1a84e78767b5a3618` |
| `data/raw/assembly_check/ucsc/hg38_chr20_279154_279356.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr20;start=279154;end=279356 | 2026-11-03T15:55:57+00:00 | 356 | `6fe6f3f0b3c5565483561b3e56dd4db6` |
| `data/raw/assembly_check/ucsc/hg19_chr21_9857745_9857947.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr21;start=9857745;end=9857947 | 2026-11-03T15:55:58+00:00 | 358 | `9e28fe42b2d76df67f5abc876127804a` |
| `data/raw/assembly_check/ucsc/hg38_chr21_9857745_9857947.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr21;start=9857745;end=9857947 | 2026-11-03T15:56:00+00:00 | 358 | `a6ea1107f4cb597f3b9f236fce62327a` |
| `data/raw/assembly_check/ucsc/hg19_chr22_18958110_18958312.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr22;start=18958110;end=18958312 | 2026-11-03T15:56:01+00:00 | 360 | `6506e3543dd75c73de79c91c56dca01c` |
| `data/raw/assembly_check/ucsc/hg38_chr22_18958110_18958312.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr22;start=18958110;end=18958312 | 2026-11-03T15:56:03+00:00 | 360 | `94108fb9e0f3546a73ecfe328be1522d` |
| `data/raw/assembly_check/ucsc/hg19_chr3_2345344_2345546.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr3;start=2345344;end=2345546 | 2026-11-03T15:56:04+00:00 | 357 | `585599f5b491cf027d5fc33f84f9458b` |
| `data/raw/assembly_check/ucsc/hg38_chr3_2345344_2345546.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr3;start=2345344;end=2345546 | 2026-11-03T15:56:05+00:00 | 357 | `b4e9256b04f0194f1aa6c84bfc05bcbc` |
| `data/raw/assembly_check/ucsc/hg19_chr4_467786_467988.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr4;start=467786;end=467988 | 2026-11-03T15:56:07+00:00 | 355 | `869b3ee2f19a10c5bd0dee5525d24384` |
| `data/raw/assembly_check/ucsc/hg38_chr4_467786_467988.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr4;start=467786;end=467988 | 2026-11-03T15:56:08+00:00 | 355 | `1f6d0ef7d8a7d516274ca13138e3c342` |
| `data/raw/assembly_check/ucsc/hg19_chr5_400024_400226.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr5;start=400024;end=400226 | 2026-11-03T15:56:10+00:00 | 355 | `847aa25e950e3604bdedbe23bc6b7c18` |
| `data/raw/assembly_check/ucsc/hg38_chr5_400024_400226.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr5;start=400024;end=400226 | 2026-11-03T15:56:11+00:00 | 355 | `5d321d6ea3d9ee6c06012b3baa6b4f53` |
| `data/raw/assembly_check/ucsc/hg19_chr6_209813_210015.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr6;start=209813;end=210015 | 2026-11-03T15:56:12+00:00 | 355 | `ff51a6001738c691eb59bc66346f8ab2` |
| `data/raw/assembly_check/ucsc/hg38_chr6_209813_210015.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr6;start=209813;end=210015 | 2026-11-03T15:56:14+00:00 | 355 | `eb869744dd451e3b3e9789457f01f473` |
| `data/raw/assembly_check/ucsc/hg19_chr7_3018102_3018304.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr7;start=3018102;end=3018304 | 2026-11-03T15:56:15+00:00 | 357 | `78794806274f8ec25f70042c28d53f5f` |
| `data/raw/assembly_check/ucsc/hg38_chr7_3018102_3018304.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr7;start=3018102;end=3018304 | 2026-11-03T15:56:17+00:00 | 357 | `31a0a55b1532fb7e0f49dd27e504917d` |
| `data/raw/assembly_check/ucsc/hg19_chr8_6522855_6523057.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr8;start=6522855;end=6523057 | 2026-11-03T15:56:18+00:00 | 357 | `e8de9b5afa06775698988b2f65106da4` |
| `data/raw/assembly_check/ucsc/hg38_chr8_6522855_6523057.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr8;start=6522855;end=6523057 | 2026-11-03T15:56:19+00:00 | 357 | `d10267749c8c3255ba487026b180d120` |
| `data/raw/assembly_check/ucsc/hg19_chr9_2015274_2015476.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr9;start=2015274;end=2015476 | 2026-11-03T15:56:21+00:00 | 357 | `140590d554f85bb382c0aa574ea7733a` |
| `data/raw/assembly_check/ucsc/hg38_chr9_2015274_2015476.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr9;start=2015274;end=2015476 | 2026-11-03T15:56:22+00:00 | 357 | `b4721322f31c6fe91650d4f06165317b` |
| `data/raw/assembly_check/ucsc/hg19_chrX_1368280_1368482.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chrX;start=1368280;end=1368482 | 2026-11-03T15:56:24+00:00 | 357 | `3dc40b31c61524c8b5833e4e221963ef` |
| `data/raw/assembly_check/ucsc/hg38_chrX_1368280_1368482.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrX;start=1368280;end=1368482 | 2026-11-03T15:56:25+00:00 | 357 | `00c712c40778ae9ff1f6a74d4759e662` |
| `data/raw/hg19ToHg38.over.chain.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | 2026-11-08T22:14:33+00:00 | 227698 | `35887f73fe5e2231656504d1f6430900` |
| `data/raw/repeats/hg19.rmsk.sql` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.sql | 2026-11-08T22:14:34+00:00 | 2097 | `d379a22c7eacafd27c9c47a29ef45ed9` |
| `data/raw/repeats/hg19.rmsk.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz | 2026-11-08T22:20:37+00:00 | 148359433 | `bfed7e64af68efec90b9d419d98a2379` |
| `data/raw/genomes/hg38.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt | 2026-11-08T22:20:38+00:00 | 720 | `57face159b2845e2f58601d0050db7f7` |
| `data/raw/genomes/hg38.2bit` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit | 2026-11-08T22:25:31+00:00 | 835393456 | `dcc3ea27079aa6dc3f9deccd7275e0f8` |
| `data/raw/genomes/verification/hg38.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes | 2026-11-08T22:38:07+00:00 | 11672 | `c42e9f75fff906c6a7143cb2fac86602` |
| `data/raw/genomes/verification/hg38.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt | 2026-11-08T22:38:07+00:00 | 720 | `57face159b2845e2f58601d0050db7f7` |
| `data/raw/genomes/verification/rheMac10.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.chrom.sizes | 2026-11-08T22:38:08+00:00 | 80592 | `d6e73c0c02e487ced58d5f74ffbe54d6` |
| `data/raw/genomes/verification/rheMac10.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/md5sum.txt | 2026-11-08T22:38:08+00:00 | 524 | `4f5ced12331a3b46b15c78354252815b` |
| `data/raw/genomes/verification/mm39.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes | 2026-11-08T22:38:08+00:00 | 1350 | `967c609a714df35c86a4f52ac85db377` |
| `data/raw/genomes/verification/mm39.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/md5sum.txt | 2026-11-08T22:38:09+00:00 | 657 | `43249e4950c03c6b59bb8ef4d9ec7842` |
| `data/raw/genomes/verification/canFam6.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/canFam6.chrom.sizes | 2026-11-08T22:38:09+00:00 | 3369 | `460096e268ee9f66625e3cd13319bffc` |
| `data/raw/genomes/verification/canFam6.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/md5sum.txt | 2026-11-08T22:38:09+00:00 | 630 | `d3cb79a221b39e348d0e99291bf59e19` |
| `data/raw/genomes/verification/bosTau9.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.chrom.sizes | 2026-11-08T22:38:10+00:00 | 59089 | `6522c148be975984d454d7c47f39f105` |
| `data/raw/genomes/verification/bosTau9.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/md5sum.txt | 2026-11-08T22:38:10+00:00 | 692 | `640c0a34f58434443426a4f0ef518d89` |
| `data/raw/genomes/verification/monDom5.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/monDom5/bigZips/monDom5.chrom.sizes | 2026-11-08T22:38:10+00:00 | 161 | `2f14bb7a9c6f8bd9c62600895288e83e` |
| `data/raw/genomes/verification/monDom5.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/monDom5/bigZips/md5sum.txt | 2026-11-08T22:38:11+00:00 | 361 | `a8acb0462e6dd97ff6c7eef628438aee` |
| `data/raw/tools/twoBitInfo` | https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/twoBitInfo | 2026-11-08T22:38:17+00:00 | 5813896 | `d601d883cbe63fedd2b85d33da2822ae` |
| `data/raw/chains/hg38ToMm39.over.chain.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm39.over.chain.gz | 2026-11-08T22:56:37+00:00 | 73017320 | `5256809040d3ac932a53d16c7b3979b6` |

<!-- appended by scripts/fetch_genomes.sh -->
| rheMac10.2bit | rhesus macaque | https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.2bit | 2026-11-09T00:32:03Z | 0.72 GB | 331f913deb7529334a09ef764cb921c3 |
| mm39.2bit | mouse | https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit | 2026-11-09T00:32:03Z | 0.67 GB | db2e9ace19fb84697ff1c5a4c709599d |
| canFam6.2bit | dog | https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/canFam6.2bit | 2026-11-09T00:32:03Z | 0.57 GB | 1508a6c203ab6493fb4d6f91b2f225cd |
| bosTau9.2bit | cow | https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.2bit | 2026-11-09T00:32:03Z | 0.66 GB | bbd2b89cef5951f39e709c946d8ead46 |
| monDom5.2bit | opossum | https://hgdownload.soe.ucsc.edu/goldenPath/monDom5/bigZips/monDom5.2bit | 2026-11-09T00:32:03Z | 0.89 GB | cd4c3f09156745a3ae4e4150cbaaa332 |
| `data/raw/geo/GSE196171/GSE196171_RAW.tar` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/GSE196171_RAW.tar | 2026-11-10T21:07:41+00:00 | 31252480 | `6fb66e9dac40962ffc5b8bc955dd96cc` |
| `data/raw/assembly_check/Table_S2_publisher.zip` | https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22001926-mmc3.zip | 2026-11-10T21:08:59+00:00 | 1470375 | `74baa0e59cef1f6af284620edd9af9b0` |
| `data/raw/assembly_check/ucsc/hg19_chr1_1235903_1236105.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr1;start=1235903;end=1236105 | 2026-11-10T21:09:02+00:00 | 357 | `6ead0e4fa05a73461ccc2d1e31960456` |
| `data/raw/assembly_check/ucsc/hg38_chr1_1235903_1236105.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=1235903;end=1236105 | 2026-11-10T21:09:04+00:00 | 357 | `5392c54a2a1f088e1c987d98272f19a6` |
| `data/raw/assembly_check/ucsc/hg19_chr10_329877_330079.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr10;start=329877;end=330079 | 2026-11-10T21:09:05+00:00 | 356 | `f285619c68bc0b41a42a2dce2b10035e` |
| `data/raw/assembly_check/ucsc/hg38_chr10_329877_330079.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=329877;end=330079 | 2026-11-10T21:09:06+00:00 | 356 | `e9b3025a5fe0a7836e0c8da528cb0620` |
| `data/raw/assembly_check/ucsc/hg19_chr11_280224_280426.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr11;start=280224;end=280426 | 2026-11-10T21:09:08+00:00 | 356 | `34ea751de88f84e57a12be73a023ca5f` |
| `data/raw/assembly_check/ucsc/hg38_chr11_280224_280426.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr11;start=280224;end=280426 | 2026-11-10T21:09:09+00:00 | 356 | `5a40aa4fec8257fadde45887e55fe4fd` |
| `data/raw/assembly_check/ucsc/hg19_chr12_176237_176439.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr12;start=176237;end=176439 | 2026-11-10T21:09:11+00:00 | 356 | `f6b3e518b63c963feb4d9896374d1e49` |
| `data/raw/assembly_check/ucsc/hg38_chr12_176237_176439.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr12;start=176237;end=176439 | 2026-11-10T21:09:12+00:00 | 356 | `a1f0eff164180713440c5a3078f8a91a` |
| `data/raw/assembly_check/ucsc/hg19_chr13_21635562_21635764.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr13;start=21635562;end=21635764 | 2026-11-10T21:09:13+00:00 | 360 | `449b5b723368b934863c9eefe31f1511` |
| `data/raw/assembly_check/ucsc/hg38_chr13_21635562_21635764.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr13;start=21635562;end=21635764 | 2026-11-10T21:09:15+00:00 | 360 | `c3a1b03fa8f4a7997c79e392f3fb9b92` |
| `data/raw/assembly_check/ucsc/hg19_chr14_20937583_20937785.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr14;start=20937583;end=20937785 | 2026-11-10T21:09:16+00:00 | 360 | `10f09d7f6e67dfe048adbc8370c0dda1` |
| `data/raw/assembly_check/ucsc/hg38_chr14_20937583_20937785.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr14;start=20937583;end=20937785 | 2026-11-10T21:09:17+00:00 | 360 | `013ebe7799cfcc649c936a878c4eac4f` |
| `data/raw/assembly_check/ucsc/hg19_chr15_23172924_23173126.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr15;start=23172924;end=23173126 | 2026-11-10T21:09:19+00:00 | 360 | `d1d587a9038a4f89a9fdabcdabdd3c46` |
| `data/raw/assembly_check/ucsc/hg38_chr15_23172924_23173126.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr15;start=23172924;end=23173126 | 2026-11-10T21:09:20+00:00 | 360 | `7fd86d037ffcc775c0b865144a7f27b2` |
| `data/raw/assembly_check/ucsc/hg19_chr16_750545_750747.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr16;start=750545;end=750747 | 2026-11-10T21:09:21+00:00 | 356 | `6d995f07dbb462d93777bcc66b39aa82` |
| `data/raw/assembly_check/ucsc/hg38_chr16_750545_750747.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr16;start=750545;end=750747 | 2026-11-10T21:09:23+00:00 | 356 | `b7f912bdd6ac748a483b9ece7099d9c1` |
| `data/raw/assembly_check/ucsc/hg19_chr17_43638_43840.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr17;start=43638;end=43840 | 2026-11-10T21:09:24+00:00 | 354 | `7abd6151ddca1e60f9e4b91b3d363ba8` |
| `data/raw/assembly_check/ucsc/hg38_chr17_43638_43840.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr17;start=43638;end=43840 | 2026-11-10T21:09:26+00:00 | 354 | `d9b14335d88b73c6f95f76fda01b3ec3` |
| `data/raw/assembly_check/ucsc/hg19_chr18_158525_158727.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr18;start=158525;end=158727 | 2026-11-10T21:09:27+00:00 | 356 | `3aeb94481a6b5f4d32840b4d138ec684` |
| `data/raw/assembly_check/ucsc/hg38_chr18_158525_158727.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr18;start=158525;end=158727 | 2026-11-10T21:09:28+00:00 | 356 | `cc24643747895be843729865467dc249` |
| `data/raw/assembly_check/ucsc/hg19_chr19_482602_482804.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr19;start=482602;end=482804 | 2026-11-10T21:09:31+00:00 | 356 | `53f2ef0fa71102230d1ccbaa0b882427` |
| `data/raw/assembly_check/ucsc/hg38_chr19_482602_482804.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr19;start=482602;end=482804 | 2026-11-10T21:09:32+00:00 | 356 | `556d5f3133284f5cbee50efdd38e1aa2` |
| `data/raw/assembly_check/ucsc/hg19_chr2_423528_423730.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr2;start=423528;end=423730 | 2026-11-10T21:09:33+00:00 | 355 | `65437e2054a6b67aae94a465eb451ea1` |
| `data/raw/assembly_check/ucsc/hg38_chr2_423528_423730.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr2;start=423528;end=423730 | 2026-11-10T21:09:35+00:00 | 355 | `8863124d0ed3bc9048f28ad2af274c03` |
| `data/raw/assembly_check/ucsc/hg19_chr20_352447_352649.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr20;start=352447;end=352649 | 2026-11-10T21:09:36+00:00 | 356 | `cb7bd953cfd41e0e831c3e9654d9d989` |
| `data/raw/assembly_check/ucsc/hg38_chr20_352447_352649.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr20;start=352447;end=352649 | 2026-11-10T21:09:38+00:00 | 356 | `3121ad8210d8d18270692f0eeac3edb2` |
| `data/raw/assembly_check/ucsc/hg19_chr21_15219147_15219349.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr21;start=15219147;end=15219349 | 2026-11-10T21:09:39+00:00 | 360 | `4078c0c5242a87e2a0f8a14891d5e2eb` |
| `data/raw/assembly_check/ucsc/hg38_chr21_15219147_15219349.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr21;start=15219147;end=15219349 | 2026-11-10T21:09:41+00:00 | 360 | `c6e268b2ae87ef207d61bce475b1136b` |
| `data/raw/assembly_check/ucsc/hg19_chr22_17601303_17601505.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr22;start=17601303;end=17601505 | 2026-11-10T21:09:42+00:00 | 360 | `6d12b14766e91f02d4331d806982ffe2` |
| `data/raw/assembly_check/ucsc/hg38_chr22_17601303_17601505.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr22;start=17601303;end=17601505 | 2026-11-10T21:09:43+00:00 | 360 | `a35ce2db2902587d6ae0638e228e7956` |
| `data/raw/assembly_check/ucsc/hg19_chr3_174151_174353.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr3;start=174151;end=174353 | 2026-11-10T21:09:45+00:00 | 355 | `15342b22246e9b58c6ba0820209d1de5` |
| `data/raw/assembly_check/ucsc/hg38_chr3_174151_174353.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr3;start=174151;end=174353 | 2026-11-10T21:09:46+00:00 | 355 | `224ad0f0ede329b7e60d0b6edf37593f` |
| `data/raw/assembly_check/ucsc/hg19_chr4_569086_569288.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr4;start=569086;end=569288 | 2026-11-10T21:09:47+00:00 | 355 | `305f718936aef548dd311fa97bd8eff3` |
| `data/raw/assembly_check/ucsc/hg38_chr4_569086_569288.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr4;start=569086;end=569288 | 2026-11-10T21:09:49+00:00 | 355 | `40640baa1b396552337f51aec290a725` |
| `data/raw/assembly_check/ucsc/hg19_chr5_57199_57401.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr5;start=57199;end=57401 | 2026-11-10T21:09:51+00:00 | 353 | `7504a16c87669c9f40fb5e54fa9407a9` |
| `data/raw/assembly_check/ucsc/hg38_chr5_57199_57401.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr5;start=57199;end=57401 | 2026-11-10T21:09:52+00:00 | 353 | `623ab4fbc3905e1e805fba0d5328dfb9` |
| `data/raw/assembly_check/ucsc/hg19_chr6_1268734_1268936.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr6;start=1268734;end=1268936 | 2026-11-10T21:09:53+00:00 | 357 | `011ca05e01ee24979116aa8ac71d1e94` |
| `data/raw/assembly_check/ucsc/hg38_chr6_1268734_1268936.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr6;start=1268734;end=1268936 | 2026-11-10T21:09:55+00:00 | 357 | `1f4c190247dedbbda4751ad088c87766` |
| `data/raw/assembly_check/ucsc/hg19_chr7_577756_577958.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr7;start=577756;end=577958 | 2026-11-10T21:09:56+00:00 | 355 | `6dadf23f9bbd3b991a2cf086bd679274` |
| `data/raw/assembly_check/ucsc/hg38_chr7_577756_577958.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr7;start=577756;end=577958 | 2026-11-10T21:09:57+00:00 | 355 | `01a14c5f2d3f40dcb9a086e6f5ca61f3` |
| `data/raw/assembly_check/ucsc/hg19_chr8_529818_530020.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr8;start=529818;end=530020 | 2026-11-10T21:09:59+00:00 | 355 | `1f6c885159017f85a8f2a53bb53c9b5d` |
| `data/raw/assembly_check/ucsc/hg38_chr8_529818_530020.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr8;start=529818;end=530020 | 2026-11-10T21:10:01+00:00 | 355 | `5f2a54581aac91826f1e2642f32daee4` |
| `data/raw/assembly_check/ucsc/hg19_chr9_228804_229006.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr9;start=228804;end=229006 | 2026-11-10T21:10:02+00:00 | 355 | `6c1809c590e9862f0924f895eb7a6e3e` |
| `data/raw/assembly_check/ucsc/hg38_chr9_228804_229006.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr9;start=228804;end=229006 | 2026-11-10T21:10:03+00:00 | 355 | `41e89f13cfb492c0815e36c6d9d82dfc` |
| `data/raw/assembly_check/ucsc/hg19_chrX_2527390_2527592.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chrX;start=2527390;end=2527592 | 2026-11-10T21:10:05+00:00 | 357 | `54c2726b6043f2c43ee321fe95bf16e9` |
| `data/raw/assembly_check/ucsc/hg38_chrX_2527390_2527592.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrX;start=2527390;end=2527592 | 2026-11-10T21:10:06+00:00 | 357 | `6f6b3c68b42b3ef96bdc88edbc7c01d0` |
| `data/raw/assembly_check/ucsc/hg19_chrY_16343681_16343883.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chrY;start=16343681;end=16343883 | 2026-11-10T21:10:07+00:00 | 359 | `c6253438c5ca77a326d4a3bea130f885` |
| `data/raw/assembly_check/ucsc/hg38_chrY_16343681_16343883.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrY;start=16343681;end=16343883 | 2026-11-10T21:10:09+00:00 | 359 | `c0221365591cdc3f508941a30451325c` |
| `data/raw/assembly_check/ucsc/hg19_chr1_956517_956719.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr1;start=956517;end=956719 | 2026-11-10T21:10:10+00:00 | 355 | `183f0c82052ba7004c30af8c2501d533` |
| `data/raw/assembly_check/ucsc/hg38_chr1_956517_956719.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=956517;end=956719 | 2026-11-10T21:10:11+00:00 | 355 | `bc4aad6c3a9611eadd6f4b26cdc4ce72` |
| `data/raw/assembly_check/ucsc/hg19_chr10_120009_120211.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr10;start=120009;end=120211 | 2026-11-10T21:10:13+00:00 | 356 | `077c1b001a5a69f046e1a778a7c4d02d` |
| `data/raw/assembly_check/ucsc/hg38_chr10_120009_120211.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=120009;end=120211 | 2026-11-10T21:10:14+00:00 | 356 | `3633a2263581dbd5ddb0a8aec82086f6` |
| `data/raw/assembly_check/ucsc/hg19_chr11_748184_748386.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr11;start=748184;end=748386 | 2026-11-10T21:10:16+00:00 | 356 | `bfe3585f44a4861af6f4dd236e3ba688` |
| `data/raw/assembly_check/ucsc/hg38_chr11_748184_748386.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr11;start=748184;end=748386 | 2026-11-10T21:10:17+00:00 | 356 | `55fdf72943105fd4f7a597719194492c` |
| `data/raw/assembly_check/ucsc/hg19_chr12_84663_84865.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr12;start=84663;end=84865 | 2026-11-10T21:10:18+00:00 | 354 | `01034282ba31987cc1490fd8769523ee` |
| `data/raw/assembly_check/ucsc/hg38_chr12_84663_84865.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr12;start=84663;end=84865 | 2026-11-10T21:10:20+00:00 | 354 | `309e02a7e1a73ec81e0afad3bede974c` |
| `data/raw/assembly_check/ucsc/hg19_chr13_21127571_21127773.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr13;start=21127571;end=21127773 | 2026-11-10T21:10:21+00:00 | 360 | `c70cdd820ebb901acea53f41e9f2c0cd` |
| `data/raw/assembly_check/ucsc/hg38_chr13_21127571_21127773.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr13;start=21127571;end=21127773 | 2026-11-10T21:10:22+00:00 | 360 | `fdcaf9b40ebcb6c751ee5e97ffdc0834` |
| `data/raw/assembly_check/ucsc/hg19_chr14_21100302_21100504.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr14;start=21100302;end=21100504 | 2026-11-10T21:10:24+00:00 | 360 | `0b2ab728e8e793819a5b3196d19610e0` |
| `data/raw/assembly_check/ucsc/hg38_chr14_21100302_21100504.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr14;start=21100302;end=21100504 | 2026-11-10T21:10:25+00:00 | 360 | `1ff3d40e8247fd3872e1ed2f4601b712` |
| `data/raw/assembly_check/ucsc/hg19_chr15_25123368_25123570.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr15;start=25123368;end=25123570 | 2026-11-10T21:10:27+00:00 | 360 | `52d59d8c19dec9181de7a70c86ded0b1` |
| `data/raw/assembly_check/ucsc/hg38_chr15_25123368_25123570.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr15;start=25123368;end=25123570 | 2026-11-10T21:10:28+00:00 | 360 | `3e6f0564a3868268aab37958e4536ec1` |
| `data/raw/assembly_check/ucsc/hg19_chr16_962794_962996.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr16;start=962794;end=962996 | 2026-11-10T21:10:30+00:00 | 356 | `08dc79fc4a3760d2954aa53e2f16a901` |
| `data/raw/assembly_check/ucsc/hg38_chr16_962794_962996.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr16;start=962794;end=962996 | 2026-11-10T21:10:31+00:00 | 356 | `32eeb2a53168ecc64a09f0998cb12b2f` |
| `data/raw/assembly_check/ucsc/hg19_chr17_465146_465348.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr17;start=465146;end=465348 | 2026-11-10T21:10:34+00:00 | 356 | `1664fa0829bfdb392a3d36144b062ae1` |
| `data/raw/assembly_check/ucsc/hg38_chr17_465146_465348.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr17;start=465146;end=465348 | 2026-11-10T21:10:35+00:00 | 356 | `4accfdf3304ed362edde417a601bccf3` |
| `data/raw/assembly_check/ucsc/hg19_chr18_961862_962064.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr18;start=961862;end=962064 | 2026-11-10T21:10:37+00:00 | 356 | `5f941a9b7d8b2bedb77686cb239cd4db` |
| `data/raw/assembly_check/ucsc/hg38_chr18_961862_962064.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr18;start=961862;end=962064 | 2026-11-10T21:10:38+00:00 | 356 | `b3e014e0c43a1f729e4f0ba4d3a209ff` |
| `data/raw/assembly_check/ucsc/hg19_chr19_277633_277835.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr19;start=277633;end=277835 | 2026-11-10T21:10:40+00:00 | 356 | `874bced0cdeb31176398dfcfb9be2c33` |
| `data/raw/assembly_check/ucsc/hg38_chr19_277633_277835.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr19;start=277633;end=277835 | 2026-11-10T21:10:41+00:00 | 356 | `fbdff10799f3b77cfcd7bcadd8235ce1` |
| `data/raw/assembly_check/ucsc/hg19_chr2_529491_529693.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr2;start=529491;end=529693 | 2026-11-10T21:10:43+00:00 | 355 | `769c2e6abd2a1767bb699f1768ddd6c0` |
| `data/raw/assembly_check/ucsc/hg38_chr2_529491_529693.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr2;start=529491;end=529693 | 2026-11-10T21:10:44+00:00 | 355 | `62cda39a914d41e1f19fcd93429d90bd` |
| `data/raw/assembly_check/ucsc/hg19_chr20_279154_279356.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr20;start=279154;end=279356 | 2026-11-10T21:10:46+00:00 | 356 | `665c1993e107ffb2084f7b4fa7ccdf25` |
| `data/raw/assembly_check/ucsc/hg38_chr20_279154_279356.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr20;start=279154;end=279356 | 2026-11-10T21:10:47+00:00 | 356 | `3307b55dd142ac9e3da568b279429718` |
| `data/raw/assembly_check/ucsc/hg19_chr21_9857745_9857947.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr21;start=9857745;end=9857947 | 2026-11-10T21:10:49+00:00 | 358 | `df344c10fd5e4b0f7238d0b3c0e8e1d7` |
| `data/raw/assembly_check/ucsc/hg38_chr21_9857745_9857947.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr21;start=9857745;end=9857947 | 2026-11-10T21:10:50+00:00 | 358 | `188b8b6c8598e398ed56ee90cdeaf267` |
| `data/raw/assembly_check/ucsc/hg19_chr22_18958110_18958312.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr22;start=18958110;end=18958312 | 2026-11-10T21:10:52+00:00 | 360 | `64a53c428b42cf0eb9b5c68c5b3dcda1` |
| `data/raw/assembly_check/ucsc/hg38_chr22_18958110_18958312.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr22;start=18958110;end=18958312 | 2026-11-10T21:10:53+00:00 | 360 | `f848d92eef9119430d5d831017f7b9f2` |
| `data/raw/assembly_check/ucsc/hg19_chr3_2345344_2345546.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr3;start=2345344;end=2345546 | 2026-11-10T21:10:54+00:00 | 357 | `7df78c1facfd36735a6253c247d98236` |
| `data/raw/assembly_check/ucsc/hg38_chr3_2345344_2345546.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr3;start=2345344;end=2345546 | 2026-11-10T21:10:56+00:00 | 357 | `0a10edf0cf3b780862ad8cee71c6c3dc` |
| `data/raw/assembly_check/ucsc/hg19_chr4_467786_467988.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr4;start=467786;end=467988 | 2026-11-10T21:10:57+00:00 | 355 | `33cd6d93a9e77255c2a8ab9454a088e1` |
| `data/raw/assembly_check/ucsc/hg38_chr4_467786_467988.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr4;start=467786;end=467988 | 2026-11-10T21:10:59+00:00 | 355 | `ef0a80aab639f01b01729751fb64aa4d` |
| `data/raw/assembly_check/ucsc/hg19_chr5_400024_400226.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr5;start=400024;end=400226 | 2026-11-10T21:11:00+00:00 | 355 | `fa3395425d1b70feda10625fc6327dd5` |
| `data/raw/assembly_check/ucsc/hg38_chr5_400024_400226.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr5;start=400024;end=400226 | 2026-11-10T21:11:02+00:00 | 355 | `fb1e53eb61fd2ec5a75edc89af3b5f2f` |
| `data/raw/assembly_check/ucsc/hg19_chr6_209813_210015.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr6;start=209813;end=210015 | 2026-11-10T21:11:04+00:00 | 355 | `f7b3b6d2f2bd34e19c680888cbb16bdf` |
| `data/raw/assembly_check/ucsc/hg38_chr6_209813_210015.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr6;start=209813;end=210015 | 2026-11-10T21:11:06+00:00 | 355 | `7b0ca45f085a019551e56a1f297d8bc4` |
| `data/raw/assembly_check/ucsc/hg19_chr7_3018102_3018304.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr7;start=3018102;end=3018304 | 2026-11-10T21:11:08+00:00 | 357 | `468f45264c1077db710a5633705009ee` |
| `data/raw/assembly_check/ucsc/hg38_chr7_3018102_3018304.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr7;start=3018102;end=3018304 | 2026-11-10T21:11:10+00:00 | 357 | `b3f02ebe51779d43aa660b3f8bb7f53e` |
| `data/raw/assembly_check/ucsc/hg19_chr8_6522855_6523057.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr8;start=6522855;end=6523057 | 2026-11-10T21:11:11+00:00 | 357 | `95322f174b2d6b9894820418e81ea733` |
| `data/raw/assembly_check/ucsc/hg38_chr8_6522855_6523057.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr8;start=6522855;end=6523057 | 2026-11-10T21:11:12+00:00 | 357 | `83d5ed1746f67f97a29f83dcd66b2690` |
| `data/raw/assembly_check/ucsc/hg19_chr9_2015274_2015476.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr9;start=2015274;end=2015476 | 2026-11-10T21:11:14+00:00 | 357 | `1005d98c76ef0f0a6cd614506a3acf17` |
| `data/raw/assembly_check/ucsc/hg38_chr9_2015274_2015476.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr9;start=2015274;end=2015476 | 2026-11-10T21:11:15+00:00 | 357 | `ae282a508b1b06d2ee84d0c59257a6f7` |
| `data/raw/assembly_check/ucsc/hg19_chrX_1368280_1368482.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chrX;start=1368280;end=1368482 | 2026-11-10T21:11:17+00:00 | 357 | `c70c28a1c15f70aa449120d7ab23bdc4` |
| `data/raw/assembly_check/ucsc/hg38_chrX_1368280_1368482.json` | https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrX;start=1368280;end=1368482 | 2026-11-10T21:11:18+00:00 | 357 | `61c4cd7b3dc0014c1a1fdc09b666c557` |
| `data/raw/hg19ToHg38.over.chain.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | 2026-11-10T21:12:49+00:00 | 227698 | `35887f73fe5e2231656504d1f6430900` |
| `data/raw/repeats/hg19.rmsk.sql` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.sql | 2026-11-10T21:12:50+00:00 | 2097 | `d379a22c7eacafd27c9c47a29ef45ed9` |
| `data/raw/repeats/hg19.rmsk.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz | 2026-11-10T21:14:29+00:00 | 148359433 | `bfed7e64af68efec90b9d419d98a2379` |
| `data/raw/genomes/hg38.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt | 2026-11-10T21:14:30+00:00 | 720 | `57face159b2845e2f58601d0050db7f7` |
| `data/raw/genomes/hg38.2bit` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit | 2026-11-10T21:47:41+00:00 | 835393456 | `dcc3ea27079aa6dc3f9deccd7275e0f8` |
| `data/raw/MA0138.2.meme` | https://jaspar.elixir.no/api/v1/matrix/MA0138.2.meme | 2026-11-10T22:12:40+00:00 | 1073 | `4de4483fff043233a0f077d2731d7275` |
| `data/raw/genomes/verification/hg38.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes | 2026-11-10T22:12:42+00:00 | 11672 | `c42e9f75fff906c6a7143cb2fac86602` |
| `data/raw/genomes/verification/hg38.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt | 2026-11-10T22:12:42+00:00 | 720 | `57face159b2845e2f58601d0050db7f7` |
| `data/raw/genomes/verification/mm39.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes | 2026-11-10T22:12:42+00:00 | 1350 | `967c609a714df35c86a4f52ac85db377` |
| `data/raw/genomes/verification/mm39.md5sum.txt` | https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/md5sum.txt | 2026-11-10T22:12:43+00:00 | 657 | `43249e4950c03c6b59bb8ef4d9ec7842` |

<!-- appended by scripts/fetch_genomes.sh -->
| mm39.2bit | mouse | https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit | 2026-11-10T23:10:18Z | 0.67 GB | db2e9ace19fb84697ff1c5a4c709599d |
