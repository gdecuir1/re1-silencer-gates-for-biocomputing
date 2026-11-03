# Data inventory

No data files are committed to this repository. This file is the record of what was
downloaded, from where, and when. Fill a row the moment you download something; a file
whose provenance you cannot reconstruct in December is a file you cannot put in the paper.

## Inputs

| File | Source | URL | Downloaded | Notes |
| --- | --- | --- | --- | --- |
| | GEO GSE196171 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171 | | MPRA per-element activity |
| | GEO GSE212253 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212253 | | Companion BRB-seq |
| | Cell Genomics supplement | https://doi.org/10.1016/j.xgen.2022.100234 | | Table S11 = ENCODE accessions |
| | UCSC | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/ | | hg38 genome |
| | UCSC | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ | | Chain files per species |
| | JASPAR | https://jaspar.elixir.no/matrix/MA0138.2 | | REST profile, MEME format |
| | ENCODE | https://www.encodeproject.org | | REST, RCOR1, SIN3A, EHMT2 in K562 |

## Queried remotely, not downloaded

phyloP and phastCons bigWigs are read over HTTP by byte range rather than downloaded. The
genome-wide files run to several gigabytes each and only the 8,436 motif windows are needed.

| Track | URL | Notes |
| --- | --- | --- |
| phyloP | | Confirm exact file name from the UCSC directory listing |
| phastCons | | |

## Target genomes

| Species | Assembly | Downloaded | Size |
| --- | --- | --- | --- |
| Rhesus macaque | | | |
| Mouse | | | |
| Dog | | | |
| Cow | | | |
| Opossum | | | |

## Phase A discovery downloads

MD5 values below are computed locally. No upstream MD5 was supplied; these values have not been compared against publisher checksums.

| File | Source URL | Downloaded (UTC) | Size (bytes) | MD5 (local) |
| --- | --- | --- | ---: | --- |
| `data/raw/geo/GSE196171/GSE196171_nonCan_K562_results.run1.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/GSE196171_nonCan_K562_results.run1.txt.gz | 2026-11-03T15:40:30+00:00 | 308452 | `2783c20b489808ce39284646caebc922` |
| `data/raw/geo/GSE196171/filelist.txt` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/filelist.txt | 2026-11-03T15:40:33+00:00 | 596 | `fe8819abe5f7f49ff9848936fe2ed4a3` |
| `data/raw/geo/GSE196171/series_metadata.soft.txt` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171&targ=self&form=text&view=full | 2026-11-03T15:40:35+00:00 | 3685 | `822bab1d23c915e706f9054b1c45b142` |
| `data/raw/geo/GSE196171/supplementary_directory.html` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/ | 2026-11-03T15:40:36+00:00 | 832 | `b0c32e31b6cf636158ea44467ac8202f` |
| `data/raw/geo/GSE196171/GSE196171_RAW.tar` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196171/suppl/GSE196171_RAW.tar | 2026-11-03T15:41:04+00:00 | 31252480 | `6fb66e9dac40962ffc5b8bc955dd96cc` |
