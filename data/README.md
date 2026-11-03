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
