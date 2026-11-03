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
