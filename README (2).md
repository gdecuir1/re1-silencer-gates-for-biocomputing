# RE1 Silencer Gates

**Mining the human genome for gene-regulatory parts to control neural biocomputing substrates.**

Status: in progress. Johns Hopkins, Applied Comparative Genomics (Fall 2026).

---

## What this is

Companies like Cortical Labs and FinalSpark now build computers that run on living human
neurons. The hardware works; the tissue is the bottleneck. Tuning those neurons means
controlling neuronal gene expression, and the parts available for that are mostly borrowed
from bacteria and viruses.

The human genome has a native option. REST (also called NRSF) is the master repressor of
neuronal genes: it binds a 21 bp motif called RE1 and shuts off nearby transcription. That
is functionally a NOT gate, with REST level as the input and the RE1-controlled gene as the
output.

The problem is that human RE1 sites are not interchangeable. This project asks whether
evolutionary age predicts how well a site would work as a part, and produces a ranked list
of candidates scored the way circuit engineers score gates.

## The question

> Do deeply conserved RE1 elements give stronger, less leaky, and more variant-tolerant
> repression than lineage-specific ones, and can conservation be used to select RE1 parts
> for neural circuits?

About a third of human RE1 elements are primate-specific and bind REST more weakly than
ancient, deeply conserved sites. If evolutionary history carries information that motif
strength alone does not, conservation becomes a cheap selection criterion for circuit parts.

## Approach

| Stage | What happens |
| --- | --- |
| Data | 8,436 RE1 elements with genome-wide MPRA measurements of repression strength, across four cell types |
| Classification | Each element lifted through UCSC hg38 chain files to five mammals; the REST motif rescored in each with JASPAR MA0138 |
| Test | Does conservation class predict repression once motif score, basal activity, TE overlap, and GC are controlled for? |
| Models | Motif score alone as the baseline, against motif + conservation, motif + cofactor ChIP-seq occupancy, and a gapped k-mer SVM on the 200 bp element |
| Validation | Held-out chromosomes, with sequence-similarity clustering so TE-derived near-duplicates cannot leak across folds |
| Output | RE1 candidates ranked on dynamic range, off-state leakiness, and variant robustness |

A neural-context comparison runs alongside: the MPRA includes SK-N-SH, a neuroblastoma line,
so conservation effects can be compared between a neural-lineage cell and K562.

## Results

In progress. This section will carry the headline figure and the model comparison table.

## Reproducing

```bash
git clone https://github.com/<user>/re1-silencer-gates.git
cd re1-silencer-gates
conda env create -f environment.yml
conda activate re1-gates
```

Input data is not committed; `data/README.md` lists every file with its source URL and
download date. Fold assignments are committed so model comparisons are reproducible.

## Data sources

| Source | Use |
| --- | --- |
| [GEO GSE196171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196171) | MPRA repression measurements and variant effects |
| [ENCODE](https://www.encodeproject.org) | REST, RCOR1, SIN3A, EHMT2 ChIP-seq |
| [UCSC hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/) | Genome, liftOver chains, phyloP and phastCons |
| [JASPAR MA0138.2](https://jaspar.elixir.no/matrix/MA0138.2) | REST binding profile |

## Key references

1. Mouri K, et al. Whole-genome functional characterization of RE1 silencers using a modified
   massively parallel reporter assay. *Cell Genomics* 2023;3(1):100234.
   [doi:10.1016/j.xgen.2022.100234](https://doi.org/10.1016/j.xgen.2022.100234)
2. Johnson R, et al. Evolution of the vertebrate gene regulatory network controlled by the
   transcriptional repressor REST. *Molecular Biology and Evolution* 2009;26(7):1491-1507.
   [doi:10.1093/molbev/msp058](https://doi.org/10.1093/molbev/msp058)
3. Kagan BJ, et al. In vitro neurons learn and exhibit sentience when embodied in a simulated
   game-world. *Neuron* 2022;110(23):3952-3969.
   [doi:10.1016/j.neuron.2022.09.001](https://doi.org/10.1016/j.neuron.2022.09.001)
4. Smirnova L, et al. Organoid intelligence (OI): the new frontier in biocomputing and
   intelligence-in-a-dish. *Frontiers in Science* 2023;1:1017235.
   [doi:10.3389/fsci.2023.1017235](https://doi.org/10.3389/fsci.2023.1017235)

## Scope

This work supplies the gene-regulatory control layer for neural biocomputing, not the
computation itself. The underlying measurements come from episomal reporter assays in cell
lines, so the output is a prioritized candidate list for experimental testing rather than
demonstrated gate behavior in neurons.

## Team

Gabrielle DeCuir and Vrinda Sehgal
