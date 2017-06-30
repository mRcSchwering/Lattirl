# Lattirl

In experimental design of single cell RNA sequencing (scRNA-seq) experiments the biological question of interest
is often confounded by batches (e.g. measure biol. group I in batch 1, then group II in batch 2).
When treating measurements from different batches as equally (i.e. pooling batches) for differential expression analysis,
then such methods lose their control over the false discovery rate (FDR).

Usually such methods completely underestimate the true FDR.
The whole problem is further discussed [here](http://b210-research.dkfz.de/computational-genome-biology/scRNAseq/).
This packages provides a model for predicting the loss of FDR control for scDD.
More specifically for scDD used with pooled cells and the Kolmogorov-Smirnov test.
This method proved to be moderately robust to batch effects and still be among the most powerful methods
(among edgeR, MAST, D3E, BPSC, M3Drop).

## Installation

Directly form github using `devtools`:

    devtools::install_github('https://github.com/mRcSchwering/Lattirl')
