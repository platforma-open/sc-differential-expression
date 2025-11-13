# Overview

Identifies differentially expressed genes between two groups of single cells (e.g., different genotypes, treatments, or conditions) using the Wilcoxon rank-sum test. The block performs statistical testing to compare gene expression levels across groups, enabling the discovery of genes that respond to experimental conditions or distinguish cell populations.

The block supports multiple comparisons by testing each numerator group against a common denominator baseline. Genes are filtered based on configurable thresholds: minimum log2 fold change and maximum adjusted p-value. Results include comprehensive differential expression tables with log2 fold changes, p-values, and adjusted p-values, as well as filtered tables containing only statistically significant genes suitable for downstream functional enrichment analysis and visualization.

The block uses scanpy v1.10.1 for differential expression analysis. When using this block in your research, cite the scanpy publication (Wolf et al. 2018) listed below.

The following publication describes the methodology used:

> Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. _Genome Biology_ **19**, 15 (2018). [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
