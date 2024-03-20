# Scripts used for analysis

## R_analysis.r 
This script contains basic code for:

* Differential expression analysis using DESeq2 R package 
* Splicing analysis using rMATS2 R package

Required libraries (installed when script is run):
* limma
* edgeR
* statmod
* DESeq2
* gplots
* biomaRt
* RColorBrewer
* tidyverse
* pheatmap

## calculate_overlaps_and_plot.R
This script contains the code used to analyse the overlap between differential expression results and genes from other databases, e.g. the SFARI database. Assuming that you have saved the DE results from our study as "DE_results.csv", and the SFARI genes as "SFARI_genes.csv", the script can be run like so:

```
Rscript calculate_overlaps_and_plot.R DE_results.csv Symbol gene-symbol SFARI_genes.csv
```

In order, the commandline arguments are:
1. DE results file in CSV/TSV format
2. The name of the column specifying the results genes being used in the overlap analysis. Must be in the same format as those in the database of interest. In this case, we use HGNC symbols. Note: only acceptable values are "Ensembl" or "Symbol" (case insensitive)
3. The name of the column specifying the database genes being used in the overlap analysis. Must be in the same format as those in the DE results (HGNC symbol or ENSEMBL ID).
4. (and 5, 6, ..., N) files with genes of interest from databases

Plots and tables will be generated.

Required libraries:
* data.table
* dplyr
* clusterProfiler
* ggplot2
* org.Hs.eg.db
* knitr

If you wish to cut down on dependencies, some of the code can be refactored (e.g., for table i/o and printing).

## run_disgenet.R
This script contains the code used to run disease enrichment analysis using the list of differentially expressed genes. The script can be run like so:

```
Rscript run_disgenet.R
```

Note: the DE results are hardcoded in as `Merlin.res.df.txt`. If you have saved the results as anything else, make sure to update the path within the script.

Required libraries:

* data.table
* disgenet2r
* dplyr
* ggplot2

## Authorship note
Code here was written by Dr Eunbi Lee
