# thermogenic_capacity_rnaseq
Code for analyzing RNAseq data for Velotta et al. This code processes gene expression data and creates a gene co-expression network from the data. It includes the raw counts table if gene expression and the co-expression network as a .RData file.

/data contains tables of read counts and metadata files
/scripts contains a single script for processing data and creating WGCNA co-expression network.
/results folder contains the .RData files containing the WGCNA network objects