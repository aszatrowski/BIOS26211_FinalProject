# BIOS 26211 Final Project: Dimension Reduction and Clustering of Gene Expression in Coronary Artery Disease
**Theresa Christiansen, Niranjan Joshi, Romy Li, Austin Szatrowski**

## Tasks
* Figure out what the numeric gene IDs are
* DR & Cluster each dataset individually
    * AS starting this
    * Label points with CAD index (`GEO20680` and `GEO20681` have 3-point indices, so I think this is a good place to start)
    * Run Cohen's $\kappa$ and Jaccard to assess correlation between clusters and true labels
* Get batch effect removal for harmonization working (see paper, there's an R package for this)
    * Then harmonize datasets and run clustering on all

## File explainer:
* The files downloaded from GEO (accessions `12288`, `20680`, and `20681`) come in a json-like matrix format where each column vector is introduced by an `!` and then the cells are tab-delimited
    * I was able to read these in with `read.table()` in R but all the categorical data (including CAD index) appears to be missing. 
    