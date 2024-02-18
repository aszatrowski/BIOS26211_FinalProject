# BIOS 26211 Final Project: Dimension Reduction and Clustering of Gene Expression in Coronary Artery Disease
**Theresa Christiansen, Niranjan Joshi, Romy Li, Austin Szatrowski**

## Tasks
* [ ] Figure out what the numeric gene IDs are
* [ ] DR & Cluster each dataset individually - AS working on this
    * [ ] GSE20680:
        * [&check;] Data wrangling
        * [&check;] PCA
        * [&check;] $k$-means clustering
        * [&check;] Jaccard index
        * [ ] Cohen's $\kappa$
    * [ ] GSE20681:
        * [ ] Data wrangling
        * [ ] PCA
        * [ ] $k$-means clustering
        * [ ] Jaccard index & Cohen's $\kappa$
    * [ ] GSE12288:
        * [ ] Data wrangling (_can't get the `.gz` to uncompress LOL_)
        * [ ] PCA
        * [ ] $k$-means clustering
        * [ ] Jaccard index & Cohen's $\kappa$
    * Label points with CAD index (`GEO20680` and `GEO20681` have 3-point indices, so I think this is a good place to start)
    * Run Cohen's $\kappa$ and Jaccard to assess correlation between clusters and true labels
* [ ] Get batch effect removal for harmonization working (see paper, there's an R package for this)
    * Then harmonize datasets and run clustering on all
* [ ] Hyperparameter tuning and consensus choice of k?
* [ ] Comparison of PCA with tSNE and UMAP?

## File explainer:
* The files downloaded from GEO (accessions `12288`, `20680`, and `20681`) come in a json-like matrix format where each column vector is introduced by an `!` and then the cells are tab-delimited
* The `*_series_matrix.txt` data are the raw json-like files; `*_no-metadata.txt` is manually preprocessed to remove the metadata and make reading easier
* I was able to read these in with `read.table()` in R but all the categorical data (including CAD index) appears to be missing; it results in a matrix table with numerical gene IDs and samples. This can be clustered but we'll obviously need to find a way to add the labels in order to validate them.
    * This file is exported as `GSE20680.csv` and should be usable.
* The one `.gz` compressed file won't uncompress for some reason. No idea why.