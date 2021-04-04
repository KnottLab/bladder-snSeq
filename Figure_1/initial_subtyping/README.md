# QC and preprocessing

Data should be processed with cellranger version 3.

These notebooks and scripts document initial QC, filtering and clustering/subtyping of nuclei.

-----

`BarcodeGeneration_HTOfilteringpart1.R` - R script that  generates barcode table from Cellranger output for input into CITEseqCount. Then takes output from CITEseqCount and filters out doublets.

`HTOfilteringpart2.m` - Matlab script that takes output from `BarcodeGeneration_HTOfilteringpart1.R` and recovers nuclei that were labelled as "negative" by Seurat script if the nuclei satisfy specific parameters set forth in the script.


The notebook `filter_PercentMito_nUMI.ipynb` shows the filtering from demux'd cells by percent reads mapping to mitochondrial genes and bandpass filtering for number of UMIs per cell.

# Clustering and subtyping

Dimensionality reduction and clustering were performed with single cell Variational Inference and the leiden algorithm.

1. Starting from a dataset in `AnnData` format called `all_cells.h5ad`:
```bash
(scrna) $ ./run_scvi_AllCells.sh all_cells.h5ad
# creates all_cells.scvi/ and all_cells.clust.h5ad
```

2. Following initial subtyping of the clusters by marker genes each major cell type was separately sub-clustered, again with scVI and the leiden algorithm.

3. Each set of subclusters were assigned a subtype name based on differentially expressed and known marker genes

4. We gather the annotations, dropping unsubtypbale and likely doublet cell clusters, and arrive at the subtype annotations
