BarcodeGeneration_HTOfilteringpart1.R - R script that  generates barcode table from Cellranger output for input into CITEseqCount. Then takes output from CITEseqCount and filters out doublets.

HTOfilteringpart2.m - Matlab script that takes output from BarcodeGeneration_HTOfilteringpart1.R and recovers nuclei that were labelled as "negative" by Seurat script if the nuclei satisfy specific parameters set forth in the script.

Figure5A_Figure5B.ipynb - python Jupyter notebook that outlines the steps for creating the Visium-derived gene expression modules (Figure 5A) as well as the spot-radius gene expression plots (Figure 5B)