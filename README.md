# bladder_snSeq


Code release: [![DOI](https://zenodo.org/badge/323445511.svg)](https://zenodo.org/badge/latestdoi/323445511)

Nature Communications paper: https://www.nature.com/articles/s41467-021-25103-7

[bibtex citation](#citation)


## Data
Data are available through GEO, the paper's supplemental materials, or by request:

single nucleus RNAseq: `GSE169379`

Visium: `GSE171351`

CODEX processed data:
- TMA1: https://figshare.com/s/4610a15363c8306dfa36
- TMA2: https://figshare.com/s/2005255a8b65de23109f
- TMA3: https://figshare.com/s/1d8c7ed76d4b3222ada4

Please contact the corresponding authors for original CODEX images


## System requirements

See `pip_freeze.txt` for software used and particular versions.

`scVI` runs optionally with GPU acceleration turned on, and receptor scoring scripts assume an environment with multiple CPU threads.

We make use of GPU accelerated implementations of UMAP embedding, nearest-neighbors finding, and leiden clustering provided by RAPIDS.ai throughout.

Analysis software was developed and run on Ubuntu 16.04 and 20.10. Untested but expected to run without modification on OSX/MacOS. Requires some modification to path strings in order to run on Windows.

CODEX image preprocessing requires MATLAB and the Image Processing Toolbox, and was run on Windows.


## Installation

Install `jupyter`, instantiate a `jupyter` server in the same directory as this README, view and run individual notebooks.

Whenever possible, all custom scripts are included in this repository and should be accessible in the applicable notebooks/scripts. [Conda](https://docs.conda.io/en/latest/miniconda.html) was used to manage the python environments. Two conda environments were used: one for the single nucleus RNA-seq and visium analysis, and another for the CODEX analysis. The details of these environments are detailed in `pip_freeze.txt` for the snSeq/visium, and in `codex_env.yml` for the CODEX analysis. A single environment, built based on `codex_env.yaml` is expected to be applicable usable for the entire repository, however this has not been tested.

Custom MATLAB code for CODEX image preprocessing are available through a separate repository: https://github.com/KnottLab/codex

An accompanying toolbox was developed to facilitate handling, visualization and analysis of preprocessed CODEX images, and is available through a separate repository: https://github.com/nathanin/micron2

To run the scripts in this repository that depend on `micron2`, clone the repository then install it in editable mode via pip:
``` bash
(env) foo@bar $ git clone https://github.com/nathanin/micron2
(env) foo@bar $ cd micron2
(env) foo@bar $ pip install -e .
```

## Demo

Each notebook shows an analysis beginning with source data, and ending, usually, with a figure panel. 


## Use

Scripts for receptor scoring are ready to run on any properly formatted `AnnData` dataset. Interaction analyses may be run after receptor scoring has been performed with some modification to the indicated notebook.


------

## Analysis

Listing of key analysis/plotting software used and notebooks/scripts containing significant custom code. Subdirectories also have separate README files with further details.


| Figure/Analysis | Description | Required software | Script/Notebook/Data |
|----------------:|:------------|:---------------------|:---------------------|
|Preprocessing| Hashing demultiplex|-|`initial_subtyping/HTOfilteringpart2.m`, `BarcodeGeneration_HTOfiltertingpart1.R`|
|Fig 1A| workflow | n/a | n/a |
|Fig 1B|UMAP|`scvi`, `scanpy`|`initial_subtyping/run_scvi_AllCells.sh`|
|Fig 1C|heatmap|`seaborn`|`initial_subtyping/all_cells.ipynb`|
|Fig 1D|epithelial compartment umap|`scvi`, `scanpy`|`initial_subtyping/epithelial_compartment.ipynb`|
|Fig 1E|gene signatures|`scanpy`|gene sets in `Figure_1/` directory|
|Fig 1F|epithelial differentiation marker dot plot|`scanpy`|-|
|Fig 1G|epithelial co-expression modules|`seaborn`|-|
|Fig 1H|SCENIC regulon scores|`pySCENIC`, `seaborn`|`Figure_1/scenic/`|
|Fig 1I|pluripotency gene scores|`scanpy`|gene sets in `Figure_1/` directory|
|---|---|---|---|
|Fig 2A|UMAP|`scvi`, `scanpy`|`Figure_2/Figure_2AB.ipynb`|
|Fig 2B|marker genes|`scanpy`|`Figure_2/Figure_2AB.ipynb`|
|Analysis|RNA velocity analysis|`velocyto`, `scVelo`|`Figure_2/Figure_2_scVelo.ipynb`|
|Fig 2C|sample A RNA velocity|`velocyto`, `scVelo`, `scanpy`|`Figure_2/Figure_2C.ipynb`|
|Fig 2D|sample A latent time|`scanpy`|`Figure_2/Figure_2DE.ipynb`|
|Fig 2E|sample A subtype density & heatmap|`scanpy`, `seaborn`|`Figure_2/Figure_2E.ipynb`|
|Analysis| normal/tumor nearest neighbor search | - |`Figure_2/Figure_2_NearestNeighbor_search.ipynb`|
|Fig 2F|subtype density & time-signature based KM-plots|`scanpy`, `seaborn`, `lifelines`|`Figure_2/Figure_2F.ipynb`|
|---|---|---|---|
|Analysis|derivation of snSeq subtype gene signatures|-|`Figure_3/gene_signatures.ipynb`|
|Fig 3A|snSeq gene signature scores in TCGA|-|`Figure_3/Figure_3A.ipynb`|
|Fig 3B|TCGA KM-plots|`lifelines`|`Figure_3/Figure_3B.ipynb`|
|Fig 3C|signature enrichment change pre-to-post NAC|-|-|
|Fig 3D|gene sets enriched post-NAC|`limma`, `shinyGO`|-|
|Analysis|Receptor scoring|-|`interactions/`|
|Analysis|Spatially coincidental receptors/ligands|-|`interactions/`|
|Fig 3E|snSeq CDH12-Fibroblast subtype interactions|`circos`|`interactions/`|
|---|---|---|---|
|Fig 4A|RMA intensity boxplots|`seaborn`|-|
|Fig 4B|PDL1, PDL2 snSeq dot plots|`scanpy`|-|
|Fig 4C|IMvigor210 pre/post chemo KM-plots|`lifelines`|same as `Figure_3B`|
|Fig 4D|IMvigor210 CDH12 association with response|`scipy`|-|
|Fig 4E|comparison of molecular signatures|`scipy`|-|
|Fig 4F|snSeq CDH12-Tcell subtype interactions|`circos`|`interactions/`|
|Fig 4G|Interaction potential of checkpoint molecules|`scanpy`|`Figure_4/Figure_4G.ipynb`|
|---|---|---|---|
|Fig 5A|visium spatial gene expression|`scanpy`, `spatialIDE`|`Figure_5/Figure_5A_Figure5B.ipynb`|
|Fig 5B|CODEX data schematic|-|-|
|Fig 5C|Inter-cellular distances|`micron2`|`Figure_5_codex/spatial_stats.ipynb`|
|Fig 5D|CODEX Voronoi diagrams|`micron2`, `matplotlib`|`Figure_5_codex/voronoi.ipynb`|
|Fig 5E|Cellular Niche subtype diversity|-|`Figure_5_codex/diversity.ipynb`|
|Fig 5F|CD8T marker intensity enrichment|`micron2`|`Figure_5_codex/marker_expression.ipynb`|
|Fig 5G|sample images|-|-|
|Fig 5H|CDH12 cell marker intensity enrichment|`micron2`|`Figure_5_codex/marker_expression.ipynb`|
|---|---|---|---|
|Fig 6A|comparison of molecular signatures|`scipy`|-|
|Fig 6A|clinical flowchart incorporating CDH12|n/a|-|
|---|---|---|---|
|Extended Data Fig 1A|MIBC snSeq QC|-|`Ext_Figure_1/QCplots.ipynb`|
|Extended Data Fig 1B|normal bladder scrublet scores|`scrublet`, custom scripts|`Ext_Figure_1/scrublet/`|
|Extended Data Fig 1C|patient subtypes|-|`Ext_Figure_1/patient_pie_charts.ipynb`|
|Extended Data Fig 1E|epithelial compartment heatmap|`seaborn`|`initial_subtyping/epithelial_compartment.ipynb`|
|Extended Data Fig 1G|fibroblast compartment UMAP|`scvi`, `scanpy`|`initial_subtyping/fibroblast_compartment.ipynb`|
|Extended Data Fig 1H|immune compartment UMAP|`scvi`, `scanpy`|`initial_subtyping/immune_compartment.ipynb`|
|Extended Data Fig 1I|fibroblast compartment heatmap|`seaborn`|`initial_subtyping/fibroblast_compartment.ipynb`|
|Extended Data Fig 1J|immune compartment heatmap|`seaborn`|`initial_subtyping/immune_compartment.ipynb`|
|---|---|---|---|
|Extended Data Fig 2A|KRT13 / KRT17 IHC|-|-|
|---|---|---|---|
|Extended Data Fig 3A|CDH12 / CDH18 IHC|-|-|
|---|---|---|---|
|Extended Data Fig 4B|ALDH1 expression in epithelial compartment cells|`scanpy`|-|
|Extended Data Fig 4C|non-tumor snSeq UMAPs|`scvi`, `scanpy`|`initial_subtyping/normal_snSeq.ipynb`|
|Extended Data Fig 4D|MIBC & non-tumor dot plot|`scanpy`|-|
|Extended Data Fig 4E|non-tumor cell type density|`velocyto`, `scVelo`, `seaborn`|`Figure_2/Figure_2E.ipynb`|
|---|---|---|---|
|Extended Data Fig 5A|GSVA scores|`gsva`, `seaborn`|-|
|Extended Data Fig 5C|GSVA score KM-plots|`gsva`, `gseapy`, `lifelines`|-|
|---|---|---|---|
|Extended Data Fig 6B|IMvigor210 ssGSEA KM-plots|`gseapy`, `lifelines`|-|
|Extended Data Fig 6C|visium QC plots|`scanpy`|-|
|Extended Data Fig 6D|visium derived gene signatures in snSeq and in-situ|`scanpy`|`Figure_5/Figure_5A_Figure5B.ipynb`|
|---|---|---|---|
|Extended Data Fig 7A|CODEX type gating|-|`Figure_5_codex/gate_celltypes.ipynb`|
|Extended Data Fig 7B|Classified subtype marker intensity|-|`Figure_5_codex/marker_expression.ipynb`|
|Extended Data Fig 7C|Identifying the number of niches to use|-|`Figure_5_codex/niche_calling.ipynb`|
|Extended Data Fig 7D|Subtype enrichment within Cellular Niches|-|`Figure_5_codex/niche_associations.ipynb`|
|---|---|---|---|
|Extended Data Fig 8,9,10|Dot plots of celltype, subtype and cellular niches on TMA spots|-|`Figure_5_codex/make_layout.ipynb`|

-----

## Citation

```
@article{gouin2021n,
  title={An N-Cadherin 2 expressing epithelial cell subpopulation predicts response to surgery, chemotherapy and immunotherapy in bladder cancer},
  author={Gouin, Kenneth H and Plummer, Jasmine T and Rosser, Charles J and Ben Cheikh, Bassem and Oh, Catherine and Chen, Stephanie S and Chan, Keith Syson and Furuya, Hideki and Tourtellotte, Warren G and Knott, Simon RV and others},
  journal={Nature Communications},
  volume={12},
  number={1},
  pages={1--14},
  year={2021},
  publisher={Nature Publishing Group}
}
```
