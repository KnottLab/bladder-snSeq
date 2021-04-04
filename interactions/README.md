# Interactions

## Contents
- [Gene sets](#build-receptor-associated-gene-lists)
- [Receptor activity scores](#execute-receptor-activity-scoring)
- [Report active interactions](#run-interaction-analysis-to-identify-active-receptor-ligand-channels)
- [Demo](#minimial-demonstration)

----

Cell-cell interactions predicted from expressed ligand-receptor pairs in groups of heterogenous cells. 

To expand upon the concept of cell-cell interaction prediction from co-expression of ligand and receptor on two homogenous cell populations, we use the aggregate expression of genes implicated as downstream of receptor activation to read out differential receptor activity. 

## Build receptor-associated gene lists

First, we use the [`NicheNet`](https://www.nature.com/articles/s41592-019-0667-5)[1] dataset to connect receptors with interacting target proteins and ultimately genes associated at a gene-regulatory-network level. We make a deviation from the NicheNet paper, and use the Receptor-Ligand table provided by Cabello-Aguilar et al. in the [`SingleCellSignalR`](https://academic.oup.com/nar/article/48/10/e55/5810485)[2] package. 

In brief, we traverse two graphs in NicheNet, beginning from a given receptor, first gathering protein-protein signaling partners, then collecting gene-regulatory associated genes of those signalling proteins.

The particular gene lists used are provided in: `receptor_genelists_2020_06_30.tar.gz` where each receptor has a list of genes as a new line-delimited `txt` file.


```
usage: query_nichenet_receptors.py [-h] [--weighted_lr_sig WEIGHTED_LR_SIG]
                                   [--weighted_gr WEIGHTED_GR] [--steps STEPS]
                                   [--n_sig_partners N_SIG_PARTNERS]
                                   receptors_fname output_dir

positional arguments:
  receptors_fname       newline-delimited list of receptors to search
  output_dir            someplace to stash the results. it will be created in
                        case it does not exist.

optional arguments:
  -h, --help            show this help message and exit
  --weighted_lr_sig WEIGHTED_LR_SIG
                        a three-column table with from, to, and weight columns
                        representing directed connections between proteins
  --weighted_gr WEIGHTED_GR
                        a three-column table with from, to, and weight columns
                        representing gene regulatory relationships
  --steps STEPS         the number of moves away from each start node to
                        consider. default=2
  --n_sig_partners N_SIG_PARTNERS
                        the maximum number of elements to take at each step.
                        default=10
```


```
[1] Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods 17, 159–162 (2020). https://doi.org/10.1038/s41592-019-0667-5

[2] Simon Cabello-Aguilar, Mélissa Alame, Fabien Kon-Sun-Tack, Caroline Fau, Matthieu Lacroix, Jacques Colinge, SingleCellSignalR: inference of intercellular networks from single-cell transcriptomics, Nucleic Acids Research, Volume 48, Issue 10, 04 June 2020, Page e55, https://doi.org/10.1093/nar/gkaa183
```

-------

## Execute receptor activity scoring


The script `score_genelists_run.py` implements receptor activity scoring. Receptor scoring assumes a readout of differential gene expression as a result of signaling downstream of receptor activation (or removing receptor activation). In contrast to a typical single cell gene set scoring, we do not assume a shared directionality in differential gene expression, and take the absolute difference to reflect this. This way, gene sets that are half activated, and half repressed do not cancel out, and we instead read out total deviation from the average expression in a background set of cells. To appropriately assess the level of differential expression in single cells, we use the average GEX from cells of the same general type as the cell in question.  


```
usage: score_genelists_run.py [-h] [--out OUT] [-j N_JOBS]
                              genelist_dir adata_path groupby

positional arguments:
  genelist_dir          A directory containing gene lists in newline delimited
                        *.txt format
  adata_path            AnnData object holding gene expression
  groupby               A column in adata.obs to use as cell groups for
                        background gene expression

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             an h5ad file to stash results. if not specified, it
                        defaults to converting the genelist_dir path into a
                        file name
  -j N_JOBS, --n_jobs N_JOBS
                        Number of parallel jobs to launch. default=8
```

-----

## Run interaction analysis to identify active receptor-ligand channels

First use `01_Filter_Ligand_Receptor_coincidence.ipynb` to generate the list `visium_channels_2.txt` that lists receptor/ligand pairs that are colocated. 

The notebook `02_Circos_interactions.ipynb` uses the subtype annotation to group the ligand expression in a set of "sender" cells, and the receptor scores in a group of "receiver" cells. To keep interpretation of these analysis managable, we restrict each plot to a single sending subtype, and one "broad" receiving population. 

The procedure is as follows:

1. Given receptor activity scores for each cell, identify with a differential expression-style analysis the active receptors in each subtype of cells.
2. Similarly, use differential expression to identify sets of ligands that are enriched in each subtype of cells.
3. Using the focused "sender" and "receiver" populations, identify the overlapping active receptors and enriched ligands.
4. Write out the results.

We call [`circos`](http://circos.ca/)[3] to visualize interacting channels.


```
[3] Circos: An information aesthetic for comparative genomics
Martin I Krzywinski, Jacqueline E Schein, Inanc Birol, Joseph Connors, Randy Gascoyne, Doug Horsman, Steven J Jones, and Marco A Marra
Genome Res. Published in Advance June 18, 2009, doi:10.1101/gr.092759.109
```


## Minimial demonstration

For demonstration purposes, a heavily subsetted version of the MIBC snSeq dataset was constructed using foreknowledge of celltypes, active receptors, and expressed ligands in each cell type. This dataset will allow us to create a facsimile of the CDH12-Fibroblast circos plot in Figure 3E.

Before running these scripts: 
- install [`ray`](https://github.com/ray-project/ray/). Ray is used to parallelize the receptor scoring procedure.
- Obtain the ligand receptor table from SingleCellSignalR.

`demo/create_demo_data.ipynb` creates files needed to run the demo. These are derived from the full MIBC snSeq dataset. We subset genes by taking genes related to the receptor-ligand channels presented in the main figures (`demo/use_receptors.txt`). Nuclei are subsetted uniformly taking at most 500 cells from each subtype.

We arrive at the following sample data:

```
AnnData object with n_obs × n_vars = 7328 × 1444
```

`demo/demo_receptor_scoring.ipynb` demonstrates running the receptor scoring script on the sample data created above.
We end up with another AnnData object containing receptor scores. Note also that all `obs` annotations are carried over to the receptor score AnnData:

```
AnnData object with n_obs × n_vars = 7328 × 23
```

Finally, the filtering and visualization notebook can be run. The notebook `demo_Circos_interactions.ipynb` will create data files for circos. After it runs, if circos is installed, run the `build.sh` script in the output directory.


