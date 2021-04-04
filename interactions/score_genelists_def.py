import numpy as np
import pandas as pd
import scanpy as sc

import os
from glob import glob
import argparse

import ray
import logging

# adata_path = '../samples_trimmed/pembroRT_TNBC_AllCells.clust.h5ad'
# groupkey = 'CellType'
# keep = ['Tcell']

# # Key to group cells for calculating a background
# groupby_key = 'CellType_v3'

# adata = sc.read_h5ad(adata_path)
# print(f'Loaded adata {adata.shape}')

# # if groupkey in adata.obs.columns:
# #   subsample = adata.obs[groupkey].isin(keep)
# #   adata = adata[subsample, :]
# #   print(f'Subsampled adata {adata.shape}')

# ## Subset genes to clean up the background calculation
# sc.pp.filter_genes(adata, min_cells = 100)
# print(f'Reduced by thresholding the number of genes: {adata.shape}')

# ## Normalize counts per cell
# sc.pp.normalize_total(adata, target_sum=10000)
# sc.pp.log1p(adata)


# # Get a list of background genes and rank them for the rank-guided background
# # selection later  
# # REF: https://github.com/theislab/scanpy/blob/master/scanpy/tools/_score_genes.py
# all_genes = set(adata.var_names)

# ## Binning method i.e. Seurat and Scanpy
# gene_avg = pd.Series(np.squeeze(np.array(np.mean(adata.X, axis=0))), index=all_genes)
# n_items = len(gene_avg) / 20 ## <-- expression bins; use fewer bins to go faster
# gene_cut = gene_avg.rank(method='min') // n_items



@ray.remote
def score_genelist(gex, var_names, cell_groups, gene_list):
  # logger = logging.getLogger()
  logging.basicConfig(level='INFO')

  list_name = os.path.splitext(os.path.basename(gene_list))[0]

  if list_name not in var_names:
    logging.info(f'WARN: gene {list_name} not found in adata.var_names. skipping')
    return None

  with open(gene_list, 'r') as f:
    genes = [l.strip() for l in f]
  og_len = len(genes)

  # Restrict to genes in the dataset
  genes = [g for g in genes if g in var_names]
  red_len = len(genes)
  logging.info(f'{list_name:<10} loaded {og_len} restricted to {red_len}')

  if red_len < 5: 
    logging.info(f'WARN got {len(genes)} genes for list {list_name}. Skipping.')
    return None

  genes = set(genes)


  X_genes = gex[:, var_names.isin(genes)].toarray()
  X_base_gene = np.squeeze(gex[:, var_names == list_name].toarray()) # Expression of the receptor, per cell.

  ## ------------------    Difference from cell type average expression method
  ## I for sure know theres a less clunky way to do this 
  ## that probably has to do with diagonals
  ## 1. Find the average expression of the genes in our set, per cell type
  # groups = adata.obs[groupby_key].values
  u_groups = np.unique(cell_groups)
  group_avgs = {}
  for u in u_groups:
    avg = np.mean(gex[cell_groups == u][:, var_names.isin(genes)].toarray(), axis=0)
    group_avgs[u] = avg # 1 x n_genes

  # Set up the difference matrix
  gene_diffs = np.zeros_like(X_genes)
  # print(f'INFO gene diffs set up ~ {gene_diffs.shape} {gene_diffs.dtype}'
  for i, g in enumerate(genes):
    # pull out the gene's expression vector for all cells
    X_gene = X_genes[:, i]

    # This variable will hold the deviation from average for each gene
    X_gene_difference = np.zeros_like(X_gene)
    
    # Loop over groups 
    for u in u_groups:
      indexer = cell_groups == u
      avg = group_avgs[u][i]
      X_gene_group = X_gene[indexer]
      diff = np.abs(avg - X_gene_group)             # <------------- Absolute value difference from mean
      # diff = np.abs(avg - X_gene_group) / avg     # <------------- Percent difference from mean
      diff[np.isnan(diff)] = 0
      X_gene_difference[indexer] = diff
    
    gene_diffs[:, i] = X_gene_difference

  # We want to make sure the output is 1D
  scores = np.squeeze(np.mean(gene_diffs, axis=1))
  scores = scores * (X_base_gene > 0).astype(np.float32) # <--- Add +1 to include all cells. make binary to not scale by expression.
  
  # print(f'INFO Returning list {list_name} {scores.shape} {scores.dtype}')


  ## ----------------------- Difference from bin-matched expression i.e. Seurat / Scanpy 
  # background_genes = set()
  # gene_ranks = gene_cut.loc[genes].values
  # for rnk in np.unique(gene_ranks):
    # rank_genes = np.array(gene_cut[gene_cut == rnk].index)
    # n_choice = min(len(rank_genes), 100)  ## <-- number of background genes per rank ; use fewer genes to go faster
    # background_genes.update(set(np.random.choice(rank_genes, n_choice, replace=False)))
  # background_genes = background_genes - genes
  # X_background = adata.X[:, adata.var_names.isin(background_genes)]
  # scores = np.squeeze(np.mean(X_genes, axis=-1) - np.mean(X_background, axis=-1))


  ## ---------------------- Straight up mean expression
  # scores = np.squeeze(np.array(np.mean(X_genes, axis=-1)))
  # scores = scores * X_base_gene

  return {list_name: scores}
