#!/usr/bin/env python

import numpy as np
import scanpy as sc
import pandas as pd

import argparse
import os

from load_data import load_data

"""
Modules have two modes: standalone from command line and pipelined

Both mode accepts a preprocessed AnnData object as input.
  - Standalone mode writes back a AnnData with new metadata
  - Pipelined mode returns the AnnData object with new metadata

Preprocessing
The preprocessing pipeline expects a raw, or raw-like dataset.
We read it, clean it, and return it.

Standard cleaning: 
  - Caclulate basic QC stats
  - Remove cells with high mitochondrial gene %

  (optional)
  - Doublet detection, removal (e.g. Scrublet)
  - Doublet removal by expected rate + nUMI count

  - Count normalization 
  - Log1p
  - Evaluate cell cycle [1]
  - set raw attribute of AnnData

  - Batch correction
    -> combat (scanpy)
    -> mnn (pymnn)
    -> scanorama  

  - Find highly variable genes, subset

REFS:
  [1] https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
  scanpy
  Seurat

"""

def add_preprocess_args(parser):
  """ Add preprocessing args so they're accessible from other scripts  """
  parser.add_argument('--mito_pct',  default=0.1,    type=float)

  parser.add_argument('--min_cells', default=20,    type=int)
  parser.add_argument('--min_genes', default=20,     type=int)
  parser.add_argument('--min_mean',  default=1/32, type=float)
  parser.add_argument('--max_mean',  default=5,      type=int)
  parser.add_argument('--min_disp',  default=0.25,   type=float)

  parser.add_argument('--mouse', action='store_true',
                      help = 'Whether to treat the sample as mouse.')

  parser.add_argument('--filter_mito', default=1, type=int,
                      help = 'Whether to filter cells by percent mitochondrial genes.')

  parser.add_argument('--estimate_cell_cycle', default=1, type=int,
                      help = 'Whether to estimate cell cycle from gene expression.')

  parser.add_argument('--count_norm', default=1, type=int,
                      help = 'Whether to count-normalize the reads per cell.')

  parser.add_argument('--log1p', default=1, type=int,
                      help = 'Whether to log-normalize counts.')

  parser.add_argument('--hvg', default=1, type=int,
                      help = 'Whether to subset highly variable genes.')

  parser.add_argument('--blacklist_genes', default=None, type=str,
                      help = 'Provided a list of genes to blacklist from analysis.')

  parser.add_argument('--require_genes', default=None, type=str,
                      help = 'Provided a list of genes to include after HVG. Use this with force_n_genes=0 to pass a pre-determined list')

  parser.add_argument('--force_n_genes', type=int,   default=None,
                      help = 'Force use of N highly variable genes for dim. reduction.')

  parser.add_argument('--batch_correction', type=str, default=None,
                      help = 'Type of batch correction procedure to run before HVG selection..')

  parser.add_argument('--batch_key', type=str, default='batch',
                      help = 'Key in `obs` to partition cells during HVG selection. \
                      https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.pp.highly_variable_genes.html')

  return parser

CELL_CYCLE_GENES = [
'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1',
'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7',
'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 
'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8', 'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2',
'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 
'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1',
'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR',
'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']


# Translated from human to mouse by KHG 
MOUSE_G2M_GENES = [
  'Bub1', 'Cenpf', 'Aurkb', 'Cdca8', 'Anln', 'Smc4', 'Nusap1', 'Cbx5', 'Rangap1', 'Ckap5', 'Hmmr',
  'Anp32e', 'Hmgb2', 'Gtse1', 'Ckap2l', 'Ttk', 'G2e3', 'Cks1b', 'Dlgap5', 'Ccnb2', 'Tubb4b', 'Cdca2',
  'Ctcf', 'Hjurp', 'Tpx2', 'Kif23', 'Cdc20', 'Top2a', 'Ncapd2', 'Nuf2', 'Mki67', 'Cks2', 'Ckap2', 'Tmpo',
  'Aurka', 'Ndc80', 'Lbr', 'Tacc3', 'Ube2c', 'Cenpe', 'Cenpa', 'Psrc1', 'Kif20b', 'Birc5', 'Cdca3', 'Gas2l3',
  'Cdk1', 'Kif2c', 'Ect2', 'Nek2', 'Kif11', 'Cdc25c'
]

MOUSE_S_GENES = [
  'Hells', 'Blm', 'Rpa2', 'Cdc45', 'Pola1', 'Gmnn', 'Rrm2', 'Wdr76', 'Tyms', 'Tipin', 'Mcm5', 'Nasp',
  'Dscc1', 'Ccne2', 'Mcm2', 'Ung', 'Cdca7', 'Chaf1b', 'Dtl', 'Gins2', 'Mcm6', 'Rad51ap1', 'Casp8ap2', 'Mcm4',
  'Fen1', 'Clspn', 'Usp1', 'Slbp', 'Msh2', 'Prim1', 'Pcna', 'Cdc6', 'Uhrf1', 'Rrm1', 'Ubr7', 'Brip1', 'Exo1',
  'E2f8', 'Rfc2', 'Pcna-ps2'
]

# def run_preprocess(adata, 
#     batch_correction=None, 
#     count_norm=True, 
#     log1p=True, 
#     hvg=True, 
#     zero_center=False, 
#     force_n_genes=None,
#     ARGS=None):
def run_preprocess(adata, ARGS=None):
  # adata.var_names_make_unique()

  # Use an unstructured annotation to track preprocessing ops
  adata.uns['preprocessing_ops'] = []

  # this_file = os.path.abspath(__file__)
  # this_dir = os.path.dirname(this_file)

  # cell_cycle_gnes = [x.strip() for x in open(f'{this_dir}/regev_lab_cell_cycle_genes.txt')]
  s_genes = CELL_CYCLE_GENES[:43]
  g2m_genes = CELL_CYCLE_GENES[43:]
  
  sc.pp.calculate_qc_metrics(adata, inplace=True)
  adata.uns['preprocessing_ops'].append('QC')

  if ARGS.filter_mito:
    mt_pattern = "mt-" if ARGS.mouse else "MT-"
    mito_genes = adata.var_names.str.startswith(mt_pattern)
    adata.obs["percent_mito"] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1
    passed_mito_filter = adata.obs["percent_mito"] < 0.1
    print(f'Mito genes filter passed by {np.sum(passed_mito_filter)} cells.')
    adata = adata[passed_mito_filter, :] ## Mitochondrial genes
    adata.uns['preprocessing_ops'].append('mito_filter')

  sc.pp.filter_cells(adata, min_genes=ARGS.min_genes)
  adata.uns['preprocessing_ops'].append('filter_cells')
  print(f'[RUN_PREPROCESSING] Filtered cells by min_genes ({ARGS.min_genes}): {adata.shape}')

  adata.raw = adata
  adata.uns['preprocessing_ops'].append('raw_set')

  # NOTE ----- do this before or after count normalization  ?
  sc.pp.filter_genes(adata, min_cells=ARGS.min_cells)
  print(f'[RUN_PREPROCESSING] Filtered genes by min_cells ({ARGS.min_cells}): {adata.shape}')

  if ARGS.blacklist_genes is not None:
    print(f'[RUN_PREPROCESSING] blacklisting genes from {ARGS.blacklist_genes}')
    blacklist = [l.strip() for l in open(ARGS.blacklist_genes, 'r')]
    adata = adata[:, ~adata.var_names.isin(blacklist)]
    print(f'[RUN_PREPROCESSING] blacklisted {len(blacklist)} genes (adata: {adata.shape})')

  if ARGS.count_norm:
    print('[RUN_PREPROCESSING] count norm')
    sc.pp.normalize_total(adata, target_sum=10000, exclude_highly_expressed=False)
    adata.uns['preprocessing_ops'].append('count_norm')
  else:
    print('[RUN_PREPROCESSING] SKIPPING count norm')

  if ARGS.log1p:
    print('[RUN_PREPROCESSING] log1p')
    sc.pp.log1p(adata)
    adata.uns['preprocessing_ops'].append('log1p')
  else:
    print('[RUN_PREPROCESSING] SKIPPING log1p')


  if ARGS.estimate_cell_cycle:
    s_genes = [x for x in s_genes if x in adata.var_names]
    g2m_genes = [x for x in g2m_genes if x in adata.var_names]
    if len(s_genes) == 0 or len(g2m_genes) == 0:
      print('WARNING: no S-phase or G2M genes found in var_names. Is this a mouse?')
      s_genes = [x for x in MOUSE_S_GENES if x in adata.var_names]
      g2m_genes = [x for x in MOUSE_S_GENES if x in adata.var_names]

    if len(s_genes) > 0 and len(g2m_genes) > 0:
      print('[RUN_PREPROCESSING] processing S and G2M scores.')
      sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
      adata.uns['preprocessing_ops'].append('cell_cycle_score')
    else:
      print('[RUN_PREPROCESSING] WARNING: No cell cycle genes in var_names.')
  else:
    print('[RUN_PREPROCESSING] SKIPPING cell cycle estimation')


  if ARGS.batch_correction is not None and ARGS.batch_correction == 'combat':
    sc.pp.combat(adata)
    adata.uns['batch_correction_method'] = 'combat'
    adata.uns['preprocessing_ops'].append('combat')


  if ARGS.hvg:
    #if ARGS.batch_key in list(adata.obs.keys()):
    #  batch_key = ARGS.batch_key
    #  # Force batch_key to be categorical
    #  batch_values = pd.Series(adata.obs[batch_key].values, dtype='category')
    #  adata.obs[f'{batch_key}_categorical'] = batch_values
    #  batch_key = f'{batch_key}_categorical' 
    #else:
    #  print(f'[RUN_PREPROCESSING] WARNING batch key {ARGS.batch_key} not in keys. Running HVG on all batches')
    #  batch_key = None
    batch_key=None

    if ARGS.force_n_genes == 0 and ARGS.require_genes is not None:
      passing_genes = set([])
    else:
      if ARGS is not None and ARGS.force_n_genes is not None:
        sc.pp.highly_variable_genes(adata, n_top_genes=ARGS.force_n_genes, subset=False, batch_key=batch_key)
        adata.uns['preprocessing_ops'].append('n_top_genes')
      else:
        sc.pp.highly_variable_genes(adata, min_mean=ARGS.min_mean, max_mean=8, 
          min_disp=ARGS.min_disp, 
          subset=False, batch_key=batch_key)
        adata.uns['preprocessing_ops'].append('hvg_by_gene_traits')

      passing_genes = set(adata.var_names.values[adata.var["highly_variable"]])

    if ARGS.require_genes is not None:
      assert os.path.exists(ARGS.require_genes)
      required_genes = set([l.strip() for l in open(ARGS.require_genes, 'r')])
      adata.uns['required_genes'] = list(required_genes)
    else:
      required_genes = set([])

    passing_genes = passing_genes.union(required_genes)
    assert len(passing_genes) > 0

    adata = adata[:, adata.var_names.isin(passing_genes)]
    print(f'[RUN_PREPROCESSING] selected hvg: {adata.shape}')

    sc.pp.filter_genes(adata, min_cells=50) 
    print(f'[RUN_PREPROCESSING] filtering genes by min_cells ( > 50 ): {adata.shape}')
    adata.uns['preprocessing_ops'].append('HVG')

  # Set a flag in adata to tell us we're writing a processed file
  adata.uns['PREPROCESSED_FLAG'] = True
  return adata


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('dataset', type=str)
  parser.add_argument('--output_adata', type=str, default=None)

  parser = add_preprocess_args(parser)

  ARGS = parser.parse_args()
  adata = load_data(ARGS.dataset)
  adata = run_preprocess(adata, ARGS)

  if ARGS.output_adata is not None:
    print(f'Saving {adata.shape} --> {ARGS.output_adata}')
    adata.write(f'{ARGS.output_adata}')

