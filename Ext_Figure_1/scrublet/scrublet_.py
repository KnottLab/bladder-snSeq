#!/usr/bin/env python

import numpy as np
import scanpy as sc
import scrublet as scr
from anndata import AnnData

import argparse
import os

from load_data import load_data
from preprocessing import run_preprocess, add_preprocess_args

"""
Modules have two modes: standalone from command line and pipelined

Both mode accepts a preprocessed AnnData object as input.
  - Standalone mode writes back a AnnData with new metadata
  - Pipelined mode returns the AnnData object with new metadata

REF:
https://github.com/AllonKleinLab/scrublet

Wrapper around scrublet, adding in a progressive doublet removal
step.

Set --rounds 1 to run the default doublet detection.

N.Ing <Nathan.Ing@cshs.org, ing.nathany@gmail.com> 

"""

def run_scrublet(adata, logfile, ARGS):
  # Scrublet works with an un-normalized count matrix
  print(f'[SCRUBLET] WARNING converting {adata.shape[0]} cells and {adata.var.highly_variable.sum()} genes to dense', file=logfile)

  mat = adata.X[:, adata.var['highly_variable'].values].toarray()
  print(mat.shape, file=logfile)

  percent = ARGS.percent / ARGS.rounds
  cumulative_doublets = np.zeros(mat.shape[0], dtype=np.bool)
  for sc_round in range(ARGS.rounds):
    # Vector thats going to stay the same size as the number cells
    cumulative_scores = np.zeros_like(cumulative_doublets).astype(np.float32)
    cumulative_scores[:] = np.nan
    this_round_doublets = np.zeros_like(cumulative_scores)
    this_round_doublets[:] = np.nan

    print(f'Running scrublet round {sc_round}', file=logfile)
    # Expected doublet rate ~ 7.5% for 10X at full saturation
    scrubber = scr.Scrublet(mat, sim_doublet_ratio=3, 
      expected_doublet_rate=ARGS.percent, 
      stdev_doublet_rate=0.025)
      

    doublet_scores, predicted_doublets = scrubber.scrub_doublets(min_counts=2, 
      min_cells=10, 
      min_gene_variability_pctl=0,  # How many genes get used in the classifier. Higher = Fewer. Use 0 if HVG already selected
      n_prin_comps=30, 
      log_transform=False,
      mean_center=True,
      synthetic_doublet_umi_subsampling=0.85,
      use_approx_neighbors=False,
      verbose=True)

    # Re-do the thresholding
    predicted_doublets = doublet_scores > 0.25

    # Add scrublet predictions to metadata; subset matrix
    cumulative_scores[np.logical_not(cumulative_doublets)] = doublet_scores
    this_round_doublets[np.logical_not(cumulative_doublets)] = predicted_doublets

    cumulative_doublets[np.logical_not(cumulative_doublets)] = predicted_doublets
    # adata.obs[f'scrublet_scores_{sc_round}'] = cumulative_scores
    # adata.obs[f'scrublet_calls_{sc_round}'] = this_round_doublets

    print(f'Called {np.sum(predicted_doublets)} doublets', file=logfile)
    print(f'Subsetting {mat.shape}', file=logfile)
    mat = mat[np.logical_not(predicted_doublets), :]
    print(f'New matrix {mat.shape}', file=logfile)

  adata.obs['scrublet'] = cumulative_doublets
  adata.obs['scrublet_score'] = doublet_scores

  print('[SCRUBLET] Scrublet run stats:', file=logfile)
  # for sc_round in range(ARGS.rounds):
  #   n_calls = np.nansum(adata.obs[f'scrublet_calls_{sc_round}'].values).astype(np.int)
  #   mn_score = np.nanmean(adata.obs[f'scrublet_scores_{sc_round}'].values)
  #   print(f'Round {sc_round} calls: {n_calls: 3d} mean score: {mn_score:3.4f}', file=logfile)

  print('[SCRUBLET] Total calls', np.nansum(adata.obs['scrublet'].values), file=logfile)
  print('[SCRUBLET] Total calls', np.nansum(adata.obs['scrublet'].values))

  # Add a UMI precentile cutoff
  if ARGS.use_umi_cutoff:
    umi_cutoff = np.quantile(adata.obs['total_counts'].values, 0.95)
    passing_counts = adata.obs['total_counts'].values < umi_cutoff
    print(f'[SCRUBLET] UMI threshold passed by {passing_counts.sum()} cells', file=logfile)
  else:
    print('[SCRUBLET] Not applying upper UMI threshold', file=logfile)
    passing_counts = np.ones(adata.shape[0], dtype=bool)

  combined_gates = np.logical_not(cumulative_doublets) * passing_counts

  print(f'[SCRUBLET] Combined scrublet + UMI gates passed by {combined_gates.sum()} cells', file=logfile)
  barcodes = adata.obs_names.values[combined_gates]

  return barcodes


if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('dataset', type=str)
  parser.add_argument('--rounds', type=int, default=1)
  parser.add_argument('--percent', type=float, default=0.1, 
    help = 'The expected overall doublet rate for this sample.')
  parser.add_argument('--use_umi_cutoff', action='store_true')

  # parser.add_argument('--output_barcodes', default=None, type=str)
  parser.add_argument('--output_adata', default=None, type=str)
  parser.add_argument('--log', default='/dev/null', type=str)

  parser = add_preprocess_args(parser)

  ARGS = parser.parse_args()
  adata = load_data(ARGS.dataset)
  adata.var_names_make_unique()
  print(f'[SCRUBLET] loaded adata {adata.shape}')

  # # Detect preprocessing run on this dataset
  # if not 'PREPROCESSED_FLAG' in adata.uns.keys():
  #   print('Detected dataset that has not been preprocessed.')
  #   adata_pp = run_preprocess(adata, ARGS=ARGS)
  # else:
  #   adata_pp = adata
  # adata

  ## We dont want to do anything to counts but we want to identify highly variable genes 
  adata.raw = adata
  sc.pp.calculate_qc_metrics(adata, inplace=True) # for count thresholding later
  sc.pp.normalize_total(adata, target_sum=10000)
  sc.pp.log1p(adata)
  sc.pp.highly_variable_genes(adata, n_top_genes=int(adata.shape[-1] * 0.1) )

  adata = AnnData(adata.raw.X, obs=adata.obs, var=adata.var)


  with open(ARGS.log, 'w+') as logfile:
    logfile.write('\nScrublet\n')
    keep_barcodes = run_scrublet(adata, logfile, ARGS)

  # if ARGS.output_barcodes is not None:
  #   doublet_barcodes = adata.obs.index.values[adata.obs['scrublet']]
  #   print(f'Writing {len(doublet_barcodes)} doublet barcodes to {ARGS.output_barcodes}')
  #   with open(ARGS.output_barcodes, 'w+') as f:
  #     for b in doublet_barcodes:
  #       f.write(f'{b}\n')

  if ARGS.output_adata is not None:
    # assert os.path.exists(os.path.dirname(ARGS.output_adata))
    # adata = adata[np.logical_not(adata.obs['scrublet'].values), :]

    all_cell_output = ARGS.output_adata.replace('.h5ad', '_all_cells.h5ad')
    print(f'Writing adata with all cells ({adata.shape}) to {all_cell_output}')
    adata.write(all_cell_output)

    adata = adata[keep_barcodes, :]
    print(f'Writing adata ({adata.shape}) to {ARGS.output_adata}')
    adata.write(ARGS.output_adata)
