#!/usr/bin/env python

import numpy as np
import scanpy as sc
import argparse
import sys
import os

from os.path import splitext


def detect_format(pth):
  assert os.path.exists(pth)

  if os.path.isdir(pth):
    # Scanpy reading 10x data expects gzipped files
    if os.path.exists(f'{pth}/matrix.mtx.gz'):
      return '10X_MTX'

  ex = splitext(pth)[-1]
  print(f'[LOAD_DATA] Checking extension {ex}')
  
  if ex == '.h5ad':
    return 'H5AD'
  elif ex == '.h5':
    return '10X_H5'
  elif ex == '.loom':
    return 'loom'
  
  # We've reached the end of known file types
  return None


def load_data(pth):
  """
  Detect format and load dataset

  Formats: 
    - 10X h5 
    - 10X mtx
    - AnnData (`h5ad`)

  """
  assert os.path.exists(pth), f'WARNING: {pth}'
  data_format = detect_format(pth)

  if data_format == 'H5AD':
    adata = sc.read_h5ad(pth)

  elif data_format == '10X_H5':
    adata = sc.read_10x_h5(pth)

  elif data_format == '10X_MTX':
    adata = sc.read_10x_mtx(pth)

  elif data_format == 'loom':
    adata = sc.read_loom(pth)
  
  else:
    print('No dataset loaded')
    sys.exit(1)

  original_varnames = np.array(adata.var_names)
  print('stashing original var names and forcing uniqueness')
  adata.uns['original_varnames'] = original_varnames
  adata.var_names_make_unique()

  print(f'[LOAD_DATA] loaded adata: {adata.shape} with {len(adata.obs.columns)} obs')
  # print(adata)

  # Tag the dataset with the original path
  # Even if the filesystem is meaningless, its better than nothing
  # if 'original_path' in adata.uns.keys():
  #   adata.uns['original_path'] = pth

  # Strip the -1 from barcodes -- why
  # barcodes = [x.replace('-1', '') for x in adata.obs.index.values]
  # adata.obs.index = barcodes

  return adata


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('dataset', type=str)
  parser.add_argument('--output_adata', default=None, type=str)
  parser.add_argument('--doublets', default=None, type=str)

  ARGS = parser.parse_args()

  adata = load_data(ARGS.dataset)

  if ARGS.doublets is not None:
    adata, doublet_adata = remove_doublets(adata, ARGS.doublets)
    doublet_file = f'{os.path.splitext(ARGS.output_adata)[0]}.doublets.h5ad'
    print(f'Saving doublets as AnnData {doublet_adata.shape} --> {doublet_file}')
    doublet_adata.write(doublet_file)

  if ARGS.output_adata is not None:
    print(f'Saving raw counts as AnnData {adata.shape} --> {ARGS.output_adata}')
    adata.write(ARGS.output_adata)
