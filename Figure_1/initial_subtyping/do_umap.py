#!/usr/bin/env python

import numpy as np
from umap import UMAP
import argparse

from load_data import load_data

try:
  import cuml
  CUML_FLAG=True
except:
  print('[DO_UMAP] WARNING failed to import cuML. GPU accelerated UMAP will not be available.')
  CUML_FLAG=False

"""
Modules have two modes: standalone from command line and pipelined

Both modes accept a preprocessed AnnData object as input.

Standalone mode writes back a AnnData with new metadata

Pipelined mode returns the AnnData object with new metadata

UMAPs with /umap-learn/cuML GPU-accelerated UMAP implementation

https://umap-learn.readthedocs.io/en/latest/

https://github.com/lmcinnes/umap

"""
def do_umap(adata, ARGS):
  latent = adata.obsm[ARGS.latent_key]

  if ARGS.gpu and CUML_FLAG:
    umap_class = cuml.UMAP
    if ARGS.metric != 'euclidean':
      print('[DO_UMAP] cuML UMAP requres euclidean distance metric.')

    emb = cuml.UMAP(
      n_neighbors = ARGS.n_neighbors,
      min_dist = ARGS.min_dist,
      spread = ARGS.spread
      ).fit_transform(latent)
  else:
    emb = UMAP(
      n_neighbors = ARGS.n_neighbors,
      min_dist = ARGS.min_dist,
      metric = ARGS.metric,
      verbose = True).fit_transform(latent)
  
  print(f'[DO_UMAP] placing embedding {emb.shape} in key {ARGS.umap_key}')
  adata.obsm[ARGS.umap_key] = emb

  print(f'[DO_UMAP] recording UMAP args')
  adata.uns['UMAP_args'] = ARGS.__dict__

  return adata


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('dataset', type=str)

  parser.add_argument('--latent_key', default='X_scVI_vanilla', type=str,
                      help = 'Key in adata.obsm to use as features for umap.')
                      
  parser.add_argument('--umap_key', default='X_scVI_umap_vanilla', type=str,
                      help = 'Key in adata.obsm to save umap embedding.')

  parser.add_argument('--gpu', action='store_true', 
                      help = 'Whether to use GPU-accelerated UMAP via RapidsAI \
                      and the cuML library. ')

  parser.add_argument('--n_neighbors', default=25, type=int)
  parser.add_argument('--min_dist', default=0.5, type=float)
  parser.add_argument('--spread', default=0.5, type=float)
  parser.add_argument('--metric', default='euclidean', type=str)

  parser.add_argument('--output_adata', default=None, type=str,
                      help = 'Path to save.')

  ARGS = parser.parse_args()
  adata = load_data(ARGS.dataset)

  adata = do_umap(adata, ARGS)
  if ARGS.output_adata is not None:
    print(f'[DO_UMAP] Writing to {ARGS.output_adata}')
    adata.write(ARGS.output_adata)


"""
UMAP metrics
------------

https://umap-learn.readthedocs.io/en/latest/parameters.html 

Minkowski style metrics
  euclidean
  manhattan
  chebyshev
  minkowski

Miscellaneous spatial metrics
  canberra
  braycurtis
  haversine

Normalized spatial metrics
  mahalanobis
  wminkowski
  seuclidean

Angular and correlation metrics
  cosine
  correlation

Metrics for binary data
  hamming
  jaccard
  dice
  russellrao
  kulsinski
  rogerstanimoto
  sokalmichener
  sokalsneath
  yule
  """