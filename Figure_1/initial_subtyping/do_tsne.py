#!/usr/bin/env python

import numpy as np
import argparse

from load_data import load_data
from MulticoreTSNE import MulticoreTSNE as TSNE

try:
  import cuml
  CUML_FLAG=True
except:
  print('[DO_TSNE] WARNING failed to import cuML. GPU accelerated TSNE will not be available.')
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
def do_tsne(adata, ARGS):
  latent = adata.obsm[ARGS.latent_key]

  if ARGS.gpu and CUML_FLAG:
    print('[DO_TSNE] Using cuML GPU-accelerated TSNE')
    umap_class = cuml.UMAP
    if ARGS.metric != 'euclidean':
      print('[DO_TSNE] cuML TSNE requres euclidean distance metric.')

    emb = cuml.TSNE(
      perplexity = ARGS.perplexity,
      learning_rate = ARGS.learning_rate,
      early_exaggeration = ARGS.early_exaggeration,
      ).fit_transform(latent)
  else:
    print('[DO_TSNE] Using MulticoreTSNE')
    emb = TSNE( perplexity = ARGS.perplexity,
      metric = ARGS.metric,
      verbose = False, n_jobs=ARGS.n_jobs).fit_transform(latent)
  
  print(f'[DO_TSNE] placing embedding {emb.shape} in key {ARGS.tsne_key}')
  adata.obsm[ARGS.tsne_key] = emb

  print(f'[DO_TSNE] recording tSNE args')
  adata.uns['tSNE_args'] = ARGS.__dict__

  return adata


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('dataset', type=str)

  parser.add_argument('--latent_key', default='X_scVI_vanilla', type=str,
                      help = 'Key in adata.obsm to use as features for tsne.')
                      
  parser.add_argument('--tsne_key', default='X_scVI_tsne_vanilla', type=str,
                      help = 'Key in adata.obsm to save tsne embedding.')

  parser.add_argument('--gpu', action='store_true', 
                      help = 'Whether to use GPU-accelerated tsne via RapidsAI \
                      and the cuML library. ')

  parser.add_argument('-j', '--n_jobs', default=12, type=int, 
                      help = 'Number of jobs for MulticoreTSNE')

  parser.add_argument('--perplexity', default=20, type=int)
  parser.add_argument('--learning_rate', default=200., type=float)
  parser.add_argument('--n_iter', default=1000, type=int)
  parser.add_argument('--metric', default='euclidean', type=str)
  parser.add_argument('--early_exaggeration', default=12, type=float)

  parser.add_argument('--output_adata', default=None, type=str,
                      help = 'Path to save.')

  ARGS = parser.parse_args()
  adata = load_data(ARGS.dataset)

  adata = do_tsne(adata, ARGS)
  if ARGS.output_adata is not None:
    print(f'[DO_TSNE] Writing to {ARGS.output_adata}')
    adata.write(ARGS.output_adata)

