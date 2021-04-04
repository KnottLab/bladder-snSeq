#!/usr/bin/env python

import argparse
from load_data import load_data

from joblib import parallel_backend

import numpy as np
import pandas as pd
import sys
import os

import scanpy as sc
from sklearn.cluster import MiniBatchKMeans, DBSCAN

"""
Modules have two modes: standalone from command line and pipelined

Both mode accepts a preprocessed AnnData object as input.
  - Standalone mode writes back a AnnData with new metadata
  - Pipelined mode returns the AnnData object with new metadata

Run PCA - based clustering with various clustering algorithms.


N.Ing <Nathan.Ing@cshs.org, ing.nathany@gmail.com> 

"""
def find_clusters(adata, method, latent_key, ARGS): 
  cluster_key = f'{latent_key}_{method}'
  
  if cluster_key in adata.obsm.keys():
    print(f'Key {cluster_key} already in adata.obsm')
    sys.exit(1)

  # Extract the latent vectors.
  latent_key = f'X_{latent_key}'
  try:
    latent = adata.obsm[latent_key]
    print(f'Loaded latent embedding of {latent.shape} from {latent_key}')
  except:
    ## Cant do this here anymore because one cannot assume the dataset has been preprocessed for PCA
    # if latent_key == 'X_pca':
    #   """
    #   https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.api.pp.pca.html

    #   scanpy.api.pp.pca(data, n_comps=50, zero_center=True, svd_solver='auto', 
    #     random_state=0, return_info=False, use_highly_variable=None, dtype='float32', 
    #     copy=False, chunked=False, chunk_size=None)
    #   """
    #   sc.pp.pca(adata)
    #   latent = adata.obsm[latent_key]
    # else:
    raise Exception(f'Error clustering. key {latent_key} not found in {adata.obsm.keys()}.')


  print(f'Running clustering {method} and adding to adata.obs["{cluster_key}"]')
  if method == 'kmeans':
    """
    https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html
    class sklearn.cluster.MiniBatchKMeans(n_clusters=8, init='k-means++', max_iter=100, 
      batch_size=100, verbose=0, compute_labels=True, random_state=None, tol=0.0, 
      max_no_improvement=10, init_size=None, n_init=3, reassignment_ratio=0.01)

    """
    clusterer = MiniBatchKMeans(n_clusters = ARGS.k, 
                                batch_size = 256,    # This might become smaller than n_samples then we got problems.
                                init_size = min(256 * 4, adata.shape[0]),
                                n_init = 25,
                                reassignment_ratio = 0.05,
                                verbose = False).fit(latent)
    clusters = clusterer.labels_
    adata.obs[cluster_key] = pd.Series(clusters.astype(np.str), index=adata.obs.index, dtype='category')

  elif method == 'dbscan':
    """
    class sklearn.cluster.DBSCAN(eps=0.5, min_samples=5, metric='euclidean', 
      metric_params=None, algorithm='auto', leaf_size=30, p=None, n_jobs=None)
    """
    clusterer = DBSCAN(eps = 1, n_jobs = ARGS.j).fit(latent)
    clusters = clusterer.labels_
    adata.obs[cluster_key] = pd.Series(clusters.astype(np.int), index=adata.obs.index, dtype='category')

  elif method == 'louvain':
    """
    Discussion on running neighbors() and UMAP with a parallel backend: 
    https://github.com/theislab/scanpy/issues/913
    https://github.com/lmcinnes/umap/issues/317 

    First:
    scanpy.pp.neighbors(adata, n_neighbors=15, n_pcs=None, use_rep=None, knn=True, 
      random_state=0, method='umap', metric='euclidean', metric_kwds={}, copy=False)

    Then:
    scanpy.api.tl.louvain(adata, resolution=None, random_state=0, restrict_to=None, 
      key_added='louvain', adjacency=None, flavor='vtraag', directed=True, 
      use_weights=False, partition_type=None, partition_kwargs=None, copy=False)
    """
    # with parallel_backend('threading', n_jobs=ARGS.j):
    sc.pp.neighbors(adata, use_rep=latent_key)
    sc.tl.louvain(adata, resolution=ARGS.resolution, key_added=cluster_key)

    # Pad single digits with 0 so they come first in alpha-order
    clusters = adata.obs[cluster_key].values
    c_padded = [f'{int(c):02d}' for c in clusters]
    adata.obs[cluster_key] = c_padded

  elif method == 'leiden':
    """
    scanpy.tl.leiden(adata, resolution=1, *, restrict_to=None, random_state=0, key_added='leiden', 
      adjacency=None, directed=True, use_weights=True, n_iterations=-1, partition_type=None, 
      copy=False, **partition_kwargs)
    """
    # with parallel_backend('threading', n_jobs=ARGS.j):
    sc.pp.neighbors(adata, use_rep=latent_key)
    sc.tl.leiden(adata, resolution=ARGS.resolution, key_added=cluster_key)

    # Pad single digits with 0 so they come first in alpha-order
    clusters = adata.obs[cluster_key].values
    c_padded = [f'{int(c):02d}' for c in clusters]
    adata.obs[cluster_key] = c_padded

  else:
    raise Exception(f'Clustering method.')

  return adata


if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('dataset', type=str)
  parser.add_argument('--output_adata', default=None, type=str)
  parser.add_argument('-k', default=10, type=int,
                      help = 'Number of initial clusters. If None, use louvain to estimate.')
  parser.add_argument('-j', default=16, type=int,
                      help = 'Number of jobs for methods that can be run in parallel.')

  parser.add_argument('--latent_key', type=str, default='pca', nargs='+',
    help = 'The obs to be used for cluster identification. Must be a key of adata.obsm. \
            An "X_" is automatrically prepended to the key. \
            commonly pca (scanpy) , scVI_vanilla, scVI_linear (scVI), Embeded_z0.8 (DESC)')

  parser.add_argument('--method', type=str, default=['kmeans'], nargs='+',
    help = 'Clustering method to be used. \n\
            kmeans (sklearn), dbscan (sklearn), louvain (scanpy), leiden (scanpy)')

  parser.add_argument('--cluster_key', type=str, default=None,
    help = 'Key to save the cluster assignments. Defaults to "<LATENT_KEY>_<METHOD>"')

  parser.add_argument('--resolution', type=float, default=0.8,
    help = 'Setting for louvain and leiden clustering.')

  parser.add_argument('--dump_csv', action='store_true',
                      help = 'Dump the obs and obsm variables to CSV following clustering.')

  ARGS = parser.parse_args()

  adata = load_data(ARGS.dataset)

  for method in ARGS.method:
    for latent_key in ARGS.latent_key:
      adata = find_clusters(adata, method, latent_key, ARGS)

  if ARGS.output_adata is not None:
    print(f'Writing adata to {ARGS.output_adata}')
    adata.write(ARGS.output_adata)
