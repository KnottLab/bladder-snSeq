#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

import os
from glob import glob
import argparse
import logging

import ray
# import multiprocessing

from scipy.sparse import csr_matrix
from score_genelists_def import score_genelist


'''
Score lists of genes 
--------------------

Read in AnnData and reset to raw counts if possible
subset cells if needed

Read in lists 
If the list has just one member, skip it
Otherwise run scoring and create a table for single cells
'''

parser = argparse.ArgumentParser()
parser.add_argument(
  'genelist_dir', 
  help='A directory containing gene lists in newline delimited *.txt format'
)
parser.add_argument(
  'adata_path',
  help='AnnData object holding gene expression'
)
parser.add_argument(
  'groupby',
  help='A column in adata.obs to use as cell groups for background gene expression'
)
parser.add_argument(
  '--out', default=None, type=str,
  help='an h5ad file to stash results. if not specified, it defaults to converting the '+\
       'genelist_dir path into a file name'
)
parser.add_argument(
  '-j', '--n_jobs', default=8, type=int,
  help='Number of parallel jobs to launch. default=8'
)
ARGS = parser.parse_args()


logger = logging.getLogger('RECEPTOR SCORE')
logger.setLevel('INFO')
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


for k, v in ARGS.__dict__.items():
  logger.info(f'{k}:\t{v}')


# // init 
ray.init(num_cpus=ARGS.n_jobs)



# // load data
adata = sc.read_h5ad(ARGS.adata_path)

sc.pp.filter_genes(adata, min_cells = 10)
logger.info(f'Reduced by thresholding the number of genes: {adata.shape}')

logger.info('Apply preprocessing: count normalizing to 10,000 and logging counts')
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

gene_lists = sorted(glob(f'{ARGS.genelist_dir}/*.txt'))
logger.info(f'Loaded {len(gene_lists)} gene lists')

cell_groups = np.array(adata.obs[ARGS.groupby])



# // make data available to multiple workers
logger.info(f'PUTting var_names, gex, and cell_groups')
var_names = adata.var_names
var_names_id = ray.put(var_names)
gex_id = ray.put(adata.X)
cell_groups_id = ray.put(cell_groups)


# p = multiprocessing.Pool(args.j)
# ret = p.map(score_genelist, gene_lists)
# p.close()
# p.join()

# defines the compute tasks
futures = [score_genelist.remote(gex_id, var_names_id, cell_groups_id, gene_list) for gene_list in gene_lists]

# this little guy triggers the big compute
ret = ray.get(futures)

logger.info(f'Returned {len(ret)}')

list_scores = {}
for d in ret:
  if d is None:
    continue
  list_scores.update(d)

x = []
cols = []
for k,v in list_scores.items():
  logger.info(f'{k}: {v.shape}')
  x.append(v)
  cols.append(k)

x = np.stack(x, axis=1)
logger.info(f'Creating scores DataFrame from x={x.shape}')

list_scores = pd.DataFrame(x, index=adata.obs_names, columns=cols)
list_scores.index.name = 'barcodes'
logger.info(f'Created cell receptor score matrix: {list_scores.shape}')
logger.info(f'Converting to AnnData')

# TODO add descriptive information in uns
rscores = AnnData(csr_matrix(list_scores.values), 
                  obs = adata.obs,
                  var = pd.DataFrame(index = list_scores.columns),
                  uns = {'receptor_scoring_info': f'Receptor scoring'}
                  )


if ARGS.out is None:
  gld = ARGS.genelist_dir
  if gld.endswith('/'):
    gld = gld[:-1]
  ad_out = gld + '_scores.h5ad'
else:
  ad_out = ARGS.out

logger.info(f'Writing to {ad_out}')
# list_scores.to_csv(df_out, float_format='%.5f')
rscores.write(ad_out)
