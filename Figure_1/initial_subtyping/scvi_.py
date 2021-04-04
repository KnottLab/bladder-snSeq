#!/usr/bin/env python 

import numpy as np
# from scvi.models import VAE, SCANVI, LDVAE
# from scvi.inference import UnsupervisedTrainer
# from scvi.dataset import AnnDatasetFromAnnData

import scvi

import torch

import argparse
import sys
import os

from load_data import load_data
from preprocessing import run_preprocess, add_preprocess_args
import scanpy as sc

from anndata import AnnData
#from utils import makelogger, getlogger

from umap import UMAP
import logging

"""
Modules have two modes: standalone from command line and pipelined

Both modes accept a preprocessed AnnData object as input.

Standalone mode writes back a AnnData with new metadata

Pipelined mode returns the AnnData object with new metadata

scVI: https://github.com/YosefLab/scVI

NOTE scVI 's data sets use the python dataclass , introduced in Python 3.7

REF:
https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/basic_tutorial.ipynb

MANUAL:
https://readthedocs.org/projects/scvi/downloads/pdf/latest/

"""

logger = logging.getLogger('scvi')
logger.setLevel('INFO')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh = logging.StreamHandler()
sh.setFormatter(formatter)
logger.addHandler(sh)

## ------------------ Parameters
LR = 1e-4
N_EPOCHS = 500
USE_CUDA = True
VAE_TYPE = 'vanilla'


def reset_to_counts(adata, ARGS):

  # TODO detect logged data

  #keep_vars = np.array([v in adata.var_names for v in adata.raw.var_names])
  adata.X = adata.raw.X[:, adata.raw.var_names.isin(adata.var_names)]

  return adata


def run_scvi(adata, ARGS):

  batch_key = 'batch' if 'batch_num' not in adata.obs.keys() else 'batch_num'

  if ARGS.reset_to_counts:
    adata = reset_to_counts(adata, ARGS)

  #dataset = AnnDatasetFromAnnData(ad=adata, batch_label=ARGS.batch_key)
  scvi.data.setup_anndata(adata, batch_key=ARGS.batch_key)

  ## https://github.com/YosefLab/scVI/blob/master/scvi/models/vae.py 
  if ARGS.vae_type == 'vanilla':
    # vae = VAE(dataset.nb_genes, 
    model = scvi.model.SCVI(adata,
      #n_batch=dataset.n_batches,
      # n_batch=ARGS.n_batch,
      n_latent=ARGS.latent,
      n_hidden=ARGS.latent * 2,
      n_layers=ARGS.layers,
      dropout_rate=0.2,
      dispersion='gene',
      #log_variational=True,
      #reconstruction_loss=ARGS.reconstruction_loss
      gene_likelihood=ARGS.reconstruction_loss
    )

  elif ARGS.vae_type == 'linear':
    logger.info('using linear decoding scvi model')
    #vae = LDVAE(dataset.nb_genes, 
    #  # n_batch=dataset.n_batches,
    #  n_batch=ARGS.n_batch,
    #  n_latent=ARGS.latent,
    #  dispersion='gene')
    model = scvi.model.LinearSCVI(adata, n_latent=ARGS.n_latent)

  # elif ARGS.vae_type == 'linear':
  #   vae = SCANVI(dataset.nb_genes, 
  #     # n_batch=dataset.n_batches,
  #     n_batch=ARGS.n_batch,
  #     n_labels=ARGS.n_labels,
  #     n_layers=ARGS.layers,
  #     n_hidden=ARGS.latent * 2,
  #     n_latent=ARGS.latent,
  #     dropout_rate=0.1,
  #     dispersion='gene',
  #     log_variational=True,
  #     reconstruction_loss='nb')

  stopping_params = {'patience': 10, 'threshold': 0}
  stopping_params['early_stopping_metric'] = 'reconstruction_error'
  stopping_params['save_best_state_metric'] = 'reconstruction_error'

  #trainer = UnsupervisedTrainer(vae, dataset, 
  #  train_size=0.9, 
  #  n_epochs_kl_warmup=ARGS.n_epochs_kl_warmup, 
  #  metrics_to_monitor=['reconstruction_error'],
  #  frequency = 2, 
  #  early_stopping_kwargs = stopping_params,
  #  # n_iter_kl_warmup=int(128*5000/400),
  #  use_cuda=True)

  #trainer.history['reconstruction_error_test_set'].append(0)
  #trainer.train(n_epochs=ARGS.epochs, lr=ARGS.lr)

  # # from Solo drop learning rate and train some more
  # trainer.early_stopping.wait = 0
  # trainer.train(n_epochs = int(0.5*ARGS.epochs), lr= 0.5*ARGS.lr)

  model.train(lr = ARGS.lr,
    train_size = 0.9,
    n_epochs=ARGS.epochs,
    n_epochs_kl_warmup = ARGS.n_epochs_kl_warmup,
    frequency = 1, 
    metrics_to_monitor=['reconstruction_error'],
    early_stopping_kwargs = stopping_params,
    # n_iter_kl_warmup=int(128*5000/400),
    # use_cuda=True)
  )

  # trainer.history['reconstruction_error_test_set'].append(0)
  # trainer.train(n_epochs=ARGS.epochs, lr=ARGS.lr)

  # Stash the trainer for downstream analyses later - depricated
  # from Solo drop learning rate and train some more
  # trainer.early_stopping.wait = 0
  # trainer.train(n_epochs = 0.5*ARGS.epochs, lr= 0.5*ARGS.lr)

  # # Stash the trainer for downstream analyses later
  # if ARGS.save_scvi is not None:
  #   save_scvi = ARGS.save_scvi
  #   if not os.path.exists(save_scvi):
  #     os.makedirs(save_scvi)
  #   torch.save(trainer.model.state_dict(), f'{save_scvi}/vae.pkl')

  # full = trainer.create_posterior(trainer.model, dataset, indices = np.arange(len(dataset)))
  # latent, batch_indices, labels = full.sequential().get_latent()

  latent = model.get_latent_representation(ad)
  print('scvi finished.')
  # print('batch_indices', batch_indices.shape)
  # print('labels', labels.shape)

  # add latent space to adata and return
  adata.obsm[f'X_scVI_{ARGS.vae_type}'] = latent

  # # Umap the latent space and add it to multidimensional obs
  # if ARGS.do_umap:
  #   emb = UMAP(n_neighbors=50, min_dist=0.25).fit_transform(latent)
  #   adata.obsm[f'X_scVI_umap_{ARGS.vae_type}'] = emb

  #   # emb = UMAP(n_components=3, n_neighbors=20, min_dist=0.5).fit_transform(latent)
  #   # adata.obsm[f'X_scVI_umap3_{ARGS.vae_type}'] = emb

  if ARGS.save_dir is not None:
    logger.info(f'Saving to {ARGS.save_dir}')
    # full.save_posterior(ARGS.save_dir)
    model.save(ARGS.save_dir)
    with open(f'{ARGS.save_dir}/scvi_args.txt', 'w+') as f:
      for k, i in ARGS.__dict__.items():
        f.write(f'{k}: {i}\n')


  # Track the used genes for later
  adata.uns['scVI_genes'] = np.array(adata.var_names)

  # Add in scVI args to log
  adata.uns['scVI_args'] = ARGS.__dict__

  return adata


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('dataset', type=str)
  parser.add_argument('--save_scvi', type=str, default = None,
                      help = 'Path to save model weights')

  parser.add_argument('--output_adata', default=None, type=str,
                      help = 'Path to save adata with scvi run')

  parser.add_argument('--save_dir', default=None, type=str,
                      help = 'Path to a directory for saving the full posterior')

  parser.add_argument('--n_batch', default=1, type=int) # unused
  parser.add_argument('--latent', default=16, type=int)
  parser.add_argument('--layers', default=1, type=int)
  parser.add_argument('--lr', default=1e-4, type=float)
  parser.add_argument('--epochs', default=300, type=int)
  parser.add_argument('--n_epochs_kl_warmup', default=100, type=int)
  parser.add_argument('--vae_type', type=str, default = 'vanilla')
  parser.add_argument('--reconstruction_loss', type=str, default = 'nb')

  parser.add_argument('--reset_to_counts', action = 'store_true',
                      help = 'Whether to force reset adata.X to adata.raw.X \
                      presumably containing raw counts data. ')

  parser.add_argument('--restore_all_genes', action = 'store_true',
                      help = 'Whether to restore all genes to the default \
                      layer (in raw count form) after scVI is finished.')

  # parser.add_argument('--do_umap', action = 'store_true')

  parser = add_preprocess_args(parser)


  ARGS = parser.parse_args()
  if os.path.exists(ARGS.save_dir):
    logger.warn('Provide save_dir that does not already exist')
    sys.exit(1)

  adata = load_data(ARGS.dataset)
  adata.raw = adata

  # Detect preprocessing run on this dataset
  if (not 'PREPROCESSED_FLAG' in adata.uns.keys()) or (not adata.uns['PREPROCESSED_FLAG']):
    logger.info('Detected dataset that has not been preprocessed.')
    adata = run_preprocess(adata, ARGS)

    # # Cut out empty cells
    # sc.pp.filter_cells(adata, min_counts=50)
    # print(f'[SCVI] filtered out empty ( < 50 counts ) cells. {adata.shape}')


  logger.info(f'Starting scVI with adata: {adata.shape}')
  adata = run_scvi(adata, ARGS)

  # Reset to raw counts (not subsetted) but carry over obs and obsm
  if ARGS.restore_all_genes:
    adata = AnnData(adata.raw.X, obs=adata.obs, obsm=adata.obsm, var=adata.raw.var, varm=adata.varm, uns=adata.uns)

  if ARGS.output_adata is not None:
    # assert os.path.exists(os.path.dirname(ARGS.output_adata))
    logger.info(f'Writing adata ({adata.shape}) --> {ARGS.output_adata}')
    adata.write(ARGS.output_adata)
