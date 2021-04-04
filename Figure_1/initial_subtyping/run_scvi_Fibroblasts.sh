#!/usr/bin/env bash

filename=$(basename -- "${1}")
filename="${filename%.*}"
output="${filename}.clust.h5ad"

echo $output

scvi_.py $1 --output_adata ${output} \
  --latent 16 \
  --layers 1 \
  --n_epochs_kl_warmup 1 \
  --epochs 200 \
  --lr 0.001 \
  --vae_type vanilla \
  --batch_key batch_num \
  --force_n_genes 1500 \
  --reset_to_counts \
  --restore_all_genes \
  --estimate_cell_cycle 0 \
  --reconstruction_loss zinb

do_umap.py ${output} --output_adata ${output} --n_neighbors 20 --min_dist 0.2 --spread 5 --gpu
do_tsne.py ${output} --output_adata ${output} --perplexity 20 --n_iter 2500 --early_exaggeration 100 --gpu

cluster_.py ${output} --output_adata ${output} --latent_key scVI_vanilla --method leiden --resolution 0.3

