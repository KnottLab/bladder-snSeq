#!/usr/bin/env bash

filename=$(basename -- "${1}")
filename="${filename%.*}"
output="${filename}.clust.h5ad"

echo $output

scvi_.py $1 --output_adata ${output} \
  --latent 256 \
  --layers 1 \
  --n_epochs_kl_warmup 2 \
  --epochs 300 \
  --lr 0.001 \
  --vae_type vanilla \
  --batch_key batch_num \
  --reset_to_counts \
  --restore_all_genes \
  --estimate_cell_cycle 0 \
  --reconstruction_loss zinb

do_umap.py ${output} --output_adata ${output} --n_neighbors 30 --min_dist 0.1 --spread 3 --gpu
do_tsne.py ${output} --output_adata ${output} --learning_rate 100 --perplexity 10 --n_iter 2500 --gpu

cluster_.py ${output} --output_adata ${output} --latent_key scVI_vanilla --method leiden --resolution 0.6

