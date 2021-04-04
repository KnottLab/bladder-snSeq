#!/usr/bin/env bash

set -e

workers=12
input=$1

database_dir=$PWD

pyscenic grn --num_workers ${workers} \
  --sparse \
  --gene_attribute var_names \
  --cell_id_attribute obs_names \
  -o ${input}.adj.csv \
  ${input} \
  ${database_dir}/hs_hgnc_curated_tfs.txt

pyscenic ctx --num_workers ${workers} \
  --sparse \
  --gene_attribute var_names \
  --cell_id_attribute obs_names \
  -o ${input}.regulons.csv \
  --expression_mtx_fname ${input} \
  --all_modules \
  --mask_dropouts \
  --thresholds 0.25 0.5 0.75 0.9 \
  --top_n_targets 100 \
  --min_genes 10 \
  --annotations_fname ${database_dir}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  ${input}.adj.csv \
  ${database_dir}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather  

pyscenic aucell --num_workers ${workers} \
  --sparse \
  --gene_attribute var_names \
  --cell_id_attribute obs_names \
  -o ${input}.aucell.csv \
  ${input} \
  ${input}.regulons.csv
