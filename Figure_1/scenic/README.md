# SCENIC

```
pyscenic 0.10.0
```

First, download the annotation files to this directory:
  - `hs_hgnc_curated_tfs.txt` 
  - `motifs-v9-nr.hgnc-m0.001-o0.0.tbl` 
  - `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather`

Next, use the gene list of highly variable genes `bladder_mibc_scenic_genes.txt` to create a record of the snSeq data in `.loom` format.

Finally, run the pipeline which calls the pyscenic cli:

```bash
$ ./pyscenic_pipeline.sh $loom_file
```

Notice that lots of RAM is required to process large datasets so keeping counts as sparse might be desirable, however the `sparse` flag may also lead to errors in the pyscenic cli. There is a patch that allows using `sparse=True` detailed here: https://github.com/aertslab/pySCENIC/issues/153
