#!/usr/bin/env bash

parallel --link -j 4 scrublet_.py {1} --output_adata scrublet/{2}.h5ad --percent 0.07 --rounds 1 --log scrublet/{2}.log :::: files.txt :::: names.txt

