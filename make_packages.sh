#!/bin/bash

# Build BSgenome package
Rscript make_BSgenome_package.R -i config.yaml --seqdir seqs

# Build TxDb package
Rscript make_TxDb_package.R -i config.yaml --annotdir annot
