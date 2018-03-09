#!/bin/bash

# Build BSgenome package
Rscript construct_BSgenome_package.R -i config.yaml --seqdir seqs

# Build TxDb package
Rscript construct_TxDb_package.R -i config.yaml --annotdir annot
