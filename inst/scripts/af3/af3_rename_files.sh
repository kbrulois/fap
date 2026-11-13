#!/bin/bash

ml R/4.1.2


exec_path=$(Rscript -e 'cat(system.file("exec", package = "fap"))')

Rscript "$exec_path/af3_rename_files_cli.R" \
--input=$1 \
--name=$2 \
--person=SE \
--location=SD
