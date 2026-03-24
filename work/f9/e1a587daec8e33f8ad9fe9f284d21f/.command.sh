#!/bin/bash -ue
eval "innit.pymicromamba shell hook --shell bash)"
micromamba activate nextflow_env

python3 input.csv         --input_data test         --file_name test         --log_output "test_LOG.txt"        --pca-output "_innit_pca"
