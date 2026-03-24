#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env

python3 innit.py         --input_data input.csv         --file_name test         --log_output "test_LOG.txt"        --pca-output "test_innit_pca"
