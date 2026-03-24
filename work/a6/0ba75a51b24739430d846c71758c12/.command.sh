#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 replicate_zero.py         --input_data input.csv         --file_name test         --quantile true         --log_file /root/projects/work/9e/71301416c4c6835abd91d5c7d8ace8/test_LOG.txt        --log_output "test_LOG.txt"         --output "test_replicate_zero.csv"         --pca_output "test_replicate_zero_pca"
