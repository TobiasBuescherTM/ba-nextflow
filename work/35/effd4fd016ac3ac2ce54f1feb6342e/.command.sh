#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 qcomics.py         --input_data /root/projects/work/f7/024c43621f5ad616e88024e298bb2f/test_replicate_zero.csv         --file_name test         --log_file /root/projects/work/f7/024c43621f5ad616e88024e298bb2f/test_LOG.txt        --log_output "test_LOG.txt"         --output "test_qcomics.csv"         --pca_output "test_qcomics_pca"
