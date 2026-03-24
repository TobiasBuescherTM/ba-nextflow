#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 qcomics.py         --input_data /root/projects/work/a6/0ba75a51b24739430d846c71758c12/test_replicate_zero.csv         --file_name test         --log_file /root/projects/work/a6/0ba75a51b24739430d846c71758c12/test_LOG.txt        --log_output "test_LOG.txt"         --output "test_qcomics.csv"         --pca_output "test_qcomics_pca"
