#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 cv_icc.py         --input_data /root/projects/work/24/75517729e6de6d1715dca6696e3dc0/test_80%.csv         --file_name test         --cv_cutoff 0.3        --icc_cutoff 0.65        --log_file /root/projects/work/24/75517729e6de6d1715dca6696e3dc0/test_LOG.txt        --log_output "test_LOG.txt"         --output "test_cv_icc.csv"         --pca_output "test_cv_icc_pca"         --correlations_path "test_cv_correlations"
