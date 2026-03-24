#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 start.py         --log_file /root/projects/work/fe/20e5e91cae966e299b70765faad74f/test_LOG.txt        --log_output "test_LOG.txt"        

Rscript imputation.r         --input_data /root/projects/work/fe/20e5e91cae966e299b70765faad74f/test_cv_icc.csv         --file_name test         --output "test_imputation.csv"   

python3 end.py         --file_name test         --log_output "test_LOG.txt"         --output "test_imputation.csv"         --pca_output "test_imputation_pca"
