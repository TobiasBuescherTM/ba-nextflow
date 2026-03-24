#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 start.py         --log_file /root/projects/work/1c/970641cd6af452fdf79d536e096d62/test_LOG.txt        --log_output "test_LOG.txt"        

Rscript imputation.r         --input_data /root/projects/work/1c/970641cd6af452fdf79d536e096d62/test_cv_icc.csv         --file_name test         --output "test_imputation.csv"   

python3 end.py         --file_name test         --log_output "test_LOG.txt"         --output "test_imputation.csv"         --pca_output "test_imputation_pca"
