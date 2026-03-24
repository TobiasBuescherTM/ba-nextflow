#!/bin/bash -ue
python3 start.py         --log_file /root/projects/work/c4/8239033891fe6c01ff97e25186c3ce/test_LOG.txt        --log_output "test_LOG.txt"        

Rscript imputation.r         --input_data /root/projects/work/c4/8239033891fe6c01ff97e25186c3ce/test_cv_icc.csv         --file_name test         --output "test_imputation.csv"   

python3 end.py         --file_name test         --log_output "test_LOG.txt"         --output "test_imputation.csv"         --pca_output "test_imputation_pca"
