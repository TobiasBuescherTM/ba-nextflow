#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env

python3 start.py         --log_file /root/projects/work/21/fa2818def01fad6e37226ba419ee13/test_LOG.txt        --log_output "test_LOG.txt"        

Rscript vsn.r         --input_data /root/projects/work/21/fa2818def01fad6e37226ba419ee13/test_imputation.csv         --file_name test         --output "test_vsn.csv"   

python3 end.py         --file_name test         --log_output "test_LOG.txt"         --output "test_vsn.csv"         --pca_output "test_vsn_pca"
