#!/bin/bash -ue
python3 start.py         --log_file /root/projects/work/98/bc142d9e7285949f0dacc63258c557/test_LOG.txt        --log_output "test_LOG.txt"        

Rscript vsn.r         --input_data /root/projects/work/98/bc142d9e7285949f0dacc63258c557/test_imputation.csv         --file_name test         --output "test_vsn.csv"   

python3 end.py         --file_name test         --log_output "test_LOG.txt"         --output "test_vsn.csv"         --pca_output "test_vsn_pca"
