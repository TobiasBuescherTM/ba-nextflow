#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 80%.py         --input_data /root/projects/work/35/effd4fd016ac3ac2ce54f1feb6342e/test_qcomics.csv         --file_name test         --log_file /root/projects/work/35/effd4fd016ac3ac2ce54f1feb6342e/test_LOG.txt        --log_output "test_LOG.txt"         --output "test_80%.csv"         --pca_output "test_80%_pca"
