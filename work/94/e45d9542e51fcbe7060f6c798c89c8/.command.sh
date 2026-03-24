#!/bin/bash -ue
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow_env
python3 80%.py         --input_data /root/projects/work/71/a204a2204a1e710cbd93a75e128437/test_qcomics.csv         --file_name test         --log_file /root/projects/work/71/a204a2204a1e710cbd93a75e128437/test_LOG.txt        --log_output "test_LOG.txt"         --output "test_80%.csv"         --pca_output "test_80%_pca"
