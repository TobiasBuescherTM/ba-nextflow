process imputation {

    container  "tobiasbueschertm/nf-imputation"

    input:
    val input_data
    val log_file
    path script
    val file_name
    path start_script
    path end_script
    

    
    output:
    path "${file_name}_LOG.txt" , emit: LOG
    path "${file_name}_imputation.csv", emit: df
    path "${file_name}_imputation_pca/*", emit: pca

    script:
    """

    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate nextflow_env
    python3 ${start_script} \
        --log_file ${log_file}\
        --log_output "${file_name}_LOG.txt" \
       
    
    Rscript ${script} \
        --input_data ${input_data} \
        --file_name ${file_name} \
        --output "${file_name}_imputation.csv"   

    python3 ${end_script} \
        --file_name ${file_name} \
        --log_output "${file_name}_LOG.txt" \
        --output "${file_name}_imputation.csv" \
        --pca_output "${file_name}_imputation_pca"

    """
}
