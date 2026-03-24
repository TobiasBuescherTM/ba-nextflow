process replicate_zero {

    container  "tobiasbueschertm/nf-replicate_zero"

    input:
    path input_data
    val log_file
    path script
    val file_name
    val quantile

    
    output:
    path "${file_name}_LOG.txt" , emit: LOG
    path "${file_name}_replicate_zero.csv", emit: df
    path "${file_name}_replicate_zero_pca/*", emit: pca

    script:
    """
    
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate nextflow_env
    python3 ${script} \
        --input_data ${input_data} \
        --file_name ${file_name} \
        --quantile ${quantile} \
        --log_file ${log_file}\
        --log_output "${file_name}_LOG.txt" \
        --output "${file_name}_replicate_zero.csv" \
        --pca_output "${file_name}_replicate_zero_pca"
        
    """
}
