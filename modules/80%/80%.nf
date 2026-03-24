process missing_percent {

    container  "tobiasbueschertm/nf-missing_percent"

    input:
    val input_data
    val log_file
    path script
    val file_name
   

    
    output:
    path "${file_name}_LOG.txt" , emit: LOG
    path "${file_name}_80%.csv", emit: df
    path "${file_name}_80%_pca/*", emit: pca

    script:
    """
    
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate nextflow_env
    python3 ${script} \
        --input_data ${input_data} \
        --file_name ${file_name} \
        --log_file ${log_file}\
        --log_output "${file_name}_LOG.txt" \
        --output "${file_name}_80%.csv" \
        --pca_output "${file_name}_80%_pca"
        
    """
}
