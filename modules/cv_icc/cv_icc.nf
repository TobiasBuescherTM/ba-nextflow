process cv_icc {

    container  "tobiasbueschertm/nf-cv_icc"

    input:
    val input_data
    val log_file
    path script
    val file_name
    val cv_cutoff
    val icc_cutoff
    

    
    output:
    path "${file_name}_LOG.txt" , emit: LOG
    path "${file_name}_cv_icc.csv", emit: df
    path "${file_name}_cv_icc_pca/*", emit: pca
    path "${file_name}_cv_correlations/*", emit: cv_correlations

    script:
    """
    
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate nextflow_env
    python3 ${script} \
        --input_data ${input_data} \
        --file_name ${file_name} \
        --cv_cutoff ${cv_cutoff}\
        --icc_cutoff ${icc_cutoff}\
        --log_file ${log_file}\
        --log_output "${file_name}_LOG.txt" \
        --output "${file_name}_cv_icc.csv" \
        --pca_output "${file_name}_cv_icc_pca" \
        --correlations_path "${file_name}_cv_correlations"
        
    """
}
