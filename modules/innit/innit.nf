process innit {

    container "tobiasbueschertm/nf-innit"
   // conda = "nextflow_env"
    input:
    path input_data
    path script
    val file_name


    output:

    path "${file_name}_LOG.txt" , emit: LOG
    path "${file_name}_innit_pca/*", emit: pca

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate nextflow_env
  
    python3 ${script} \
        --input_data ${input_data} \
        --file_name ${file_name} \
        --log_output "${file_name}_LOG.txt"\
        --pca-output "${file_name}_innit_pca"


    """






}
