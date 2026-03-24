nextflow.enable.dsl=2

// Include modules
include {innit} from "./modules/innit/innit.nf"
include {replicate_zero} from "./modules/replicate_zero/replicate_zero.nf"
include {qcomics} from "./modules/qcomics/qcomics.nf"
include {missing_percent} from "./modules/80%/80%.nf"
include {cv_icc} from "./modules/cv_icc/cv_icc.nf"
include {imputation} from "./modules/imputation/imputation.nf"
include {vsn} from "./modules/vsn/vsn.nf"
 

params {
    //init
    input_data: Path 
    file_name: String 
    innit_script: Path 

    //replicate_zero
    replicate_zero_script: Path
    quantile: Boolean 

    //qcomics
    qcomics_script: Path 
    
    //missing_percent
    missing_percent_script: Path 

    //cv_icc
    cv_icc_script: Path 
    cv_cutoff: Float 
    icc_cutoff: Float 

    //imputation
    imputation_script: Path 
    imputation_start_script: Path 
    imputation_end_script: Path 

    //vsn
    vsn_script: Path 
    vsn_start_script: Path 
    vsn_end_script: Path 
   
}

workflow {

    main:

    innit(
        params.input_data
        ,params.innit_script
        ,params.file_name

    )


    replicate_zero(
        params.input_data
        ,innit.out.LOG
        ,params.replicate_zero_script
        ,params.file_name
        ,params.quantile
    )

    qcomics(
        replicate_zero.out.df
        ,replicate_zero.out.LOG
        ,params.qcomics_script
        ,params.file_name
    )

    missing_percent(
        qcomics.out.df
        ,qcomics.out.LOG
        ,params.missing_percent_script
        ,params.file_name
        
    )

    cv_icc(
        missing_percent.out.df
        ,missing_percent.out.LOG
        ,params.cv_icc_script
        ,params.file_name
        ,params.cv_cutoff
        ,params.icc_cutoff
    )

    imputation(
        cv_icc.out.df
        ,cv_icc.out.LOG
        ,params.imputation_script
        ,params.file_name
        ,params.imputation_start_script
        ,params.imputation_end_script
    )

    vsn(
       imputation.out.df
        ,imputation.out.LOG
        ,params.vsn_script
        ,params.file_name
        ,params.vsn_start_script
        ,params.vsn_end_script
    )

    publish:
    //innit
    log_innit_out = innit.out.LOG
    innit_pca = innit.out.pca
    //replicate_zero
    log_replicate_zero_out = replicate_zero.out.LOG
    replicate_zero_df = replicate_zero.out.df
    replicate_zero_pca = replicate_zero.out.pca
    //qcomics
    log_qcomics_out = qcomics.out.LOG
    qcomics_df = qcomics.out.df
    qcomics_pca = qcomics.out.pca
    //missing_percent
    log_missing_percent_out = missing_percent.out.LOG
    missing_percent_df = missing_percent.out.df
    missing_percent_pca = missing_percent.out.pca
    //cv_icc
    log_cv_icc_out = cv_icc.out.LOG
    cv_icc_df = cv_icc.out.df
    cv_icc_pca = cv_icc.out.pca
    cv_icc_correlations = cv_icc.out.cv_correlations
    //imputation
    log_imputation_out = imputation.out.LOG
    imputation_df = imputation.out.df
    imputation_pca = imputation.out.pca
    //vsn
    log_vsn_out = vsn.out.LOG
    vsn_df = vsn.out.df
    vsn_pca = vsn.out.pca
    final_df = vsn.out.df

}


output {

      //innit
    log_innit_out {
        path "outputs/${params.file_name}"
        mode "copy"

    }
    innit_pca {
        path "outputs/${params.file_name}/innit"
        mode "copy"
    }

    //replicate_zero
    log_replicate_zero_out {
        path "outputs/${params.file_name}"
        mode "copy"
    }
    replicate_zero_df {
        path "outputs/${params.file_name}/replicate_zero"
        mode "copy"
    }
    replicate_zero_pca {
        path "outputs/${params.file_name}/replicate_zero"
        mode "copy"
    }

    //qcomics
    log_qcomics_out {
        path "outputs/${params.file_name}"
        mode "copy"
    }
    qcomics_df {
        path "outputs/${params.file_name}/qcomics"
        mode "copy"
    }
    qcomics_pca {
        path "outputs/${params.file_name}/qcomics"
        mode "copy"
    }

    //missing_percent
    log_missing_percent_out {
        path "outputs/${params.file_name}"
        mode "copy"
    }
    missing_percent_df {
        path "outputs/${params.file_name}/missing_percent"
        mode "copy"
    }
    missing_percent_pca {
        path "outputs/${params.file_name}/missing_percent"
        mode "copy"
    }

    //missing_percent
    log_cv_icc_out {
        path "outputs/${params.file_name}"
        mode "copy"
    }
    cv_icc_df {
        path "outputs/${params.file_name}/cv_icc"
        mode "copy"
    }
    cv_icc_pca {
        path "outputs/${params.file_name}/cv_icc"
        mode "copy"
    }
    cv_icc_correlations {
        path "outputs/${params.file_name}/cv_icc"
        mode "copy"
    }

    //imputation
    log_imputation_out {
        path "outputs/${params.file_name}"
        mode "copy"
    }
    imputation_df {
        path "outputs/${params.file_name}/imputation"
        mode "copy"
    }
    imputation_pca {
        path "outputs/${params.file_name}/imputation"
        mode "copy"
    }

    //vsn
    log_vsn_out {
        path "outputs/${params.file_name}"
        mode "copy"
    }
    vsn_df {
        path "outputs/${params.file_name}/vsn"
        mode "copy"
    }
    vsn_pca {
        path "outputs/${params.file_name}/vsn"
        mode "copy"
    }

   final_df {
        path "outputs/${params.file_name}"
        mode "copy"
    }
    
}
