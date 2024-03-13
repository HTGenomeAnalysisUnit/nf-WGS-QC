process PREPARE_SUMMARY {
    label 'lowcores'
    publishDir "${params.summary_output}", mode: 'copy'
    
    input:
        file(qc_data)
        file(align_fails)
        file(var_fails)
        file(align_warnings)
        file(var_warnings)
        file(previous_var_rds)
        file(previous_align_rds)
    
    output:
        file("qc_align_summary.tsv")
        file("qc_vars_summary.tsv")
        file("project_qc_summary.tsv")
        file("qc_align.RDS")
        file("qc_vars.RDS")
        file("prepare_summary.log")

    script:
    def align_fails_opt = align_fails.name == 'NO_FILE' ? '' : "--align_fails $align_fails"
    def var_fails_opt = var_fails.name == 'NO_FILE' ? '' : "--var_fails $var_fails"
    def align_warnings_opt = align_warnings.name == 'NO_FILE' ? '' : "--align_warnings $align_warnings"
    def var_warnings_opt = var_warnings.name == 'NO_FILE' ? '' : "--var_warnings $var_warnings"
    def force_update_opt = params.force_update ? "--force_update" : ''
    def update_app_opt = params.update_app ? "--update_app" : ''
    def by_tranche_opt = params.by_tranche ? "--by_tranche" : ''
    def var_rds_opt = previous_var_rds.exists() ? "--var_qc_rds $previous_var_rds" : ''
    def align_rds_opt = previous_align_rds.exists() ? "--align_qc_rds $previous_align_rds" : ''    
    def proj_qc_dir = params.by_tranche ? './' : qc_data[0].baseName
    """
    prepare_project.R \
        --project_qc_dir ${proj_qc_dir} \
        --project_name ${params.project_id} \
        --app_dir ${params.app_dir} \
        --app_name ${params.app_name} \
        $by_tranche_opt \
        $align_fails_opt \
        $var_fails_opt \
        $align_warnings_opt \
        $var_warnings_opt \
        $var_rds_opt \
        $align_rds_opt \
        $force_update_opt \
        $update_app_opt
    
    cp .command.log prepare_summary.log
    """
}