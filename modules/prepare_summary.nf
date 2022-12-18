process PREPARE_SUMMARY {
    publishDir "${params.summary_output}", mode: 'move'
    
    input:
        file(qc_data)
        file(align_filters)
        file(var_filters)
        file(previous_var_rds)
        file(previous_align_rds)
    
    output:
        file("qc_align_summary.tsv")
        file("qc_vars_summary.tsv")
        file("project_qc_summary.tsv")
        file("qc_align.RDS")
        file("qc_vars.RDS")

    script:
    def align_filters_opt = align_filters.name == 'NO_FILE' ? '' : "--align_filters $align_filters"
    def var_filters_opt = var_filters.name == 'NO_FILE' ? '' : "--var_filters $var_filters"
    def force_update_opt = params.force_update ? "--force_update" : ''
    def update_app_opt = params.update_app ? "--update_app" : ''
    def var_rds_opt = previous_var_rds.exists() ? "--var_qc_rds $previous_var_rds" : ''
    def align_rds_opt = previous_align_rds.exists() ? "--align_qc_rds $previous_align_rds" : ''    
    """
    prepare_project.R \
        --project_qc_dir ./ \
        --project_name ${params.project_id} \
        --app_dir ${params.app_dir} \
        --app_name ${params.app_name} \
        --by_tranche \
        $align_filters_opt \
        $var_filters_opt \
        $var_rds_opt \
        $align_rds_opt \
        $force_update_opt \
        $update_app_opt
    """
}