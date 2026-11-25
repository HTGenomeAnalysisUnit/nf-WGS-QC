process PREPARE_QC_REPORT {
    tag "${project_id}"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/numpy_polars_pyyaml:6ee80bfda2a803db'

    input:
    val project_id
    path align_datavzrd_yml
    path align_metrics_tsv
    path align_qc_flags_tsv
    path align_qc_criteria_tsv
    path var_datavzrd_yml
    path var_metrics_tsv
    path var_qc_flags_tsv
    path var_qc_criteria_tsv

    output:
    path ("*.datavzrd.yml"), emit: datavzrd_yml
    path ("*.project_qc_summary.js"), emit: report_js
    path ("*.project_qc_summary.tsv"), emit: qc_summary_tsv
    path ("*.qc_flags_summary.tsv"), emit: flags_summary_tsv
    path ("*.qc_flags.tsv"), emit: flags_tsv
    tuple val("${task.process}"), val('prepar_qc_report'), eval("create_qc_report.py --version"), topic: versions, emit: versions_prepareqcreport

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${project_id}"
    """
    create_qc_report.py \\
        ${args} \\
        --project_id ${project_id} \\
        --align_datavzrd_yml ${align_datavzrd_yml} \\
        --align_metrics_tsv ${align_metrics_tsv} \\
        --align_qc_flags_tsv ${align_qc_flags_tsv} \\
        --align_qc_criteria_tsv ${align_qc_criteria_tsv} \\
        --var_datavzrd_yml ${var_datavzrd_yml} \\
        --var_metrics_tsv ${var_metrics_tsv} \\
        --var_qc_flags_tsv ${var_qc_flags_tsv} \\
        --var_qc_criteria_tsv ${var_qc_criteria_tsv} \\
        --output_prefix ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${project_id}"
    """
    echo ${args}
    
    touch ${prefix}.datavzrd.yml
    touch ${prefix}.project_qc_summary.tsv
    touch ${prefix}.qc_flags_summary.tsv
    touch ${prefix}.project_qc_summary.js
    touch ${prefix}.qc_flags.tsv
    """
}
