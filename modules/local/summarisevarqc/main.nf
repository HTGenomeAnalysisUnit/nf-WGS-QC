process SUMMARISE_VARIANT_QC {
    tag "${project_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/numpy_polars_pyyaml:6ee80bfda2a803db'

    input:
    val project_id
    val report_var_qual
    val report_var_dp
    path qc_flags
    path bcftools_stats_region
    path bcftools_idxstats
    path bcftools_stats
    path charr_contamination

    output:
    path ("*-var_qc.tsv"), emit: qc_report_tsv
    path ("*-var_qc.json"), emit: qc_report_json
    path ("*-var_qc.qc_flags.tsv"), emit: qc_flags
    path ("*-var_qc.qc_criteria.tsv"), emit: qc_criteria
    path ("*-var_qc.qc_var_datavzrd.yml"), emit: datavzrd_yml
    tuple val("${task.process}"), val('summarise_var_qc'), eval("summarise_var_qc.py --version"), topic: versions, emit: versions_summarisevarqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${project_id}"
    """
    summarise_var_qc.py \\
        ${args} \\
        --workers ${task.cpus} \\
        --qc_flags ${qc_flags} \\
        --var_qual_report ${report_var_qual} \\
        --var_dp_report ${report_var_dp} \\
        --input_bcftools_idxstats ${bcftools_idxstats} \\
        --input_bcftools_stats ${bcftools_stats} \\
        --input_bcftools_stats_region ${bcftools_stats_region} \\
        --input_charr ${charr_contamination} \\
        --output_prefix ${prefix}-var_qc
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${project_id}"
    """
    echo ${args}
    
    touch ${prefix}-var_qc.tsv
    touch ${prefix}-var_qc.json
    touch ${prefix}-var_qc.qc_flags.tsv
    touch ${prefix}-var_qc.qc_criteria.tsv
    touch ${prefix}-var_qc.qc_var_datavzrd.yml
    """
}
