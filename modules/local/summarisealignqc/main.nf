process SUMMARISE_ALIGN_QC {
    tag "${project_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/numpy_polars_pyyaml:6ee80bfda2a803db'

    input:
    val project_id
    val target_depths
    val target_base_qualities
    path align_qc_flags
    path somalier_ancestry
    path somalier_pairs
    path somalier_samples
    path mosdepth_inputs
    path samtools_idxstats_inputs
    path samtools_stats_inputs
    path verifybamid_inputs
    path qrs_inputs

    output:
    path ("*-align_qc.tsv"), emit: qc_report_tsv
    path ("*-align_qc.json"), emit: qc_report_json
    path ("*-align_qc.qc_flags.tsv"), emit: qc_flags
    path ("*-align_qc.qc_criteria.tsv"), emit: qc_criteria
    path ("*-align_qc.qc_align_datavzrd.yml"), emit: datavzrd_yml
    tuple val("${task.process}"), val('summarise_align_qc'), eval("summarise_align_qc.py --version"), topic: versions, emit: versions_summarisealignqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${project_id}"
    """
    summarise_align_qc.py \\
        ${args} \\
        --workers ${task.cpus} \\
        --target_depths ${target_depths} \\
        --target_base_qualities ${target_base_qualities} \\
        --somalier_ancestry ${somalier_ancestry} \\
        --somalier_pairs ${somalier_pairs} \\
        --somalier_samples ${somalier_samples} \\
        --input_mosdepth ${mosdepth_inputs} \\
        --input_samtools_idxstats ${samtools_idxstats_inputs} \\
        --input_samtools_stats ${samtools_stats_inputs} \\
        --input_verifybamid ${verifybamid_inputs} \\
        --input_qrs ${qrs_inputs} \\
        --qc_flags ${align_qc_flags} \\
        --output_prefix ${prefix}-align_qc
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${project_id}"
    """
    echo ${args}
    
    touch ${prefix}-align_qc.tsv
    touch ${prefix}-align_qc.json
    touch ${prefix}-align_qc.qc_flags.tsv
    touch ${prefix}-align_qc.qc_criteria.tsv
    touch ${prefix}-align_qc.qc_align_datavzrd.yml
    """
}
