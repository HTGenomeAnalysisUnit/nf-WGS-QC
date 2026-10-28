process QUICK_READS_STATS {
    tag "$meta.id"
    label 'process_single'
    label 'process_mid_memory'

    // Does not support conda at the moment
    container 'ghcr.io/htgenomeanalysisunit/quick_reads_stats:0.1.2'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.per_sample_run.tsv'), emit: per_sample_run
    tuple val(meta), path('*.rundetails.tsv'), emit: rundetails
    tuple val(meta), path('*.read_stats.tsv'), emit: read_stats
    tuple val("${task.process}"), val('quick_reads_stats'), eval("qrs --help | head -1 | sed 's/.\\+VERSION //'"), topic: versions, emit: versions_qrs

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QUICK_READS_STATS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    qrs \\
        $args \\
        --out ${prefix} ${bam}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QUICK_READS_STATS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.per_sample_run.tsv
    touch ${prefix}.rundetails.tsv
    touch ${prefix}.read_stats.tsv
    """
}
