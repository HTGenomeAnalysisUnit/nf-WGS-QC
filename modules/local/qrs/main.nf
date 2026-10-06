process QUICK_READS_STATS {
    tag "$meta.id"
    label 'process_single'

    // Does not support conda at the moment
    container 'ghcr.io/htgenomeanalysisunit/quick_reads_stats:0.1.1'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path '*.per_sample_run.tsv', emit: per_sample_run
    path '*.rundetails.tsv', emit: rundetails
    path '*.read_stats.tsv', emit: read_stats
    tuple val("${task.process}"), val('quick_reads_stats'), eval("qrs --help | head -1 | sed 's/.\+VERSION //'"), topic: versions, emit: versions_qrs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    qrs \\
        $args \\
        readname_stats -s 5000000 \\
        --out ${prefix} ${bam}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.per_sample_run.tsv
    touch ${prefix}.rundetails.tsv
    touch ${prefix}.read_stats.tsv
    """
}
