process COMPUTE_CHARR {
    tag "$meta.id"
    label 'process_single'

    container 'ghcr.io/htgenomeanalysisunit/sce-vcf:0.1.4'

    input:
    tuple val(meta), file(vcf), file(vcf_index)
    val(af_field)

    output:
    tuple val(meta), path("*.charr.tsv"), emit: charr_stats
    tuple val("${task.process}"), val('compute_charr'), eval("sceVCF --help 2>1 | head 1 | sed 's/.\\+v//'"), topic: versions, emit: versions_computecharr

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "COMPUTE_CHARR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sceVCF ${args} --out ${prefix}.charr.tsv -a $af_field $vcf
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "COMPUTE_CHARR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.charr.tsv
    """
}
