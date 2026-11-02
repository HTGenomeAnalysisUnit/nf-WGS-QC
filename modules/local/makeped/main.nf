process MAKE_PED {
    tag "${project_id}_ped"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/pandas:3.1.0--f6c858a02e385828'

    input:
    path manifest
    val project_id

    output:
    path "*.ped", emit: ped_file
    tuple val("${task.process}"), val('make_ped'), eval("make_ped.py --version | cut -d' ' -f2"), topic: versions, emit: versions_makeped

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    make_ped.py \\
        $args \\
        --out ${project_id}.ped \\
        $manifest
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    echo $args
    
    touch ${project_id}.ped
    """
}
