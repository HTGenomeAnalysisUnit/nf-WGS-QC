process MAKE_PED {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(project_metadata)
    path(ext2int_ids_map)

    output:
    path('ped_file.ped')

    script:
    """
    ./make_ped.py $project_metadata $ext2int_ids_map
    """
}