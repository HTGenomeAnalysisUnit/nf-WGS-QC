process MULTIQC {
    label 'singlethread'
    publishDir "${params.reportdir}/${tranche_name}", mode: 'copy'

    input:
    tuple val(tranche_name), file(coverage_files), file(samtools_stats), file(somalier_relate), file(somalier_ancestry)
    file(multiqc_config)

    output:
    path 'multiqc_report.html', emit: html_report
    path 'multiqc_data', emit: data_dir

    script:
    """
    multiqc .
    """
}