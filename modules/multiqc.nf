process MULTIQC {
    label 'singlethread'
    publishDir "$params.outdir", mode: 'copy'

    input:
    file(qc_files)
    file(multiqc_config)

    output:
    path 'multiqc_report.html', emit: html_report
    path 'multiqc_data', emit: data_dir

    script:
    """
    multiqc .
    """
}