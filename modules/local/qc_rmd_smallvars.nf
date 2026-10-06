process QC_RMD_SMALLVARS {
    publishDir "${params.reportdir}/${tranche_name}", mode: 'move'
    
    input:
        tuple val(tranche_name), file(bcftools_stats), file(bcftools_idx_stats), file(somalier_ancestry)
    
    output:
        file("${tranche_name}_QC_smallvars.html")

    script:
    """
    cp /opt/app/*.R* \$PWD
    R -e "rmarkdown::render('QC_report_smallvars.Rmd', \
        output_file='${tranche_name}_QC_smallvars.html', \
        params=list( \
            doc_title='${tranche_name} smallvars QC', \
            vcf_stats='$bcftools_stats', \
            chr_counts='$bcftools_idx_stats', \
            ancestry='$somalier_ancestry' \
        ))"
    """
}