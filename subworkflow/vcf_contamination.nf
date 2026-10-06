workflow VCF_CONTAMINATION {
    take:
        vcf_files
        annotation_zip

    main:
        BCFTOOLS_FILTER_NORM(vcf_files)
        ANNOTATE_AF(BCFTOOLS_FILTER_NORM.out, annotation_zip)
        COMPUTE_CHARR(ANNOTATE_AF.out)
        //CONCAT_STATS(COMPUTE_CHARR.out.tsv_stats.collect())

    emit:
        tsv_stats = COMPUTE_CHARR.out
}

process BCFTOOLS_FILTER_NORM {
    label 'singlethread'

    input:
    tuple val(sample), file(vcf_file), file(vcf_index)

    output:
    tuple val(sample), file("${sample}.norm.vcf.gz"), file("${sample}.norm.vcf.gz.csi")

    script:
    """
    bcftools view -f "PASS,." -Ou $vcf_file \
        | bcftools norm -m- -Oz > ${sample}.norm.vcf.gz
    bcftools index ${sample}.norm.vcf.gz
    """ 
}   

process ANNOTATE_AF {
    label 'singlethread'

    input:
    tuple val(sample), file(vcf_file), file(vcf_index)
    file(annotation_zip)

    output:
    tuple val(sample), file("${sample}.anno.vcf.gz"), file("${sample}.anno.vcf.gz.csi")

    script:
    """
    echtvar anno -e $annotation_zip $vcf_file ${sample}.anno.vcf.gz
    tabix -p vcf --csi ${sample}.anno.vcf.gz
    """ 
}

process COMPUTE_CHARR {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), file(vcf_file), file(vcf_index)

    output:
    path("${sample}.cont.tsv"), emit: tsv_stats

    script:
    """
    sceVCF --out ${sample}.cont.tsv -a gnomad_af $vcf_file
    """
}

/*
process CONCAT_STATS {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(charr_stats_tsv)

    output:
    path('contamination_stats.tsv')

    script:
    """
    headfile=\$(ls *.cont.tsv | head -1)
    head -1 \$headfile > contamination_stats.tsv 
    for f in *.cont.tsv
    do
        tail -n+2 \$f >> contamination_stats.tsv
    done
    """
}
*/