process GET_SAMPLE_ID_BAM {
    label 'small_on_local'
    
    input:
    tuple path(bam_file), path(bam_index)

    output:
    tuple stdout, path(bam_file), path(bam_index), emit: bam_files
    tuple val(bam_file.simpleName), stdout, emit: id_map
    //internal_id, external_id

    script:
    """
    samtools view -H $bam_file \
        | grep "@RG" | tr "\\t" "\\n" | grep "SM:" | sed 's/SM://' | sort -u | tr -d '\\n'
    """
}

process GET_SAMPLE_ID_VCF {
    label 'small_on_local'

    input:
    tuple path(vcf_file), path(vcf_index)

    output:
    tuple stdout, path(vcf_file), path(vcf_index), emit: vcf_files
    tuple val(vcf_file.simpleName), stdout, emit: id_map
    //internal_id, external_id

    script:
    """
    zgrep -m1 "#CHROM" $vcf_file | cut -f10 | tr -d '\\n'
    """
}
