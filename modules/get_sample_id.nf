process GET_SAMPLE_ID_BAM {
    input:
    tuple path(bam_file), path(bam_index)

    output:
    tuple stdout, path(bam_file), path(bam_index), emit: bam_files
    tuple val(bam_file.simpleName), stdout, emit: id_map
    //internal_id, external_id

    script:
    """
    samtools view -H $bam_file \
        | grep "@RG" | cut -f6 | sed 's/SM://' | sort -u
    """
}

process GET_SAMPLE_ID_VCF {
    input:
    tuple path(vcf_file), path(vcf_index)

    output:
    tuple stdout, path(vcf_file), path(vcf_index), emit: vcf_files
    tuple val(vcf_file.simpleName), stdout, emit: id_map
    //internal_id, external_id

    script:
    """
    zgrep -m1 "#CHROM" $vcf_file | cut -f10
    """
}
