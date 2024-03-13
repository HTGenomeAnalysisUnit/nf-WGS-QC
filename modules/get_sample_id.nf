process GET_SAMPLE_ID_BAM {
    label 'small_on_local'
    
    input:
    tuple path(bam_file), path(bam_index)

    output:
    tuple stdout, path(bam_file), path(bam_index), emit: bam_files
    //internal_id, external_id

    script:
    """
    samtools view -H $bam_file \
        | grep "^@RG" | tr "\\t" "\\n" | grep "SM:" | sed 's/SM://' | sort -u | tr -d '\\n'
    """
}

process GET_SAMPLE_ID_VCF {
    label 'small_on_local'

    input:
    tuple path(vcf_file), path(vcf_index)

    output:
    tuple stdout, path(vcf_file), path(vcf_index), emit: vcf_files
    //internal_id, external_id

    script:
    """
    zgrep -m1 "#CHROM" $vcf_file | cut -f10 | tr -d '\\n'
    """
}

process GET_BAM_MANIFEST {
    label 'small_on_local'
    
    input:
    tuple path(bam_file), path(bam_index)

    output:
    tuple val(bam_file.simpleName), stdout
    //internal_id, external_id

    script:
    """
    sampleid=\$(samtools view -H $bam_file | grep "^@RG" | tr "\\t" "\\n" | grep "SM:" | sed 's/SM://' | sort -u)
    file_date=\$(date -r $bam_file "+%Y%m%d")
    file_path=\$(readlink -f $bam_file)
    echo -e "\$sampleid\\t\$file_date\\t\$file_path" | tr -d '\\n'
    """
}

process GET_VCF_MANIFEST {
    label 'small_on_local'

    input:
    tuple path(vcf_file), path(vcf_index)

    output:
    tuple val(vcf_file.simpleName), stdout
    //internal_id, external_id

    script:
    """
    sampleid=\$(zgrep -m1 "#CHROM" $vcf_file | cut -f10 | tr -d '\\n')
    file_date=\$(date -r $vcf_file "+%Y%m%d")
    file_path=\$(readlink -f $vcf_file)
    echo -e "\$sampleid\\t\$file_date\\t\$file_path" | tr -d '\\n'
    """
}
