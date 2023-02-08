process COUNT_BASES {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), file(bam_file), file(bai_file)

    output:
    path "${sample}.totbases"

    script:
    """
    tot_bases=\$(samtools view $bam_file | awk '{bases += length(\$10)}; END {print bases}')
    genome_size=\$(samtools view -H $bam_file | grep "@SQ" | grep -v "HLA\\|decoy\\|EBV" | cut -f3 | sed 's/LN://g' | awk '{ bases += \$1 }; END { print bases }')

    echo -e "sample\tN_bases\tGenome_size" > ${sample}.totbases
    echo -e "${sample}\t\${tot_bases}\t\${genome_size}" >> ${sample}.totbases
    """
}

/*
process CONCAT_READSTATS {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(read_stats_tsv)

    output:
    tuple file('per_sample_readstats.tsv'), file('per_sample_seqruns.tsv')

    script:
    """
    headfile=\$(ls *.per_sample_run.tsv | head -1)
    head -1 \$headfile > per_sample_readstats.tsv 
    for f in *.per_sample_run.tsv
    do
        tail -n+2 \$f >> per_sample_readstats.tsv
    done

    headfile=\$(ls *.rundetails.tsv | head -1)
    head -1 \$headfile > per_sample_seqruns.tsv 
    for f in *.rundetails.tsv
    do
        tail -n+2 \$f >> per_sample_seqruns.tsv
    done
    """
}
*/