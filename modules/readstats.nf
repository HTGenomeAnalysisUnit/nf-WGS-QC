process COMPUTE_READSTATS {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample}.GC_stats/${sample}.GC_stats.tsv", saveAs: { filename -> "GC_stats/${sample}.GC_stats.tsv" }
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sample), file(bam_file), file(bai_file)

    output:
    tuple file("${sample}.per_sample_run.tsv"), file("${sample}.rundetails.tsv"), emit: tsv_stats
    path "${sample}.GC_stats/*", emit: gc_stats

    script:
    """
    readname_stats --out $sample $bam_file
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