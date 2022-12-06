// samtools operations

// MERGE
// merge multiple BAMs
/*
INPUT
sample name, list of bam files
A tuple of files to be merged per sample
like: (sample_id, [part1.bam, part2.bam])

OUTPUT
sample name, merged.bam

PUBLISH
resulting BAM not published since it's usually passed to dedup
*/

process SAMTOOLS_MERGE {
    label 'lowcores'

    input:
    tuple val(sample), file(myfiles)

    output:
    tuple val(sample), file('merged.bam'), emit: merged_bam

    script:
    def input_args = myfiles.collect{ "$it" }.join(" ")
    """
    samtools merge -n -f -@ ${task.cpus} merged.bam $input_args
    """
}

// SORT
/*
This will dynamically allocate memory: 4Gb * task.cpus

INPUT
sample name, bam file

OUTPUT
sample name, sorted.bam, sorted.bam.bai

PUBLISH
${sample}_sorted.bam and .bai

PARAMS
outdir: publish dir
*/
process SAMTOOLS_SORT {
    memory { 4.GB * task.cpus }
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.bam*"

    input:
    tuple val(sample), file(bam_file)
    
    output:
    tuple val(sample), file("${sample}.${bam_file.getName().replaceFirst(/\.bam$/, ".sorted.bam")}"), file("${sample}.${bam_file.getName().replaceFirst(/\.bam$/, ".sorted.bam.bai")}"), emit: sorted_bam
    tuple val(sample), val("${file(params.outdir)}/${sample}.${bam_file.getName().replaceFirst(/\.bam$/, ".sorted.bam")}"), emit: bam_path

    script:
    def prefix = bam_file.getName().replaceFirst(/\.bam$/, "")
    """
    samtools sort -@ ${task.cpus} -m 4G -T bam_tmp -o ${sample}.${prefix}.sorted.bam $bam_file
    samtools index -@ ${task.cpus} ${sample}.${prefix}.sorted.bam
    """
}

// STATS
// Run both flagstat and idxstat
/*
INPUT
sample name, bam file, bai file

OUTPUT
sample.flagstat, sample.mapstat

PUBLISH
${sample}.flagstat, ${sample}.mapstats

PARAMS
threads: N parallel threads
outdir: publish dir
*/
process SAMTOOLS_STATS {
    label 'lowcores'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), file(bam_file), file(bai_file)

    output:
    file("*.{flagstat,idxstat}")
    
    script:
    """
    samtools flagstat -@ ${task.cpus} $bam_file > ${sample}.flagstat
    samtools idxstats $bam_file > ${sample}.idxstat
    """
}

process SAMTOOLS_SAM2BAM {
    label 'lowcores'
    if (params.publishbam) {
        publishDir "${params.outdir}", mode: 'copy'
    }

    input:
    tuple val(sample), file(sam_file)

    output:
    tuple val(sample), file("${sam_file.getName().replaceFirst(/\.sam$/, ".bam")}")
    
    script:
    def prefix = sam_file.getName().replaceFirst(/\.sam$/, "")
    """
    samtools view -@ ${task.cpus} -b $sam_file > ${prefix}.bam
    """
}