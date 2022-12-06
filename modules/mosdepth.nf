// Whole genome coverage
/*
INPUT
tuple sample name, bam file, bai index

OUTPUT
folder aligned SAM

PUBLISH
no publish, SAM is usually passed to samtools for BAM conversion

PARAMS
threads: N parallel threads
outdir: publish dir
*/

process coverage_WGS {
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> if (filename.endsWith('.summary.txt')) { "summary/$filename" } else { "$filename" } }
    label 'lowcores'

    input:
    tuple val(sample), file(bam_file), file(bai_file)

    output:
    path "*.{dist.txt,bed.gz}", emit: dist_bed_gz_files
    path "*.summary.txt", emit: summary_files

    script:
    """
    export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
    export MOSDEPTH_Q1=LESS_THAN_5  # 1..4
    export MOSDEPTH_Q2=LOW_COV  # 5..9
    export MOSDEPTH_Q3=CALLABLE  # 10..150
    export MOSDEPTH_Q4=HIGH_COV # 150..

    MOSDEPTH_PRECISION=4 mosdepth -n -t ${task.cpus} --quantize 0:1:5:10:150: ${sample} $bam_file
    """
}

// Region specific coverage
/*
INPUT
tuple sample name, bam file, bai index
bed file: specific regions for coverage

OUTPUT
folder aligned SAM

PUBLISH
no publish, SAM is usually passed to samtools for BAM conversion

PARAMS
threads: N parallel threads
outdir: publish dir
*/
process coverage_WES {
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> if (filename.endsWith('.summary.txt')) { "summmary/$filename" } else { "$filename" } }
    label 'lowcores'

    input:
    tuple val(sample), file(bam_file), file(bai_file)
    file(regions)

    output:
    path "*.{dist.txt,bed.gz}", emit: dist_bed_gz_files
    path "*.summary.txt", emit: summary_files

    script:
    """
    export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
    export MOSDEPTH_Q1=LESS_THAN_5  # 1..4
    export MOSDEPTH_Q2=LOW_COV  # 5..9
    export MOSDEPTH_Q3=CALLABLE  # 10..150
    export MOSDEPTH_Q4=HIGH_COV # 150..

    MOSDEPTH_PRECISION=4 mosdepth -n -t ${task.cpus} --quantize 0:1:5:10:150: --by $regions --thresholds 1,5,10,20 ${sample} $bam_file
    """
}