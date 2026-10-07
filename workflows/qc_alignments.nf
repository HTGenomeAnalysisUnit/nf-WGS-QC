/*
output_somalier = file("${params.outdir}/somalier")
output_readstats = file("${params.outdir}/readstats")
output_coverage = file("${params.outdir}/coverage")
output_mapflagstats = file("${params.outdir}/mapflag_stats")
*/

include { SAMTOOLS_STATS } from "../modules/nf-core/samtools/stats"
include { SAMTOOLS_IDXSTATS } from "../modules/nf-core/samtools/idxstats"
include { MOSDEPTH } from "../modules/nf-core/mosdepth" 
include { SOMALIER_EXTRACT } from "../modules/nf-core/somalier/extract" 
include { QUICK_READS_STATS } from "../modules/local/qrs" 

workflow RUN_ALIGNEMENTS_QC {
    take:
        bamfiles // [meta, file(bam), file(bai)]
        regions_file // file(bed)
        genome_fasta // file(genome.fasta)
        genome_fai // file(genome.fasta.fai)
        somalier_varset // file(somalier_variants.vcf.gz)
        
    main:
        mosdepth_input = bamfiles.combine(regions_file)
        MOSDEPTH(mosdepth_input, ['genome_fasta', genome_fasta], params.mosdepth_quantize_bins)
        SAMTOOLS_STATS(bamfiles, ['ref_genome', genome_fasta, genome_fai])
        SAMTOOLS_IDXSTATS(bamfiles)
        SOMALIER_EXTRACT(bamfiles, ['genome_fasta', genome_fasta], ['genome_fai', genome_fai], ['somalier_varset', somalier_varset])
        QUICK_READS_STATS(bamfiles)

    emit:
        output_coverage = MOSDEPTH.out
        output_mapflagstats = SAMTOOLS_STATS.out
        output_somalier = SOMALIER_EXTRACT.out
        output_readstats = QUICK_READS_STATS.out

    publish:
        output_coverage
        output_mapflagstats
        output_somalier
        output_readstats
}