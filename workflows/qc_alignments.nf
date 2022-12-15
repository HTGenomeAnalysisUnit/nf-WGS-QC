nextflow.enable.dsl=2

output_somalier = file("${params.outdir}/somalier")
output_readstats = file("${params.outdir}/readstats")
output_coverage = file("${params.outdir}/coverage")
output_mapflagstats = file("${params.outdir}/mapflag_stats")

include { SAMTOOLS_STATS } from "${projectDir}/modules/samtools" addParams( outdir: "${output_mapflagstats}")
if (params.mode == 'WGS') {
    include { coverage_WGS as COVERAGE } from "${projectDir}/modules/mosdepth" addParams( outdir: "${output_coverage}")
} else if (params.mode == 'WES') {
    include { coverage_WES as COVERAGE } from "${projectDir}/modules/mosdepth" addParams( outdir: "${output_coverage}")
}
include { SOMALIER_EXTRACT } from "${projectDir}/modules/somalier" addParams( outdir: "${output_somalier}")
include { COMPUTE_READSTATS } from "${projectDir}/modules/readstats" addParams( outdir: "${output_readstats}")

workflow RUN_ALIGNEMENTS_QC {
    take:
        /* a tuple containing bam files
        [sample, file.bam, file.bam.bai] */ 
        bamfiles
        regions
        somalier_varset
        genome_ref
        
    main:
        if (params.mode == 'WGS') {
            COVERAGE(bamfiles)
        } else if (params.mode == 'WES') {
            regions_file = file(regions, checkIfExists: true)
            COVERAGE(bamfiles, regions_file)
        }
        SAMTOOLS_STATS(bamfiles)
        SOMALIER_EXTRACT(bamfiles.map{return tuple(it[1], it[2])}, somalier_varset, genome_ref)
        COMPUTE_READSTATS(bamfiles)
}