nextflow.enable.dsl=2

output_somalier = file("${params.outdir}/somalier")
output_readstats = file("${params.outdir}/readstats")
output_fastp = file("${params.outdir}/fastp")
output_coverage = file("${params.outdir}/coverage")
output_dedup = file("${params.outdir}/dedup")
output_mapflagstats = file("${params.outdir}/mapflag_stats")
output_reports = file("${params.outdir}/reports")

include { SAMTOOLS_STATS } from "${projectDir}/modules/samtools" addParams( outdir: "${output_mapflagstats}")
if (params.mode == 'WGS') {
    include { coverage_WGS as COVERAGE } from "${projectDir}/modules/mosdepth" addParams( outdir: "${output_coverage}")
} else if (params.mode == 'WES') {
    include { coverage_WES as COVERAGE } from "${projectDir}/modules/mosdepth" addParams( outdir: "${output_coverage}")
}
include { SOMALIER_EXTRACT } from "${projectDir}/modules/somalier" addParams( outdir: "${output_somalier}")
include { SOMALIER_RELATE } from "${projectDir}/modules/somalier" addParams( reportdir: "${output_reports}", outdir: "${output_somalier}")
include { SOMALIER_ANCESTRY } from "${projectDir}/modules/somalier" addParams( reportdir: "${output_reports}", outdir: "${output_somalier}")
include { READSTATS } from "${projectDir}/modules/readstats" addParams( outdir: "${output_readstats}")
include { MULTIQC } from "${projectDir}/modules/multiqc" addParams( outdir: "$output_reports")

workflow RUN_ALIGNEMENTS_QC {
    take:
        /* a tuple containing bam files
        [sample, file.bam, file.bam.bai] */ 
        bamfiles
        fastp_file
        samblaster_file
        regions
        somalier_varset
        somalier_ancestry_labels
        somalier_1kg_data
        ped_file
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
        SOMALIER_RELATE(SOMALIER_EXTRACT.out.collect(), ped_file)
        SOMALIER_ANCESTRY(SOMALIER_EXTRACT.out.collect(), somalier_ancestry_labels, somalier_1kg_data)
        READSTATS(bamfiles)

        qc_files = fastp_file
            .mix(samblaster_file,
                SAMTOOLS_STATS.out,
                COVERAGE.out.dist_bed_gz_files,
                SOMALIER_RELATE.out,
                SOMALIER_ANCESTRY.out.ancestry_tsv)
            .collect()
        MULTIQC(qc_files, file("${projectDir}/conf/multiqc_config.yaml"))

    emit:
        multiqc_report = MULTIQC.out.html_report
        somalier_ancestry = SOMALIER_ANCESTRY.out.ancestry_tsv
}