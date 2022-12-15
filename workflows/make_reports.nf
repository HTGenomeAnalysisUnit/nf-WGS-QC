include { MAKE_PED }  from "${projectDir}/modules/make_ped" addParams( outdir: "${outdir}")
include { SOMALIER_RELATE } from "${projectDir}/modules/somalier" addParams( reportdir: "${output_reports}", outdir: "${output_somalier}")
include { SOMALIER_ANCESTRY } from "${projectDir}/modules/somalier" addParams( reportdir: "${output_reports}", outdir: "${output_somalier}")
include { MULTIQC } from "${projectDir}/modules/multiqc" addParams( outdir: "$output_reports")

output_reports = file("${params.outdir}/reports")

workflow MAKE_REPORTS_WF {
    somalier_labels = file(params.somalier_labels, checkIfExists: true)
    somalier_1kgdata = file(params.somalier_1kgdata, checkIfExists: true)
    somalier_files = Channel.fromPath(params.somalier_profiles)

            fastp_files = params.fastp_json ? Channel.fromPath(params.fastp) : Channel.empty()
        samblaster_files = params.samblaster_log ? Channel.fromPath(params.samblaster) : Channel.empty()

    // Generate a PED file from sample IDs if one is not provided
    if (params.ped_file) {
        ped_file = file(params.ped_file)
    } else {
        project_metadata = file(params.project_metadata, checkIfExists: true)
        MAKE_PED(project_metadata, ext2int_map_ch)
        ped_file = MAKE_PED.out
    }

    SOMALIER_RELATE(somalier_files, ped_file)
    SOMALIER_ANCESTRY(somalier_files, somalier_ancestry_labels, somalier_1kg_data)

    qc_files = fastp_file
            .mix(samblaster_file,
                SAMTOOLS_STATS.out,
                COVERAGE.out.dist_bed_gz_files,
                SOMALIER_RELATE.out,
                SOMALIER_ANCESTRY.out.ancestry_tsv)
            .collect()
    MULTIQC(qc_files, file("${projectDir}/conf/multiqc_config.yaml"))

    //FOr variants
    report_input_ch = BCFTOOLS_STATS.out.join(BCFTOOLS_IDX_STATS.out)
    QC_RMD_SMALLVARS(report_input_ch, somalier_ancestry)

    //Generate tranches data

    //Update app

    //
}