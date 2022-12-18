// Expects tarnches subfolders in the main input folder

output_reports = file("${params.summary_output}/per_tranche_reports")

include { MAKE_PED }  from "${projectDir}/modules/make_ped"
include { SOMALIER_RELATE } from "${projectDir}/modules/somalier" addParams( reportdir: "${output_reports}")
include { SOMALIER_ANCESTRY } from "${projectDir}/modules/somalier" addParams( reportdir: "${output_reports}")
include { MULTIQC } from "${projectDir}/modules/multiqc" addParams( reportdir: "${output_reports}")
include { QC_RMD_SMALLVARS } from "${projectDir}/modules/qc_rmd_smallvars" addParams( reportdir: "${output_reports}")
include { PREPARE_SUMMARY } from "${projectDir}/modules/prepare_summary"

workflow MAKE_REPORTS_WF {
    somalier_labels = file(params.somalier_labels, checkIfExists: true)
    somalier_1kg_data = file(params.somalier_1kgdata, checkIfExists: true)
    manifest_file = file(params.project_manifest, checkIfExists: true)
    align_filters_file = file(params.align_filters)
    var_filters_file = file(params.var_filters)
    var_qc_rds = file(params.previous_var_qc_rds)
    align_qc_rds = file(params.previous_align_qc_rds)

    // Generate a PED file from sample IDs if one is not provided
    if (params.ped_file) {
        ped_file = file(params.ped_file)
    } else {
        //read the ext2int_map file
        project_metadata = file(params.project_metadata, checkIfExists: true)
        MAKE_PED(project_metadata, manifest_file)
        ped_file = MAKE_PED.out
    }

    //Make channels per tranche for each data type
    tranches_ch = Channel
        .fromPath("${params.main_qc_dir}/*", type: 'dir')
        .multiMap { 
            somalier: tuple(it.baseName, file("$it/qc_alignment/somalier/*"))
            coverage: tuple(it.baseName, file("$it/qc_alignment/coverage/**"))
            sam_stats: tuple(it.baseName, file("$it/qc_alignment/mapflag_stats/*"))
            var_stats: tuple(it.baseName, file("$it/qc_variants/variant_stats/*.vars.stats"), file("$it/qc_variants/variant_stats/*.varsPerContig.tsv"))
            full_tranche: tuple(it.baseName, it)
        }

    //Run Somalier per tranche
    SOMALIER_RELATE(tranches_ch.somalier, ped_file)
    SOMALIER_ANCESTRY(tranches_ch.somalier, somalier_labels, somalier_1kg_data)

    //Report for alignment QC
    qc_files = tranches_ch.coverage
            .join(tranches_ch.sam_stats)
            .join(SOMALIER_RELATE.out)
            .join(SOMALIER_ANCESTRY.out.ancestry_tsv)
            //.groupTuple()
            //.view()
            
    MULTIQC(qc_files, file("${projectDir}/conf/multiqc_config.yaml"))

    //Generate summarised data per tranche
    if (params.make_summary) {
    full_tranche_data = tranches_ch.full_tranche
        .join(SOMALIER_RELATE.out)
        .join(SOMALIER_ANCESTRY.out.ancestry_tsv)
        .map { tuple(it[1],it[2],it[3]) }
        .flatten()
        .collect()
        
    //Update summary and app
    PREPARE_SUMMARY(full_tranche_data, file(params.align_filters, checkIfExists: true), file(params.var_filters, checkIfExists: true), var_qc_rds, align_qc_rds)
    }
}