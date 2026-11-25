include { samplesheetToList } from 'plugin/nf-schema'
include { MAKE_PED } from '../modules/local/makeped'
include { SOMALIER_RELATE } from '../modules/local/somalier/relate'
include { SOMALIER_ANCESTRY } from '../modules/local/somalier/ancestry'
include { SUMMARISE_ALIGN_QC } from '../modules/local/summarisealignqc'
include { SUMMARISE_VARIANT_QC } from '../modules/local/summarisevarqc'
include { PREPARE_QC_REPORT } from '../modules/local/prepareqcreport'
include { DATAVZRD } from '../modules/nf-core/datavzrd'

workflow MAKE_REPORTS {

    main:
    somalier_labels = file(params.somalier_labels, checkIfExists: true)
    somalier_1kg_data = file(params.somalier_1kgdata, checkIfExists: true)
    align_qc_flags = file(params.alignment_qc_criteria, checkIfExists: true)
    var_qc_flags = file(params.variant_qc_criteria, checkIfExists: true)

    // Collect manifest files across input folders
    // Alignment QC manifests
    align_qc_folders = params.project_qc_folder ? channel.fromPath("${params.project_qc_folder}/qc_alignments", type: 'dir', checkIfExists: true) : channel.empty()

    if (params.alignment_qc_inputs) {
        align_qc_folders = align_qc_folders.mix(channel.fromList(samplesheetToList(params.alignment_qc_inputs, "${projectDir}/assets/schema_input_reports.json")))
    }

    align_manifests_ch = align_qc_folders.multiMap { qc_alignment_dir ->
        mosdepth: file("${qc_alignment_dir}/published_files_manifest/coverage_summary.csv", checkIfExists: true)
        qrs: file("${qc_alignment_dir}/published_files_manifest/qrs_*.csv", checkIfExists: true)
        samtools_idxstats: file("${qc_alignment_dir}/published_files_manifest/samtools_idxstats.csv", checkIfExists: true)
        samtools_stats: file("${qc_alignment_dir}/published_files_manifest/samtools_stats.csv", checkIfExists: true)
        somalier: file("${qc_alignment_dir}/published_files_manifest/somalier_extract.csv", checkIfExists: true)
        verifybamid: file("${qc_alignment_dir}/published_files_manifest/verifybamid_result.csv", checkIfExists: true)
        bam_manifest: file("${qc_alignment_dir}/input_bam_manifest.tsv", checkIfExists: true)
    }

    align_collected_manifests_ch = align_manifests_ch
        .getNames()
        .collectEntries { name ->
            def collected_ch = align_manifests_ch
                .getProperty(name)
                .collectFile(
                    name: "${name}_inputs.tsv",
                    keepHeader: true,
                )

            [(name): collected_ch]
        }

    // Variants QC manifests
    var_qc_folders = params.project_qc_folder ? channel.fromPath("${params.project_qc_folder}/qc_variants", type: 'dir', checkIfExists: true) : channel.empty()

    if (params.variant_qc_inputs) {
        var_qc_folders = var_qc_folders.mix(channel.fromList(samplesheetToList(params.variant_qc_inputs, "${projectDir}/assets/schema_input_reports.json")))
    }

    var_manifests_ch = var_qc_folders.multiMap { qc_variant_dir ->
        def bcftools_stats_region = file("${qc_variant_dir}/published_files_manifest/bcftools_stats_regions.csv").exists() ? file("${qc_variant_dir}/published_files_manifest/bcftools_stats_regions.csv") : []
        bcftools_stats_region: bcftools_stats_region
        bcftools_idxstats: file("${qc_variant_dir}/published_files_manifest/bcftools_idxstats.csv", checkIfExists: true)
        bcftools_stats: file("${qc_variant_dir}/published_files_manifest/bcftools_stats.csv", checkIfExists: true)
        charr_contamination: file("${qc_variant_dir}/published_files_manifest/charr_contamination.csv", checkIfExists: true)
    }

    var_collected_manifests_ch = var_manifests_ch
        .getNames()
        .collectEntries { name ->
            def collected_ch = var_manifests_ch
                .getProperty(name)
                .collectFile(
                    name: "${name}_inputs.tsv",
                    keepHeader: true,
                )

            [(name): collected_ch]
        }

    // -- GENERATE SOMALIER REPORTS --
    // Generate a PED file from sample IDs if one is not provided
    if (params.ped_file) {
        ped_file = file(params.ped_file, checkIfExists: true)
    }
    else {
        //make a PED file from the BAM manifest
        MAKE_PED(align_collected_manifests_ch.somalier, params.project_id)
        ped_file = MAKE_PED.out.ped_file
    }

    // Perform Somalier analyses on collected samples
    somalier_profiles_input = align_collected_manifests_ch.somalier
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            [[id: params.project_id], row.output_file]
        }
        .groupTuple()
    somalier_relate_input = somalier_profiles_input.combine(ped_file)
    SOMALIER_RELATE(somalier_relate_input, [])
    SOMALIER_ANCESTRY(somalier_profiles_input, ['ancestry_ref_data', somalier_labels, somalier_1kg_data])

    somalier_ancestry_tsv = SOMALIER_ANCESTRY.out.tsv.map { _meta, tsv -> tsv }
    somalier_relate_pairs_tsv = SOMALIER_RELATE.out.pairs_tsv.map { _meta, tsv -> tsv }
    somalier_relate_samples_tsv = SOMALIER_RELATE.out.samples_tsv.map { _meta, tsv -> tsv }

    // -- GENERATE ALIGN AND VARIANTS QC SUMMARY --
    SUMMARISE_ALIGN_QC(
        params.project_id,
        params.report_target_depths,
        params.report_target_baseq,
        align_qc_flags,
        somalier_ancestry_tsv,
        somalier_relate_pairs_tsv,
        somalier_relate_samples_tsv,
        align_collected_manifests_ch.mosdepth,
        align_collected_manifests_ch.samtools_idxstats,
        align_collected_manifests_ch.samtools_stats,
        align_collected_manifests_ch.verifybamid,
        align_collected_manifests_ch.qrs,
    )

    SUMMARISE_VARIANT_QC(
        params.project_id,
        params.report_var_qual,
        params.report_var_dp,
        var_qc_flags,
        var_collected_manifests_ch.bcftools_stats_region,
        var_collected_manifests_ch.bcftools_idxstats,
        var_collected_manifests_ch.bcftools_stats,
        var_collected_manifests_ch.charr_contamination,
    )

    // -- COLLECT QC METRICS AND MAKE PROJECT QC CONFIG --
    PREPARE_QC_REPORT(
        params.project_id,
        SUMMARISE_ALIGN_QC.out.datavzrd_yml,
        SUMMARISE_ALIGN_QC.out.qc_report_tsv,
        SUMMARISE_ALIGN_QC.out.qc_flags,
        SUMMARISE_ALIGN_QC.out.qc_criteria,
        SUMMARISE_VARIANT_QC.out.datavzrd_yml,
        SUMMARISE_VARIANT_QC.out.qc_report_tsv,
        SUMMARISE_VARIANT_QC.out.qc_flags,
        SUMMARISE_VARIANT_QC.out.qc_criteria,
    )

    // -- GENERATE DATAVZRD REPORT --
    // First need to collect all necessary input files into a single list
    datavzrd_input_files = SUMMARISE_ALIGN_QC.out.datavzrd_yml
        .mix(
            SUMMARISE_ALIGN_QC.out.qc_report_tsv,
            SUMMARISE_ALIGN_QC.out.qc_criteria,
            SUMMARISE_VARIANT_QC.out.qc_report_tsv,
            SUMMARISE_VARIANT_QC.out.qc_criteria,
            PREPARE_QC_REPORT.out.report_js,
            PREPARE_QC_REPORT.out.qc_summary_tsv,
            PREPARE_QC_REPORT.out.flags_summary_tsv,
            PREPARE_QC_REPORT.out.flags_tsv,
        )
        .collect()
        .toList()

    datavzrd_input_ch = PREPARE_QC_REPORT.out.datavzrd_yml
        .combine(datavzrd_input_files)
        .map { yml_file, input_files ->
            [[id: params.project_id], yml_file, input_files]
        }

    DATAVZRD(datavzrd_input_ch)

    align_manifests_out_ch = align_collected_manifests_ch.values().inject(channel.empty()) { acc, ch -> acc.mix(ch) }
    var_manifests_out_ch = var_collected_manifests_ch.values().inject(channel.empty()) { acc, ch -> acc.mix(ch) }

    emit:
    ped_file = ped_file
    align_manifests = align_manifests_out_ch
    var_manifests = var_manifests_out_ch

    somalier_relate_pairs = SOMALIER_RELATE.out.pairs_tsv
    somalier_relate_samples = SOMALIER_RELATE.out.samples_tsv
    somalier_relate_html = SOMALIER_RELATE.out.html
    somalier_ancestry_tsv = SOMALIER_ANCESTRY.out.tsv
    somalier_ancestry_html = SOMALIER_ANCESTRY.out.html

    align_report_tsv = SUMMARISE_ALIGN_QC.out.qc_report_tsv
    align_report_json = SUMMARISE_ALIGN_QC.out.qc_report_json
    align_qc_flags = SUMMARISE_ALIGN_QC.out.qc_flags
    align_qc_criteria = SUMMARISE_ALIGN_QC.out.qc_criteria

    var_report_tsv = SUMMARISE_VARIANT_QC.out.qc_report_tsv
    var_report_json = SUMMARISE_VARIANT_QC.out.qc_report_json
    var_qc_flags = SUMMARISE_VARIANT_QC.out.qc_flags
    var_qc_criteria = SUMMARISE_VARIANT_QC.out.qc_criteria

    project_qc_datavzrd_yml = PREPARE_QC_REPORT.out.datavzrd_yml
    project_qc_summary_tsv = PREPARE_QC_REPORT.out.qc_summary_tsv
    project_qc_flags_summary = PREPARE_QC_REPORT.out.flags_summary_tsv
    project_qc_flags_tsv = PREPARE_QC_REPORT.out.flags_tsv

    datavzrd_report = DATAVZRD.out.report
}
