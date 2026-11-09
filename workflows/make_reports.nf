include { samplesheetToList         } from 'plugin/nf-schema'
include { MAKE_PED } from '../modules/local/makeped'
include { SOMALIER_RELATE } from '../modules/local/somalier/relate'
include { SOMALIER_ANCESTRY } from '../modules/local/somalier/ancestry'

workflow MAKE_REPORTS {
    main:
    somalier_labels = file(params.somalier_labels, checkIfExists: true)
    somalier_1kg_data = file(params.somalier_1kgdata, checkIfExists: true)
    // align_rules_file = file(params.alignment_qc_criteria)
    // var_rules_file = file(params.variant_qc_criteria)

    // Collect manifest files across input folders
    // Alignment QC manifests
    align_manifests_ch = channel
            .fromList(samplesheetToList(params.alignment_qc_inputs, "${projectDir}/assets/schema_input_reports.json"))
            .multiMap{ qc_alignment_dir ->
                mosdepth: file("${qc_alignment_dir}/published_files_manifest/coverage_summary.csv", checkIfExists: true)
                qrs_sample_runs: file("${qc_alignment_dir}/published_files_manifest/qrs_run_files.csv", checkIfExists: true)
                samtools_idxstats: file("${qc_alignment_dir}/published_files_manifest/samtools_idxstats.csv", checkIfExists: true)
                samtools_stats: file("${qc_alignment_dir}/published_files_manifest/samtools_stats.csv", checkIfExists: true)
                somalier: file("${qc_alignment_dir}/published_files_manifest/somalier_extract.csv", checkIfExists: true)
                verifybamid: file("${qc_alignment_dir}/published_files_manifest/verifybamid_result.csv", checkIfExists: true)
                bam_manifest: file("${qc_alignment_dir}/input_bam_manifest.tsv", checkIfExists: true)
            }

    align_collected_manifests_ch = align_manifests_ch.getNames()
        .collectEntries { name ->
            def collected_ch = align_manifests_ch.getProperty(name)
                .collectFile(
                    name: "${name}_inputs.tsv",
                    storeDir: "${params.outdir}/qc_reports/input_files/align_qc",
                    keepHeader: true)

            [(name): collected_ch]
        }

    // Variants QC manifests
    var_manifests_ch = channel
            .fromList(samplesheetToList(params.variant_qc_inputs, "${projectDir}/assets/schema_input_reports.json"))
            .multiMap{ qc_variant_dir ->
                def bcftools_stats_region = file("${qc_variant_dir}/published_files_manifest/bcftools_stats_regions.csv").exists() ? file("${qc_variant_dir}/published_files_manifest/bcftools_stats_regions.csv") : null
                bcftools_stats_region: bcftools_stats_region
                bcftools_idxstats: file("${qc_variant_dir}/published_files_manifest/bcftools_idxstats.csv", checkIfExists: true)
                bcftools_stats: file("${qc_variant_dir}/published_files_manifest/bcftools_stats.csv", checkIfExists: true)
                charr_contamination: file("${qc_variant_dir}/published_files_manifest/charr_contamination.csv", checkIfExists: true)
            }

    var_collected_manifests_ch = var_manifests_ch.getNames()
        .collectEntries { name ->
            def collected_ch = var_manifests_ch.getProperty(name)
                .collectFile(
                    name: "${name}_inputs.tsv",
                    storeDir: "${params.outdir}/qc_reports/input_files/variant_qc",
                    keepHeader: true)

            [(name): collected_ch]
        }

    // -- GENERATE SOMALIER REPORTS --
    // Generate a PED file from sample IDs if one is not provided
    if (params.ped_file) {
        ped_file = file(params.ped_file, checkIfExists: true)
    } else {
        //make a PED file from the BAM manifest
        MAKE_PED(align_collected_manifests_ch.somalier, params.project_id)
        ped_file = MAKE_PED.out.ped_file
    }

    // Perform Somalier analyses on collected samples
    somalier_profiles_input = align_collected_manifests_ch.somalier
        .splitCsv(sep: '\t', header: true)
        .map{ row ->
            [[id: params.project_id], row.output_file] 
        }
        .groupTuple()
    somalier_relate_input = somalier_profiles_input.combine(ped_file)
    SOMALIER_RELATE(somalier_relate_input, [])
    SOMALIER_ANCESTRY(somalier_profiles_input, ['ancestry_ref_data', somalier_labels, somalier_1kg_data])

    // -- GENERATE THE PROJECT QC SUMMARY --

    /*
    //Make channels per tranche for each data type
    if (params.by_tranche) {
    tranches_ch = Channel
        .fromPath("${params.main_qc_dir}/*", type: 'dir')
        .multiMap { 
            somalier: tuple(it.baseName, file("$it/qc_alignment/somalier/*"))
            coverage: tuple(it.baseName, file("$it/qc_alignment/coverage/**"))
            sam_stats: tuple(it.baseName, file("$it/qc_alignment/mapflag_stats/*"))
            var_stats: tuple(it.baseName, file("$it/qc_variants/variant_stats/*.vars.stats"), file("$it/qc_variants/variant_stats/*.varsPerContig.tsv"))
            full_tranche: tuple(it.baseName, it)
        }
    } else {
    tranches_ch = Channel
        .fromPath("${params.main_qc_dir}")
        .multiMap { 
            somalier: tuple(params.project_id, file("$it/qc_alignment/somalier/*"))
            coverage: tuple(params.project_id, file("$it/qc_alignment/coverage/**"))
            sam_stats: tuple(params.project_id, file("$it/qc_alignment/mapflag_stats/*"))
            var_stats: tuple(params.project_id, file("$it/qc_variants/variant_stats/*.vars.stats"), file("$it/qc_variants/variant_stats/*.varsPerContig.tsv"))
            full_tranche: tuple(params.project_id, it)
        }
    }

    SOMALIER_RELATE(tranches_ch.somalier, ped_file)
    SOMALIER_ANCESTRY(tranches_ch.somalier, somalier_labels, somalier_1kg_data)

    //Report for alignment QC
    if (params.run_multiqc){
        qc_files = tranches_ch.coverage
                .join(tranches_ch.sam_stats)
                .join(SOMALIER_RELATE.out)
                .join(SOMALIER_ANCESTRY.out.ancestry_tsv)
                //.groupTuple()
                //.view()
            
        MULTIQC(qc_files, file("${projectDir}/conf/multiqc_config.yaml"))
    }

    //Generate summarised data per tranche
    if (params.make_summary) {
    full_tranche_data = tranches_ch.full_tranche
        .join(SOMALIER_RELATE.out)
        .join(SOMALIER_ANCESTRY.out.ancestry_tsv)
        .map { tuple(it[1],it[2],it[3]) }
        .flatten()
        .collect()
        
    //Update summary and app
    PREPARE_SUMMARY(full_tranche_data, align_rules_file, var_rules_file, align_warn_file, var_warn_file, var_qc_rds, align_qc_rds)
    }
    */

    emit:
    ped_file = ped_file

    somalier_relate_pairs = SOMALIER_RELATE.out.pairs_tsv
    somalier_relate_samples = SOMALIER_RELATE.out.samples_tsv
    somalier_relate_html = SOMALIER_RELATE.out.html
    somalier_ancestry_tsv = SOMALIER_ANCESTRY.out.tsv
    somalier_ancestry_html = SOMALIER_ANCESTRY.out.html
}

//FIX SEX COMPUTE. 
//MALE = X_het < 5 && scaled_X_depth < 1.2
//PLOIDY CHECK AS MEAN_COV Y / MEAN_COV X

//REPORT N FAIL SAMPLES FOR EACH FAIL REASON IN PROJECT SUMMARY REPORT.
//CAN WE COUNT TOTAL SEQUENCED BASES FROM THE BAM FILES?
