include { completionSummary } from "./modules/local/pipeline_utils"
include { validateParameters; paramsSummaryLog; paramsHelp } from 'plugin/nf-schema'
include { ALIGNMENTS_QC } from "./workflows/qc_alignments"
include { VARIANTS_QC } from "./workflows/qc_variants"
// include { MAKE_REPORTS_WF } from "${projectDir}/workflows/make_reports"

workflow {
    main:
        // Print help message if requested
        if (params.help) {
            log.info paramsHelp()
            exit 0
        }

        // Log pipeline parameters
        log.info paramsSummaryLog(
            monochrome_logs: params.monochrome_logs
        )
        
        // Validate universal pipeline parameters
        validateParameters(parameters_schema: 'assets/schema_shared_params.json', monochrome_logs: params.monochrome_logs)
        
        // -- BAM FILES QC --
        if (params.input_bam) {   
            validateParameters(parameters_schema: 'assets/schema_alignment_qc_params.json', monochrome_logs: params.monochrome_logs)
            ALIGNMENTS_QC() 
        } 

        // -- VCF FILES QC --
        if (params.input_vcf) {
            validateParameters(parameters_schema: 'assets/schema_variant_qc_params.json', monochrome_logs: params.monochrome_logs)
            VARIANTS_QC()
        }

        // -- CLOSING OPERATIONS --
        // Store complete input params in a YAML file
        def pipeline_inputs_dir = "${params.outdir}/pipeline_inputs"
        channel
            .fromList(params.entrySet())
            .map { entry -> "${entry.key}: ${entry.value}" }
            .collectFile(name: "params.${params.trace_report_suffix}.yml", storeDir: pipeline_inputs_dir, newLine: true)

        // Store tool versions in a YAML file
        channel.topic('versions')
            .map { process, tool, version ->
                [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
            }
            .groupTuple(by: 0)
            .map { process, tool_versions ->
                tool_versions.unique().sort()
                "${process}:\n${tool_versions.join('\n')}"
            }
            .collectFile(name: "collated_versions.${params.trace_report_suffix}.yml", newLine: true, storeDir: pipeline_inputs_dir)

        // Log completion status
        workflow.onComplete {
            completionSummary()
        }

    publish:
        // -- ALIGNMENTS QC OUTPUTS --
        // coverage outputs
        coverage_main_files = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.coverage_global.mix(
                ALIGNMENTS_QC.out.coverage_summary, 
                ALIGNMENTS_QC.out.coverage_regions
            )
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
            
        coverage_detailed_files = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.coverage_regions_bed.mix(
                ALIGNMENTS_QC.out.coverage_per_base_d4,
                ALIGNMENTS_QC.out.coverage_per_base_bed,
                ALIGNMENTS_QC.out.coverage_quantized_bed, 
                ALIGNMENTS_QC.out.coverage_thresholds_bed
            )
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
        coverage_index_files =  (
            params.input_bam 
            ? ALIGNMENTS_QC.out.coverage_per_base_csi.mix(
                ALIGNMENTS_QC.out.coverage_regions_csi,
                ALIGNMENTS_QC.out.coverage_quantized_csi,
                ALIGNMENTS_QC.out.coverage_thresholds_csi
            )
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }

        // samtools outputs
        samtools_stats = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.samtools_stats
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
        samtools_idxstats = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.samtools_idxstats
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }

        // somalier outputs
        somalier_extract = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.somalier_extract
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }

        // quick read stats output
        qrs_run_files = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.qrs_per_sample_run.mix(ALIGNMENTS_QC.out.qrs_rundetails)
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
        qrs_read_stats = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.qrs_read_stats
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }

        // verifybamid outputs
        verifybamid_result = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.verifybamid_result
            : channel.empty()
        ).map{ meta, output_file, log_file -> meta + [output_file: output_file, log_file: log_file] }

        // bam manifest file
        bam_manifest = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.bam_manifest
            : channel.empty()
        )
        sample_to_user_map_bam = (
            params.input_bam 
            ? ALIGNMENTS_QC.out.sample_to_user_map
            : channel.empty()
        )

        // -- VCF QC OUTPUTS --
        charr_stats = (
            params.input_vcf 
            ? VARIANTS_QC.out.charr_stats
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
        bcftools_stats = (
            params.input_vcf 
            ? VARIANTS_QC.out.bcftools_stats
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
        bcftools_stats_regions = (
            params.input_vcf 
            ? VARIANTS_QC.out.bcftools_stats_regions
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
        bcftools_idxstats = (
            params.input_vcf 
            ? VARIANTS_QC.out.bcftools_idxstats
            : channel.empty()
        ).map{ meta, output_file -> meta + [output_file: output_file] }
        vcf_manifest = (
            params.input_vcf 
            ? VARIANTS_QC.out.vcf_manifest
            : channel.empty()
        )
        sample_to_user_map_vcf = (
            params.input_vcf 
            ? VARIANTS_QC.out.sample_to_user_map
            : channel.empty()
        )
}

output {
    coverage_main_files {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/coverage/summary'
        index {
            path 'qc_alignments/published_files_manifest/coverage_summary.csv'
            header true
            sep '\t'
        }
    }
    coverage_detailed_files {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/coverage/detailed'
        index {
            path 'qc_alignments/published_files_manifest/coverage_detailed.csv'
            header true
            sep '\t'
        }
    }
    coverage_index_files {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/coverage/detailed'
    }
    samtools_stats {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/samtools/stats'
        index {
            path 'qc_alignments/published_files_manifest/samtools_stats.csv'
            header true
            sep '\t'
        }
    }
    samtools_idxstats {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/samtools/idxstats'
        index {
            path 'qc_alignments/published_files_manifest/samtools_idxstats.csv'
            header true
            sep '\t'
        }
    }
    somalier_extract {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/somalier/extract'
        index {
            path 'qc_alignments/published_files_manifest/somalier_extract.csv'
            header true
            sep '\t'
        }
    }
    qrs_run_files {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/qrs/run_files'
        index {
            path 'qc_alignments/published_files_manifest/qrs_run_files.csv'
            header true
            sep '\t'
        }
    }
    qrs_read_stats {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/qrs/read_stats'
        index {
            path 'qc_alignments/published_files_manifest/qrs_read_stats.csv'
            header true
            sep '\t'
        }
    }
    verifybamid_result {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments/verifybamid'
        index {
            path 'qc_alignments/published_files_manifest/verifybamid_result.csv'
            header true
            sep '\t'
        }
    }
    bam_manifest {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments'
    }
    sample_to_user_map_bam {
        enabled (params.input_bam ? true : false)
        path 'qc_alignments'
    }

    //-- VARIANT QC OUTPUTS --
    charr_stats {
        enabled (params.input_vcf ? true : false)
        path 'qc_variants/charr_contamination'
        index {
            path 'qc_variants/published_files_manifest/charr_contamination.csv'
            header true
            sep '\t'
        }
    }
    bcftools_stats {
        enabled (params.input_vcf ? true : false)
        path 'qc_variants/bcftools/stats'
        index {
            path 'qc_variants/published_files_manifest/bcftools_stats.csv'
            header true
            sep '\t'
        }
    }
    bcftools_stats_regions {
        enabled (params.input_vcf ? true : false)
        path 'qc_variants/bcftools/stats'
        index {
            path 'qc_variants/published_files_manifest/bcftools_stats_regions.csv'
            header true
            sep '\t'
        }
    }
    bcftools_idxstats {
        enabled (params.input_vcf ? true : false)
        path 'qc_variants/bcftools/idxstats'
        index {
            path 'qc_variants/published_files_manifest/bcftools_idxstats.csv'
            header true
            sep '\t'
        }
    }
    vcf_manifest {
        enabled (params.input_vcf ? true : false)
        path 'qc_variants'
    }
    sample_to_user_map_vcf {
        enabled (params.input_vcf ? true : false)
        path 'qc_variants'
    }
}

workflow MAKE_REPORTS {
    // MAKE_REPORTS_WF()
}
