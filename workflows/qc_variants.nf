include { samplesheetToList } from 'plugin/nf-schema'
include { BCFTOOLS_STATS    } from "../modules/local/bcftools/stats"
include { BCFTOOLS_STATS as BCFTOOLS_STATS_WITH_REGIONS } from "../modules/local/bcftools/stats"
include { BCFTOOLS_IDXSTATS } from "../modules/local/bcftools/idxstats"
include { BCFTOOLS_INDEX    } from "../modules/nf-core/bcftools/index"
include { BCFTOOLS_VIEW     } from '../modules/nf-core/bcftools/view'
include { VCF_MANIFEST      } from "../modules/local/vcfmanifest"
include { VCF_CONTAMINATION } from "../subworkflows/local/vcf_contamination"

workflow VARIANTS_QC {
    main:
        annotation_zip_file = file(params.af_annotation_zip, checkIfExists: true)
        genome_fasta        = file(params.genome_fasta, checkIfExists: true)

        // Configure input channels
        input_ch = channel
            .fromList(samplesheetToList(params.input_vcf, "${projectDir}/assets/schema_input_vcf.json"))
            .multiMap{ meta, vcf_file, samplesheet_regions_bed -> 
                // Check the filename string: Path.endsWith compares whole path components.
                def file_type = vcf_file.name.endsWith('.g.vcf') || vcf_file.name.endsWith('.g.vcf.gz') ? 'gvcf' : 'vcf'
                def m = [ id: meta.id as String, sex: meta.sex as String, batch_id: meta.batch_id as String, file_type: file_type ] // ensure values are treated as string
                // When params regions_bed is provided, this is used when a regions BED file is not specified in the samplesheet
                def regions_bed_file = samplesheet_regions_bed ? samplesheet_regions_bed : params.regions_bed ? file(params.regions_bed, checkIfExists: true) : []
                vcf: tuple(m, vcf_file, files("${vcf_file}.{csi,tbi}")) 
                regions: tuple(m, regions_bed_file)
            }

        // Create indexes where necessary
        vcf_processing_ch = input_ch.vcf
            .branch{ meta, vcf, csi ->
                indexed: csi
                    return [meta, vcf, file(csi[0], checkIfExists: true)]
                to_index: !csi
                    return [meta, vcf]
            }
        BCFTOOLS_INDEX(vcf_processing_ch.to_index)
        
        vcf_files_with_created_index = vcf_processing_ch.to_index.join(BCFTOOLS_INDEX.out.index)
        vcf_files_indexed = vcf_files_with_created_index.mix(vcf_processing_ch.indexed)
        
        // now filter variants to generate the subset suitable for stats and contamination analysis
        BCFTOOLS_VIEW(vcf_files_indexed, [], [], [])
        processed_vcf_ch = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.index)

        // get sample ids as stored in the VCF file and make a manifest file
        // This is a very cheap operation so we collect all files and make it in a single process
        vcf_manifest_input_ch = processed_vcf_ch.map{meta, vcf, csi -> [meta.id, vcf, csi]}.toList().transpose().toList()
        VCF_MANIFEST(vcf_manifest_input_ch)

        // take the sample id extracted from vcf files and make them into a channel
        vcf_id_2_user_id = VCF_MANIFEST.out.sample_to_user_map
            .splitCsv(sep: '\t', skip: 1, header: false)

        // Here we run each component and then add the VCF sample ID to the metadata for downstream analysis
        // This metadata operation is made here to avoid processing waiting for the manifest extraction to complete
        // -- contamination analysis --
        VCF_CONTAMINATION(processed_vcf_ch, annotation_zip_file, genome_fasta)
        vcf_contaimination_out = VCF_CONTAMINATION.out.charr_stats
            .map{meta, charr_stats -> [meta.id, meta, charr_stats]}
            .join(vcf_id_2_user_id)
            .map{ _sample_id, meta, charr_stats, vcf_sample_id ->
                def m = meta + [vcf_sample_id: vcf_sample_id as String]
                return [m, charr_stats]
            }
        // -- bcftools stats --
        // Include the regions file for stats analysis
        bcftools_stats_input = processed_vcf_ch.join(input_ch.regions)
            .branch{ _meta, _vcf, _csi, regions_bed_file ->
                with_regions: regions_bed_file
                without_regions: !regions_bed_file
            }
        // Run stats on specific regions if provided
        BCFTOOLS_STATS_WITH_REGIONS(bcftools_stats_input.with_regions, ['targets', []], ['samples', []], ['exons', []], ['fasta', []])
        bcftools_stats_regions_out = BCFTOOLS_STATS_WITH_REGIONS.out.stats
            .map{meta, stats -> [meta.id, meta, stats]}
            .join(vcf_id_2_user_id)
            .map{ _sample_id, meta, stats, vcf_sample_id ->
                def m = meta + [vcf_sample_id: vcf_sample_id as String]
                return [m, stats]
            }
        // Run also on all processed variants
        bcftools_stats_input_all = processed_vcf_ch.join(input_ch.regions)
        BCFTOOLS_STATS(bcftools_stats_input_all, ['targets', []], ['samples', []], ['exons', []], ['fasta', []])
        bcftools_stats_out = BCFTOOLS_STATS.out.stats
            .map{meta, stats -> [meta.id, meta, stats]}
            .join(vcf_id_2_user_id)
            .map{ _sample_id, meta, stats, vcf_sample_id ->
                def m = meta + [vcf_sample_id: vcf_sample_id as String]
                return [m, stats]
            }
        // -- bcftools idxstats --
        BCFTOOLS_IDXSTATS(processed_vcf_ch)
        bcftools_idxstats_out = BCFTOOLS_IDXSTATS.out.idxstats
            .map{meta, idxstats -> [meta.id, meta, idxstats]}
            .join(vcf_id_2_user_id)
            .map{ _sample_id, meta, idxstats, vcf_sample_id ->
                def m = meta + [vcf_sample_id: vcf_sample_id as String]
                return [m, idxstats]
            }
    
    emit:
        // qc module outputs
        charr_stats = vcf_contaimination_out
        bcftools_stats = bcftools_stats_out
        bcftools_stats_regions = bcftools_stats_regions_out
        bcftools_idxstats = bcftools_idxstats_out

        // vcf manifest file
        vcf_manifest = VCF_MANIFEST.out.vcf_manifest
        sample_to_user_map = VCF_MANIFEST.out.sample_to_user_map
}
