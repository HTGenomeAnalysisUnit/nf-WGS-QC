include { samplesheetToList     } from 'plugin/nf-schema'
include { BAM_MANIFEST          } from "../modules/local/bammanifest"
include { SAMTOOLS_STATS        } from "../modules/nf-core/samtools/stats"
include { SAMTOOLS_IDXSTATS     } from "../modules/nf-core/samtools/idxstats"
include { SAMTOOLS_INDEX        } from "../modules/nf-core/samtools/index" 
include { MOSDEPTH              } from "../modules/nf-core/mosdepth" 
include { SOMALIER_EXTRACT      } from "../modules/nf-core/somalier/extract" 
include { QUICK_READS_STATS     } from "../modules/local/qrs" 

workflow ALIGNMENTS_QC {        
    main:
        // Make files value
        somalier_varset     = file(params.somalier_vars, checkIfExists: true) //VCF of somalier reference variants
        // somalier_labels     = file(params.somalier_labels, checkIfExists: true) //ancestry labels for 1kg sample 
        // somalier_1kgdata    = file(params.somalier_1kgdata, type: 'dir', checkIfExists: true) //path to a folder containing .somalier files for the 1000G data
        genome_fasta        = file(params.genome_fasta, checkIfExists: true)
        genome_fai          = file("${params.genome_fasta}.fai", checkIfExists: true)

        // Make mosdepth quantize bins into a list and ensure quantize labels are aligned
        mosdepth_quantize_bins_list = params.mosdepth_quantize_bins.split(':')
        mosdepth_quantize_labels_list = params.mosdepth_quantize_labels.split(',')
        if (mosdepth_quantize_bins_list.size() != mosdepth_quantize_labels_list.size()) {
            throw new IllegalArgumentException("The number of mosdepth quantize bins (mosdepth_quantize_bins) must match the number of quantize labels (mosdepth_quantize_labels)")
        }

        // Configure input channels
        input_ch = channel
            .fromList(samplesheetToList(params.input_bam, "${projectDir}/assets/schema_input_bam.json"))
            .multiMap{ meta, bam_file, samplesheet_regions_bed -> 
                def m = [ id: meta.id as String, sex: meta.sex as String, batch_id: meta.batch_id as String ] // ensure values are treated as string
                // When params regions_bed is provided, this is used when a regions BED file is not specified in the samplesheet
                def regions_bed_file = samplesheet_regions_bed ? samplesheet_regions_bed : params.regions_bed ? file(params.regions_bed, checkIfExists: true) : []
                bam: tuple(m, bam_file, files("${bam_file}.{bai,crai}")) 
                regions: tuple(m, regions_bed_file)
            }

        // Create indexes where necessary
        bam_processing_ch = input_ch.bam
            .branch{ meta, bam, bai ->
                indexed: bai
                    return [meta, bam, file(bai[0], checkIfExists: true)]
                to_index: !bai
                    return [meta, bam]
            }
        SAMTOOLS_INDEX(bam_processing_ch.to_index)
        
        bam_files_with_created_index = bam_processing_ch.to_index.join(SAMTOOLS_INDEX.out.index)
        bam_files_indexed = bam_files_with_created_index.mix(bam_processing_ch.indexed)

        // get sample ids as stored in the BAM file and make a manifest file
        // This is a very cheap operation so we collect all files and make it in a single process
        bam_manifest_input_ch = bam_files_indexed.map{meta, bam, bai -> [meta.id, bam, bai]}.toList().transpose().toList()
        BAM_MANIFEST(bam_manifest_input_ch)

        // take the sample id extracted from bam files and make a channel
        bam_id_2_user_id = BAM_MANIFEST.out.sample_to_user_map
            .splitCsv(sep: '\t', skip: 1, header: false)
        
        // Here we run each component and then add the VCF sample ID to the metadata for downstream analysis
        // This metadata operation is made here to avoid processing waiting for the manifest extraction to complete
        // -- MOSDEPTH --
        // Include the regions file for coverage analysis
        mosdepth_input = bam_files_indexed.join(input_ch.regions)
        MOSDEPTH(mosdepth_input, ['genome_fasta', genome_fasta], mosdepth_quantize_labels_list)
        mosdepth_out = MOSDEPTH.out.getNames()
            .findAll { name -> !name.startsWith('versions_') }
            .collectEntries { name ->
                def output_ch = MOSDEPTH.out.getProperty(name)

                def processed_ch = output_ch
                    .map { meta, output_file ->
                        tuple(meta.id, meta, output_file)
                    }
                    .join(bam_id_2_user_id)
                    .map { _sample_id, meta, output_file, bam_sample_id ->
                        def m = meta + [bam_sample_id: bam_sample_id as String]
                        tuple(m, output_file)
                    }

                [(name): processed_ch]
            }

        // -- SAMTOOLS STATS --
        SAMTOOLS_STATS(bam_files_indexed, ['ref_genome', genome_fasta, genome_fai])
        samtools_stats_out = SAMTOOLS_STATS.out.stats
            .map{meta, stats -> [meta.id, meta, stats]}
            .join(bam_id_2_user_id)
            .map{ _sample_id, meta, stats, bam_sample_id ->
                def m = meta + [bam_sample_id: bam_sample_id as String]
                return [m, stats]
            }

        // -- SAMTOOLS IDXSTATS --
        SAMTOOLS_IDXSTATS(bam_files_indexed)
        samtools_idxstats_out = SAMTOOLS_IDXSTATS.out.idxstats
            .map{meta, idxstats -> [meta.id, meta, idxstats]}
            .join(bam_id_2_user_id)
            .map{ _sample_id, meta, idxstats, bam_sample_id ->
                def m = meta + [bam_sample_id: bam_sample_id as String]
                return [m, idxstats]
            }

        // -- SOMALIER EXTRACT --
        SOMALIER_EXTRACT(bam_files_indexed, ['genome_fasta', genome_fasta], ['genome_fai', genome_fai], ['somalier_varset', somalier_varset])
        somalier_extract_out = SOMALIER_EXTRACT.out.extract
            .map{meta, somalier_profile -> [meta.id, meta, somalier_profile]}
            .join(bam_id_2_user_id)
            .map{ _sample_id, meta, somalier_profile, bam_sample_id ->
                def m = meta + [bam_sample_id: bam_sample_id as String]
                return [m, somalier_profile]
            }

        // -- QUICK READS STATS --
        QUICK_READS_STATS(bam_files_indexed)
        qrs_out = QUICK_READS_STATS.out.getNames()
            .findAll { name -> !name.startsWith('versions_') }
            .collectEntries { name ->
                def output_ch = QUICK_READS_STATS.out.getProperty(name)

                def processed_ch = output_ch
                    .map { meta, output_file ->
                        tuple(meta.id, meta, output_file)
                    }
                    .join(bam_id_2_user_id)
                    .map { _sample_id, meta, output_file, bam_sample_id ->
                        def m = meta + [bam_sample_id: bam_sample_id as String]
                        tuple(m, output_file)
                    }

                [(name): processed_ch]
            }

    emit:
        // coverage outputs
        coverage_global = mosdepth_out.global_txt
        coverage_summary = mosdepth_out.summary_txt
        coverage_regions = mosdepth_out.regions_txt
        coverage_per_base_d4 = mosdepth_out.per_base_d4
        coverage_per_base_bed = mosdepth_out.per_base_bed
        coverage_per_base_csi = mosdepth_out.per_base_csi
        coverage_regions_bed = mosdepth_out.regions_bed
        coverage_regions_csi = mosdepth_out.regions_csi
        coverage_quantized_bed = mosdepth_out.quantized_bed
        coverage_quantized_csi = mosdepth_out.quantized_csi
        coverage_thresholds_bed = mosdepth_out.thresholds_bed
        coverage_thresholds_csi = mosdepth_out.thresholds_csi

        // samtools outputs
        samtools_stats = samtools_stats_out
        samtools_idxstats = samtools_idxstats_out

        // somalier outputs
        somalier_extract = somalier_extract_out

        // quick read stats output
        qrs_per_sample_run = qrs_out.per_sample_run
        qrs_rundetails = qrs_out.rundetails
        qrs_read_stats = qrs_out.read_stats

        // bam manifest file
        bam_manifest = BAM_MANIFEST.out.bam_manifest
        sample_to_user_map = BAM_MANIFEST.out.sample_to_user_map
}