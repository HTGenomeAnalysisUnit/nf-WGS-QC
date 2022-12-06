nextflow.enable.dsl=2

main_outdir = params.outdir ? params.outdir : "QC"
outdir = params.project_id ? "${main_outdir}/${params.project_id}" : "${main_outdir}"
file(outdir).mkdirs()
outdir_align = "${outdir}/qc_alignment"
outdir_vars = "${outdir}/qc_variants"

include { RUN_ALIGNEMENTS_QC } from "${projectDir}/workflows/qc_alignments" addParams( outdir: "${outdir_align}")
include { RUN_VARIANTS_QC } from "${projectDir}/workflows/qc_variants" addParams( outdir: "${outdir_vars}")
include { MAKE_PED }  from "${projectDir}/modules/make_ped" addParams( outdir: "${outdir}")
include { GET_SAMPLE_ID_BAM; GET_SAMPLE_ID_VCF } from "${projectDir}/modules/get_sample_id"

workflow {
        log.info """\
    =============================================
         WGS QC    -    N F   P I P E L I N E    
    =============================================
    bam files   : ${params.bamfiles ? params.vcffiles : "NO BAM files"}
    VCF files   : ${params.vcffiles ? params.vcffiles : "NO VCF files"}
    VCF index   : ${params.vcf_index_type}
    mode        : ${params.mode}
    ped_file    : ${params.ped_file ? params.ped_file : "imputed"}
    genome      : ${params.genome}
    regions     : ${params.regions}
    outdir      : ${params.outdir}
    =============================================
    """
    .stripIndent()

    project_metadata = file(params.project_metadata, checkIfExists: true)

    // Set input files channels
    //bamfiles/vcffiles is a glob of multiple files
    if (params.bamfiles) {
        bamfiles = Channel
            .fromPath(params.bamfiles, checkIfExists: true)
            .map { tuple(it, file("${it}.bai", checkIfExists: true)) }
    }
    if (params.vcffiles) {
        vcf_files = Channel
            .fromPath(params.vcffiles, checkIfExists: true)
            .map{ tuple(it, file("${it}.${params.vcf_index_type}", checkIfExists: true)) }
    }

    // Get sample IDs from BAM if available otherwise from VCF files
    header_channel = Channel.of("internal_id\texternal_id")
    if (params.bamfiles) {
        GET_SAMPLE_ID_BAM(bam_files)
        ext2int_map_ch = GET_SAMPLE_ID_BAM.out.id_map
            .collect {it.join("\t")}
            .collectFile(name: "${project_outdir}/ext2int_map.tsv", newLine: true)
    } else if (params.vcffiles) {
        GET_SAMPLE_ID_VCF(bam_files)
        ext2int_map_ch = GET_SAMPLE_ID_BAM.out.id_map
            .collect {it.join("\t")}
            .collectFile(name: "${project_outdir}/ext2int_map.tsv", newLine: true)
    }

    // Generate a PED file from sample IDs if one is not provided
    if (params.ped_file) {
        ped_file = file(params.ped_file)
    } else {
        MAKE_PED(project_metadata, ext2int_map_ch)
        ped_file = MAKE_PED.out
    }

    //Perform the QC on BAMs and/or VCFs
    if (params.bamfiles) {
        checkPathParamList = [
                params.somalier_vars, params.somalier_labels, 
                params.somalier_1kgdata, params.ped_file, params.genome,
                "${params.genome}.fai"
                ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
        
        fastp_files = params.fastp_json ? Channel.fromPath(params.fastp) : Channel.empty()
        samblaster_files = params.samblaster_log ? Channel.fromPath(params.samblaster) : Channel.empty()
        somalier_varset = file(params.somalier_vars)
        somalier_labels = file(params.somalier_labels)
        somalier_1kgdata = file(params.somalier_1kgdata)
        genome_ref = tuple(file(params.genome), file("${params.genome}.fai"))

        RUN_ALIGNEMENTS_QC(GET_SAMPLE_ID_BAM.out.bam_files, fastp_files, samblaster_files, params.regions, somalier_varset, somalier_labels, somalier_1kgdata, ped_file, genome_ref) 

        somalier_ancestry = RUN_ALIGNEMENTS_QC.out.somalier_ancestry
    } else {
        somalier_ancestry = file('noancestry')
    }

    if (params.vcffiles) {
        vars_qc_input_ch = vcf_files.combine(somalier_ancestry)
        annotation_zip_file = file(params.af_annotation_zip, checkIfExists: true)
        RUN_VARIANTS_QC(vcf_files, somalier_ancestry, annotation_zip_file)
    }
}

workflow.onComplete { 
	log.info ( workflow.success ? """\
    
    PIPELINE COMPLETED!    
    ============================================
    QC folder    --> ${params.outdir}
    Alignements  --> ${outdir_align}
    Variants     --> ${outdir_vars}   
    """.stripIndent() : "Oops .. something went wrong" )
}