nextflow.enable.dsl=2

main_outdir = params.outdir ? params.outdir : "QC"
outdir = params.project_id ? "${main_outdir}/${params.project_id}" : "${main_outdir}"
file(outdir).mkdirs()
outdir_align = "${outdir}/qc_alignment"
outdir_vars = "${outdir}/qc_variants"
outdir_bamindex = "${outdir}/bam_index"
outdir_vcfindex = "${outdir}/vcf_index"

include { RUN_ALIGNEMENTS_QC } from "${projectDir}/workflows/qc_alignments" addParams( outdir: "${outdir_align}")
include { RUN_VARIANTS_QC } from "${projectDir}/workflows/qc_variants" addParams( outdir: "${outdir_vars}")
include { GET_SAMPLE_ID_BAM; GET_SAMPLE_ID_VCF } from "${projectDir}/modules/get_sample_id"
include { GET_BAM_MANIFEST; GET_VCF_MANIFEST } from "${projectDir}/modules/get_sample_id"
include { MAKE_REPORTS_WF } from "${projectDir}/workflows/make_reports"
include { SAMTOOLS_INDEX } from "${projectDir}/modules/samtools" addParams( outdir: "${outdir_bamindex}")
include { BCFTOOLS_INDEX } from "${projectDir}/modules/bcftools" addParams( outdir: "${outdir_vcfindex}")

workflow {
        log.info """\
    =============================================
         WGS QC    -    N F   P I P E L I N E    
    =============================================
    bam files   : ${params.bamfiles ? params.bamfiles : "NO BAM files"}
    VCF files   : ${params.vcffiles ? params.vcffiles : "NO VCF files"}
    VCF index   : ${params.vcf_index_type}
    mode        : ${params.mode}
    genome      : ${params.genome}
    ${params.regions ? "regions: ${params.regions}" : ''}
    outdir      : ${outdir}
    =============================================
    """
    .stripIndent()

    // Set input files channels
    //bamfiles/vcffiles is a glob of multiple files
    bam_input_ch = Channel.empty()
    vcf_input_ch = Channel.empty()
    if (params.bamlist) {
        bam_input_paths = file(params.bamlist)
            .readLines()
        bam_input_ch = Channel.fromList(bam_input_paths)
            .map { it -> file("${it.split('\t')[1]}") }  
    }
    if (params.vcflist) {
        vcf_input_paths = file(params.vcflist)
            .readLines()
        vcf_input_ch = Channel.fromList(vcf_input_paths)
            .map { it -> file("${it.split('\t')[1]}") }  
    }
    if (params.bamfiles) {
        bam_input_ch = Channel
            .watchPath(params.bamfiles, 'create,modify')
    }
    if (params.vcffiles) {
        vcf_input_ch = Channel
            .watchPath(params.vcffiles, 'create,modify')
    }

    //Pre-process: create indexes, get sample ids, make manifest
    if (params.make_bam_index) {
        bam_files = SAMTOOLS_INDEX(bam_input_ch)
    } else {
        bam_files = bam_input_ch
            .map { bam_file -> tuple(bam_file, file("${bam_file}.bai", checkIfExists: true) }
    }
    GET_SAMPLE_ID_BAM(bam_files)
    GET_BAM_MANIFEST(bam_files)

    if (params.make_vcf_index) {
        vcf_files = BCFTOOLS_INDEX(vcf_input_ch)
    } else {
        vcf_files = vcf_input_ch
            .map { vcf_file -> tuple(vcf_file, file("${vcf_file}.${params.vcf_index_type}", checkIfExists: true) }
    }
    GET_SAMPLE_ID_VCF(vcf_files)
    GET_VCF_MANIFEST(vcf_files)

    // Get sample IDs from BAM if available otherwise from VCF files
    if (params.bamfiles || params.bamlist) {
        manifest_bam = file("${outdir}/bam_manifest.tsv")
        ext2int_map_ch = GET_BAM_MANIFEST.out
            .map { "${it.join('\t')}" }
            .subscribe { manifest_bam.append("$it\n") }
            //.collectFile(name: "${outdir}/ext2int_map_bam.tsv", newLine: true)
    }
    if (params.vcffiles || params.vcflist) {
        manifest_vcf = file("${outdir}/vcf_manifest.tsv")
        ext2int_map_ch = GET_VCF_MANIFEST.out
            .map { "${it.join('\t')}" }
            .subscribe { manifest_vcf.append("$it\n") }
            //.collectFile(name: "${outdir}/ext2int_map_vcf.tsv", newLine: true)
    }

    //Perform the QC on BAMs and/or VCFs
    if (params.bamfiles || params.bamlist) {
        checkPathParamList = [
                params.somalier_vars, params.genome, "${params.genome}.fai"
                ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
        
        somalier_varset = file(params.somalier_vars)
        genome_ref = tuple(file(params.genome), file("${params.genome}.fai"))

        RUN_ALIGNEMENTS_QC(GET_SAMPLE_ID_BAM.out.bam_files, params.regions, somalier_varset, genome_ref) 
    } 

    if (params.vcffiles || params.vcflist) {
        annotation_zip_file = file(params.af_annotation_zip, checkIfExists: true)
        RUN_VARIANTS_QC(GET_SAMPLE_ID_VCF.out.vcf_files, annotation_zip_file)
    }
}

workflow MAKE_REPORTS {
    MAKE_REPORTS_WF()
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

