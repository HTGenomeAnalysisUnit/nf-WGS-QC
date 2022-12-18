nextflow.enable.dsl=2

main_outdir = params.outdir ? params.outdir : "QC"
outdir = params.project_id ? "${main_outdir}/${params.project_id}" : "${main_outdir}"
file(outdir).mkdirs()
outdir_align = "${outdir}/qc_alignment"
outdir_vars = "${outdir}/qc_variants"

include { RUN_ALIGNEMENTS_QC } from "${projectDir}/workflows/qc_alignments" addParams( outdir: "${outdir_align}")
include { RUN_VARIANTS_QC } from "${projectDir}/workflows/qc_variants" addParams( outdir: "${outdir_vars}")
include { GET_SAMPLE_ID_BAM; GET_SAMPLE_ID_VCF } from "${projectDir}/modules/get_sample_id"
include { GET_BAM_MANIFEST; GET_VCF_MANIFEST } from "${projectDir}/modules/get_sample_id"
include { MAKE_REPORTS_WF } from "${projectDir}/workflows/make_reports"

workflow {
        log.info """\
    =============================================
         WGS QC    -    N F   P I P E L I N E    
    =============================================
    bam files   : ${params.bamfiles ? params.vcffiles : "NO BAM files"}
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
    if (params.bamfiles) {
        bam_files = Channel
            .watchPath(params.bamfiles)
            .map { tuple(it, file("${it}.bai")) }
        GET_SAMPLE_ID_BAM(bam_files)
        GET_BAM_MANIFEST(bam_files)
    }
    if (params.vcffiles) {
        vcf_files = Channel
            .watchPath(params.vcffiles)
            .map{ tuple(it, file("${it}.${params.vcf_index_type}")) }
        GET_SAMPLE_ID_VCF(vcf_files)
        GET_VCF_MANIFEST(vcf_files)
    }

    // Get sample IDs from BAM if available otherwise from VCF files
    if (params.bamfiles) {
        manifest_bam = file("${outdir}/bam_manifest.tsv")
        ext2int_map_ch = GET_BAM_MANIFEST.out
            .map { "${it.join('\t')}" }
            .subscribe { manifest_bam.append("$it\n") }
            //.collectFile(name: "${outdir}/ext2int_map_bam.tsv", newLine: true)
    }
    if (params.vcffiles) {
        manifest_vcf = file("${outdir}/vcf_manifest.tsv")
        ext2int_map_ch = GET_VCF_MANIFEST.out
            .map { "${it.join('\t')}" }
            .subscribe { manifest_vcf.append("$it\n") }
            //.collectFile(name: "${outdir}/ext2int_map_vcf.tsv", newLine: true)
    }

    //Perform the QC on BAMs and/or VCFs
    if (params.bamfiles) {
        checkPathParamList = [
                params.somalier_vars, params.genome, "${params.genome}.fai"
                ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
        
        somalier_varset = file(params.somalier_vars)
        genome_ref = tuple(file(params.genome), file("${params.genome}.fai"))

        RUN_ALIGNEMENTS_QC(GET_SAMPLE_ID_BAM.out.bam_files, params.regions, somalier_varset, genome_ref) 
    } 

    if (params.vcffiles) {
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