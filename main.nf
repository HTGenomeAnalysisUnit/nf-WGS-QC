nextflow.enable.dsl=2

outdir = params.outdir ? params.outdir : "QC"
file(outdir).mkdirs()
outdir_align = "${outdir}/qc_alignment"
outdir_vars = "${outdir}/qc_variants"

include { RUN_ALIGNEMENTS_QC } from "${projectDir}/workflows/qc_alignments" addParams( outdir: "${outdir_align}")
include { RUN_VARIANTS_QC } from "${projectDir}/workflows/qc_variants" addParams( outdir: "${outdir_vars}")

workflow {
        log.info """\
    =============================================
         WGS QC    -    N F   P I P E L I N E    
    =============================================
    bam files   : ${params.bamfiles}
    VCF files   : ${params.vcffiles}
    VCF index   : ${params.vcf_index_type}
    mode        : ${params.mode}
    ped_file    : ${params.ped_file}
    genome      : ${params.genome}
    regions     : ${params.regions}
    outdir      : ${params.outdir}
    =============================================
    """
    .stripIndent()

    if (params.bamfiles) {
        bamfiles = Channel
            .fromPath(params.bamfiles, checkIfExists: true)
            .splitCsv(sep: '\t', header: true)
            .map{ tuple(it.sample,file(it.bam, checkIfExists: true),file("${it.bam}.bai", checkIfExists: true)) }

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
        ped_file = file(params.ped_file)
        genome_ref = tuple(file(params.genome), file("${params.genome}.fai"))

        RUN_ALIGNEMENTS_QC(bamfiles, fastp_files, samblaster_files, params.regions, somalier_varset, somalier_labels, somalier_1kgdata, ped_file, genome_ref) 

        somalier_ancestry = RUN_ALIGNEMENTS_QC.out.somalier_ancestry
    } else {
        somalier_ancestry = file('noancestry')
    }

    if (params.vcffiles) {
        vcf_files = Channel
        .fromPath(params.vcffiles, checkIfExists: true)
        .splitCsv(sep: '\t', header: true)
        .map{ it -> tuple(it.sample,file(it.vcf, checkIfExists: true),file("${it.vcf}.${params.vcf_index_type}", checkIfExists: true)) }

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