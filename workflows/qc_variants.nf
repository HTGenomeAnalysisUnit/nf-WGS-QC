nextflow.enable.dsl=2

// params.outdir = "QC/variants", output folder for QC files */
output_contamination = file("${params.outdir}/contamination_stats")
output_stats = file("${params.outdir}/variant_stats")
output_reports = file("${params.outdir}/reports")

include { BCFTOOLS_STATS } from "${projectDir}/modules/bcftools" addParams( outdir: "${output_stats}")
include { BCFTOOLS_IDX_STATS } from "${projectDir}/modules/bcftools" addParams( outdir: "${output_stats}")
include { VCF_CONTAMINATION } from "${projectDir}/modules/vcf_contamination" addParams( outdir: "${output_contamination}")

workflow RUN_VARIANTS_QC {
    take:
        vcf_files // [sample, file.vcf.gz, file.vcf.gz.{tbi,csi}]
        annotation_zip // zip file of gnomad for echtvar annotation

    main:
        VCF_CONTAMINATION(vcf_files, annotation_zip)
        BCFTOOLS_STATS(vcf_files)
        BCFTOOLS_IDX_STATS(vcf_files)
}