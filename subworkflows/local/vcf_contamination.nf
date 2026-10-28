include { ECHTVAR_ANNO } from '../../modules/nf-core/echtvar/anno'
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index'
include { COMPUTE_CHARR } from '../../modules/local/computecharr'

workflow VCF_CONTAMINATION {
    take:
        vcf_files // [meta, vcf, tbi]
        annotation_zip // path(annotation_zip)

    main:
        if (params.annotate_af) {
            // ECHTVAR_ANNO expects [meta, vcf], without the existing index.
            vcf_to_annotate = vcf_files.map { meta, vcf, _index -> tuple(meta, vcf) }
            ECHTVAR_ANNO(vcf_to_annotate, [[id: 'af_annotation'], [annotation_zip]], 'af_annotated.vcf.gz')
            BCFTOOLS_INDEX(ECHTVAR_ANNO.out.vcf)
            compute_charr_input = ECHTVAR_ANNO.out.vcf.join(BCFTOOLS_INDEX.out.index)
        } else {
            compute_charr_input = vcf_files
        }
        COMPUTE_CHARR(compute_charr_input, params.charr_af_field)

    emit:
        charr_stats = COMPUTE_CHARR.out.charr_stats
}
