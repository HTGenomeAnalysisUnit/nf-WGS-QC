include { ECHTVAR_ANNO      } from '../../modules/nf-core/echtvar/anno'
include { BCFTOOLS_INDEX    } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_NORM     } from '../../modules/nf-core/bcftools/norm'
include { COMPUTE_CHARR     } from '../../modules/local/computecharr'

workflow VCF_CONTAMINATION {
    take:
        vcf_files // [meta, vcf, tbi]
        annotation_zip // path(annotation_zip)
        genome_fasta // path(genome_fasta)

    main:
        if (params.annotate_af) {
            // TODO NEED TO IMPLEMENT BCFTOOLS NORM HERE BEFORE PASSING FILES TO ECHTVAR_ANNO
            BCFTOOLS_NORM(vcf_files, ['genome_fasta', genome_fasta])
            // ECHTVAR_ANNO expects [meta, vcf], without the existing index.
            ECHTVAR_ANNO(BCFTOOLS_NORM.out.vcf, [[id: 'af_annotation'], [annotation_zip]], 'af_annotated.vcf.gz')
            BCFTOOLS_INDEX(ECHTVAR_ANNO.out.vcf)
            compute_charr_input = ECHTVAR_ANNO.out.vcf.join(BCFTOOLS_INDEX.out.index)
        } else {
            compute_charr_input = vcf_files
        }
        COMPUTE_CHARR(compute_charr_input, params.charr_af_field)

    emit:
        charr_stats = COMPUTE_CHARR.out.charr_stats
}
