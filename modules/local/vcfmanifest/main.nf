process VCF_MANIFEST {
    tag "${params.project_id}_vcf_manifest"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(sample_ids), path(vcf, stageAs: 'dir??/*'), path(vcf_index, stageAs: 'dir??/*') //[list(val(sample_id)), list(file(vcf)), list(file(csi))]

    output:
    path('input_vcf_manifest.tsv'), emit: vcf_manifest
    path('user_to_vcf_sampleid_map.tsv'), emit: sample_to_user_map
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), emit: versions_vcfmanifest, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // in this way we rename the directories to match the sample IDs
    def rename_dirs = sample_ids.withIndex()
        .collect { id, i -> "mv dir${String.format('%02d', i + 1)} ${id}" }
        .join('\n')
    """
    ${rename_dirs}

    echo -e "user_sampleid\\tvcf_sampleid\\tfile_date\\tfile_path" > input_vcf_manifest.tsv
    find -maxdepth 2 \\( -name '*.vcf.gz' -o -name '*.bcf' \\) -print \\
    | while IFS= read -r vcf_file
    do
        user_sampleid=\$(basename \$(dirname \$vcf_file))
        vcf_sampleid=\$(bcftools view -h \$vcf_file | tail -1 | cut -f10- | tr '\t' '\n')
        n_samples=\$(echo "\$vcf_sampleid" | wc -l)
        if [ "\$n_samples" -ne 1 ]; then
            echo "ERROR: The VCF file \$vcf_file contains \$n_samples samples. We expect a single sample VCF."
            exit 1
        fi
        file_date=\$(date -r \$vcf_file "+%Y%m%d")
        file_path=\$(readlink -f \$vcf_file)
        echo -e "\$user_sampleid\\t\$vcf_sampleid\\t\$file_date\\t\$file_path" >> input_vcf_manifest.tsv
    done
    awk '{OFS=FS="\t"}; {print \$1,\$2}' input_vcf_manifest.tsv > user_to_vcf_sampleid_map.tsv
    """

    stub:
    def rename_dirs = sample_ids.withIndex()
        .collect { id, i -> "mv dir${String.format('%02d', i + 1)} ${id}" }
        .join('\n')
    """
    ${rename_dirs}

    n=0
    echo -e "user_sampleid\\tvcf_sampleid\\tfile_date\\tfile_path" > input_vcf_manifest.tsv
    find -maxdepth 2 \\( -name '*.vcf.gz' -o -name '*.bcf' \\) -print \\
    | while IFS= read -r vcf_file
    do
        n=\$((n+1))
        user_sampleid=\$(basename \$(dirname \$vcf_file))
        vcf_sampleid=\$(bcftools view -h \$vcf_file | tail -1 | cut -f10- | tr '\t' '\n')
        n_samples=\$(echo "\$vcf_sampleid" | wc -l)
        if [ "\$n_samples" -ne 1 ]; then
            echo "ERROR: The VCF file \$vcf_file contains \$n_samples samples. We expect a single sample VCF."
            exit 1
        fi
        file_date=\$(date -r \$vcf_file "+%Y%m%d")
        file_path=\$(readlink -f \$vcf_file)
        echo -e "\$user_sampleid\\t\$vcf_sampleid\\t\$file_date\\t\$file_path" >> input_vcf_manifest.tsv
    done
    awk '{OFS=FS="\t"}; {print \$1,\$2}' input_vcf_manifest.tsv > user_to_vcf_sampleid_map.tsv
    """
}
