process BAM_MANIFEST {
    tag "${params.project_id}_bam_manifest"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'

    input:
    tuple val(sample_ids), path(bam, stageAs: 'dir??/*'), path(bai, stageAs: 'dir??/*') //[list(val(sample_id)), list(file(bam)), list(file(bai))]

    output:
    path('input_bam_manifest.tsv'), emit: bam_manifest
    path('user_to_bam_sampleid_map.tsv'), emit: sample_to_user_map
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_bammanifest, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // in this way we rename the directories to match the sample IDs
    def rename_dirs = sample_ids.withIndex()
        .collect { id, i -> "mv dir${String.format('%02d', i + 1)} ${id}" }
        .join('\n')

    """
    ${rename_dirs}

    echo -e "user_sampleid\\tbam_sampleid\\tfile_date\\tfile_path" > input_bam_manifest.tsv
    find -maxdepth 2 \\( -name '*.bam' -o -name '*.cram' \\) -print \\
    | while IFS= read -r bam_file
    do
        user_sampleid=\$(basename \$(dirname \$bam_file))
        bam_sampleid=\$(samtools view -H \$bam_file | grep "^@RG" | tr "\\t" "\\n" | grep "SM:" | sed 's/SM://' | sort -u)
        file_date=\$(date -r \$bam_file "+%Y%m%d")
        file_path=\$(readlink -f \$bam_file)
        echo -e "\$user_sampleid\\t\$bam_sampleid\\t\$file_date\\t\$file_path" >> input_bam_manifest.tsv
    done
    awk '{OFS=FS="\t"}; {print \$1,\$2}' input_bam_manifest.tsv > user_to_bam_sampleid_map.tsv
    """

    stub:
    def rename_dirs = sample_ids.withIndex()
        .collect { id, i -> "mv dir${String.format('%02d', i + 1)} ${id}" }
        .join('\n')
    """
    ${rename_dirs}

    n=0
    echo -e "user_sampleid\\tbam_sampleid\\tfile_date\\tfile_path" > input_bam_manifest.tsv
    find -maxdepth 2 \\( -name '*.bam' -o -name '*.cram' \\) -print \\
    | while IFS= read -r bam_file
    do
        n=\$((n+1))
        user_sampleid=\$(basename \$(dirname \$bam_file))
        bam_sampleid=\$(samtools view -H \$bam_file | grep "^@RG" | tr "\\t" "\\n" | grep "SM:" | sed 's/SM://' | sort -u)
        file_date=\$(date -r \$bam_file "+%Y%m%d")
        file_path=\$(readlink -f \$bam_file)
        echo -e "\$user_sampleid\\t\$bam_sampleid\\t\$file_date\\t\$file_path" >> input_bam_manifest.tsv
    done
    awk '{OFS=FS="\t"}; {print \$1,\$2}' input_bam_manifest.tsv > user_to_bam_sampleid_map.tsv
    """
}
