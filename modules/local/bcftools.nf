/*
Most of these processes expect files and a prefix as input.
One can obtain a prefix for input using something like
datasets = Channel
                .fromPath(params.input)
                .map { file -> tuple(file.baseName, file) }

of to remove a specific extension
file.getName().replaceFirst(/\.vcf\.gz$/, "")
*/

process BCFTOOLS_NORM {
    label 'singlecore'

    input:
        tuple file(vcf_file), file(vcf_index)
        tuple file(genome), file(genome_fai)

    output:
        tuple file("${prefix}.norm.bcf"), file("${prefix}.norm.bcf.csi"), emit: bcf_file
    
    script:
    def prefix = vcf_file.getName().replaceFirst(/(\.bcf|\.vcf\.gz)$/, "")
    """
    bcftools norm -m- -f $genome \
        -Ob -o ${prefix}.norm.bcf $vcf_file

    bcftools index --csi ${prefix}.norm.bcf
    """
}

process BCFTOOLS_FILTER {
    label 'singlecore'
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy'
    }

    input:
        tuple file(vcf_file), file(vcf_index), val(filter_list)
        // filter_list is a list like [[filter_tag, filter_string], [filter_tag, filter_string], ...] 
        // a filter_string is like 'QUAL < 20 || N_PASS(GQ >= 20) == 0'
        val(operation)
        val(output_ext)

    output:
        tuple file("${out_filename}"), file("${out_filename}.csi"), emit: filtered_bcf

    script:
    def filter_logic = operation == 'include' ? "-i" : "-e"
    def output_format_arg = output_ext == 'vcf.gz' ? "z" : "b"
    def filters = filter_list.collate(2,2).collect{ it[0]+"\t"+it[1] }.join("\n")
    prefix = vcf_file.getName().replaceFirst(/(\.bcf|\.vcf\.gz)$/, "")
    out_filename = params.publish_filename ? "${params.publish_filename}.${output_ext}" : "${prefix}.filters.${output_ext}"
"""
function run_filter () {
    cmd="bcftools filter -m + -s \$3 $filter_logic \$4 -Ob -o temp\$2.bcf \$1"
    echo \$cmd
    eval \$cmd
}

cat <<EOT >> filters.txt
$filters
EOT

n=0
inputfile="$vcf_file"
while read -r filter_tag filter_string
do
run_filter \$inputfile \$n \$filter_tag "'\$filter_string'"
inputfile=temp\$n.bcf
n=\$((n + 1))
done < filters.txt

bcftools view -O$output_format_arg -o ${out_filename} \$inputfile 
bcftools index --csi ${out_filename}
"""
}

/*
process BCFTOOLS_FILTER {
    label 'singlecore'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(prefix), file(vcf_file), file(vcf_index)
        val(filters)
        // a list of strings like 'filterTAG,QUAL < 20 || N_PASS(GQ >= 20) == 0'
        val(operation)

    output:
        tuple file("${prefix}.PASS.vcf.gz"), file("${prefix}.PASS.vcf.gz.csi"), emit: pass_vcf
        tuple file("${prefix}.processed.vcf.gz"), file("${prefix}.processed.vcf.gz.csi"), emit: processed_vcf

    script:
    def filter_strings = filters.collect{ "$it" }.join('\n')
    def filter_logic = operation == 'include' ? "-i" : "-e"
    """
    echo '$filter_strings' > filter_definitions.txt
    export n=0
    mv $vcf_file ${prefix}.\$n.vcf.gz
    while read -r filter
    do
        export old=\$n
        export n=\$(( \$n + 1 ))
        filter_name=\$(echo \$filter | cut -d"," -f1)
        filter_string="'"\$(echo \$filter | cut -d"," -f2)"'"
        
        cmd="bcftools filter -m + -s \$filter_name \
        $filter_logic \$filter_string \
        -Oz -o ${prefix}.\$n.vcf.gz \
        ${prefix}.\$old.vcf.gz"

        eval \$cmd  
    done < filter_definitions.txt
    
    mv ${prefix}.\$n.vcf.gz ${prefix}.processed.vcf.gz
    bcftools index --csi ${prefix}.processed.vcf.gz

    bcftools view -f PASS \
    -Oz -o ${prefix}.PASS.vcf.gz \
    ${prefix}.processed.vcf.gz 

    bcftools index --csi ${prefix}.PASS.vcf.gz
    """
}
*/

process BCFTOOLS_MERGE {
    label 'lowcores'
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy'
    }

    input:
        tuple file(vcf_file), file(vcf_index)
        val(out_filename)
        val(merge_mode)

    output:
        tuple file("${out_filename}.bcf"), file("${out_filename}.bcf.csi"), emit: bcf_file
    
    script:
    def input_vcfs = vcf_file.collect{ "$it" }.join(" ")
    """
    bcftools merge --threads ${task.cpus} \
        -F x -m $merge_mode \
        -Ob -o ${out_filename}.bcf \
        $input_vcfs

    bcftools index --csi ${out_filename}.bcf
    """
}

process BCFTOOLS_STATS {
    label 'lowcores'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample), file(vcf_file), file(vcf_index)

    output:
        tuple val(sample), file("${sample}.vars.stats")

    script:
    def make_index = vcf_index.exists() ? '' : "bcftools index $vcf_file"
    """
    $make_index
    bcftools stats -f "PASS,." -s- --threads ${task.cpus} $vcf_file > ${sample}.vars.stats
    """
}

process BCFTOOLS_INDEX {
    label 'lowcores'
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.csi'

    input:
        file(vcf_file)

    output:
        tuple file(vcf_file), file("${vcf_file}.csi")

    script:
    """
    bcftools index -f $vcf_file
    """
}

process BCFTOOLS_IDX_STATS {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample), file(vcf_file), file(vcf_index)

    output:
        tuple val(sample), file("${sample}.varsPerContig.tsv"), emit: vars_stats
        tuple val(sample), file("${sample}.PASS.varsPerContig.tsv"), emit: pass_vars_stats

    script:
    def make_index = vcf_index.exists() ? '' : "bcftools index $vcf_file"
    """
    $make_index
    bcftools view -f "PASS,." -Ob -o PASS_vars.bcf $vcf_file
    bcftools index PASS_vars.bcf
    bcftools index --stats PASS_vars.bcf > ${sample}.PASS.varsPerContig.tsv
    bcftools index --stats $vcf_file > ${sample}.varsPerContig.tsv
    """
}

process BCFTOOLS_FILLTAGS {
    label 'singlethread'
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy'
    }

    input:
        tuple file(vcf_file), file(vcf_index)
        val(tags) //string of comma-separated tags for the plugin or all

    output:
        tuple file("${prefix}.fill.vcf.gz"), file("${prefix}.fill.vcf.gz.csi")

    script:
    def prefix = vcf_file.getName().replaceFirst(/(\.bcf|\.vcf\.gz)$/, "")
    """
    bcftools +fill-tags -Oz -o ${prefix}.fill.vcf.gz $vcf_file -- -t $tags
    bcftools index --csi ${prefix}.fill.vcf.gz
    """
}

process BCFTOOLS_VIEW {
    label 'lowcores'
    if (params.publish) { 
        publishDir "${params.outdir}", mode: 'copy'
    }

    input:
        tuple file(vcf_file), file(vcf_index)
        val(suffix) //string to add to output file
        val(output_ext) //string to set output vcf.gz or bcf
        file(regions_bed) //bed file of regions to extract or file('NO_REGION')
        val(region) //string region to extract or 'NO_REGION'
        val(filter) //string of comma-separated tags to apply to FILTER or 'NO_FILTER'
        val(include) //string to use for include -i or 'NO_INCLUDE'
    
    output:
        tuple val(suffix), file("$out_filename"), file("${out_filename}.csi")

    script:
    def filter_arg = filter != 'NO_FILTER' ? "-f ${filter}" : ""
    def region_string_arg = region != 'NO_REGION' ? "-r ${region}" : ""
    def region_bed_arg = regions_bed.name != 'NO_REGION' ? "-R ${regions_bed}" : ""
    def include_arg = include != 'NO_INCLUDE' ? "-i '${include}'" : ""
    def output_format_arg = output_ext == 'vcf.gz' ? "z" : "b"
    prefix = vcf_file.getName().replaceFirst(/(\.bcf|\.vcf\.gz)$/, "")
    out_filename = params.publish_filename ? "${params.publish_filename}.${output_ext}" : "${prefix}.${suffix}.${output_ext}"
    """
    bcftools view --threads ${task.cpus} \
        -O$output_format_arg -o $out_filename \
        $filter_arg \
        $region_string_arg \
        $region_bed_arg \
        $include_arg \
        $vcf_file 
    bcftools index --csi $out_filename
    """
}

process BCFTOOLS_CONCAT {
    label 'lowcores'
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy'
    }

    input:
        file(vcf_file)
        val(output_ext) //string to set output vcf.gz or bcf
        val(out_filename)
        
    output:
        tuple file("${out_filename}.${output_ext}"), file("${out_filename}.${output_ext}.csi"), emit: bcf_file
    
    script:
    def input_vcfs = vcf_file.collect{ "$it" }.join(" ")
    def output_format_arg = output_ext == 'vcf.gz' ? "z" : "b"
    """
    bcftools concat --threads ${task.cpus} -Ob -o temp.bcf $input_vcfs
    bcftools sort -m 32G -O$output_format_arg -o ${out_filename}.${output_ext} temp.bcf

    bcftools index --csi ${out_filename}.${output_ext}
    """
}

process BCFTOOLS_FIXREF {
    label 'lowcores'

    input:
        tuple file(vcf_file), file(vcf_index)
        tuple file(genome_fa), file(genome_fai)

    output:
        tuple file("${prefix}.fixref.bcf"), file("${prefix}.fixref.bcf.csi"), emit: bcf_file
    
    script:
    def prefix = vcf_file.getName().replaceFirst(/\.bcf$/, "")
    """
    bcftools +fixref --threads ${task.cpus} -Ob -o ${prefix}.fixref.bcf $vcf_file -- -f $genome_fa
    bcftools index --csi ${prefix}.fixref.bcf
    """
}

process BCFTOOLS_ROH {
    label 'singlethread'

    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy'
    }
    
    input:
        tuple file(vcf_file), file(vcf_index), val(sample_id)
        tuple file(af_file), file(af_index)
        file genetic_maps //genetic map files 1 per chromosome named {CHROM}.genetic_map.txt

    output:
        file("${sample_id}_ROH.txt")

    script:
    def af_option = af_file.name != 'NO_FILE' ? "--AF-file ${af_file}" : "--AF-tag AF"
    """
    bcftools roh \
        --threads ${task.cpus} \
        -s $sample_id \
        -m {CHROM}.genetic_map.txt \
        -Or -o ${sample_id}_ROH.txt \
        $af_option $vcf_file
    """
}