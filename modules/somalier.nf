// EXTRACT
// Create somalier profiles
/*
INPUT
sample name, bam file, bai file
var sites file: var sites to extract as .vcf.gz
genome_fasta, genome_fai_idx

OUTPUT
${sample}.somalier

PUBLISH
${sample}.somalier

PARAMS
outdir: publish dir
*/

process SOMALIER_EXTRACT {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.somalier"

    input:
    tuple file(bam_file), file(bai_file)
    file(var_sites)
    tuple file(genome_fasta), file(genome_index)

    output:
    file("*.somalier")

    script:
    """  
    somalier extract -d ./ -f $genome_fasta -s $var_sites $bam_file
    """
}

// RELATE
// Estimate sex and relationship across samples using profiles
/*
INPUT
series of somalier profile files (.somalier)
ped file

OUTPUT
cohort_relate.*

PUBLISH
cohort_relate.*

PARAMS
outdir: publish dir
*/
process SOMALIER_RELATE {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.reportdir}", mode: 'copy', pattern: '*.html'
    
    input:
    file(somalier)
    file(ped_file)

    output:
    file 'somalier-relate.*'

    script:
    """
    somalier relate -p $ped_file -o \$TMPDIR/somalier-relate *.somalier
    cp \$TMPDIR/somalier-relate* ./
    """
}

// ANCESTRY
// Estimate ancestry components from profiles based on 1000G data
/*
INPUT
series of somalier profile files (.somalier)
tsv file with 1000g ancestry labels
folder containing somalier profiles for 1000G samples (.somalier)

OUTPUT
cohort_ancestry.*

PUBLISH
cohort_ancestry.*

PARAMS
outdir: publish dir
*/
process SOMALIER_ANCESTRY {
    label 'singlethread'
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.reportdir}", mode: 'copy', pattern: '*.html'
    
    input:
    file(somalier)
    file(ancestry_labels)
    path(data_folder)

    output:
    path 'somalier.somalier-ancestry.html', emit: html_file
    path 'somalier.somalier-ancestry.tsv', emit: ancestry_tsv

    script:
    """
    somalier ancestry -o \$TMPDIR/ --labels $ancestry_labels ${data_folder}/*.somalier ++ *.somalier
    cp \$TMPDIR/somalier.somalier-ancestry* ./
    """
}