#!/bin/bash
#SBATCH --job-name=WGSQC-metrics
#SBATCH --mail-type=END
#SBATCH --mail-user=edoardo.giacopuzzi@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --output=logs/WGSQC-metrics_%A.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --chdir "/processing_data/shared_datasets/molisani/raw_data/nf-manual-wf"

module load nextflow/21.10.6

#Create a list of new samples to process
ls /processing_data/shared_datasets/molisani/raw_data/wgs/bam/*/*.bam > filelist_bam
grep -v --file processed_bam filelist_bam > bams_to_run

ls /processing_data/shared_datasets/molisani/raw_data/wgs/vcf/*/*.vcf.gz > filelist_vcf
grep -v --file processed_vcf filelist_vcf > vcfs_to_run

#run nextflow here
nextflow run /project/alfredo/pipelines/nf-WGS-QC/main.nf -profile ht_cluster -c wgs_manual_qc.conf

#when completed succesfully add samples to processed ones
if [ $? == 0 ]
then
	cat bams_to_run >> processed_bam
	cat vcfs_to_run >> processed_vcf
else
	exit 1
fi
