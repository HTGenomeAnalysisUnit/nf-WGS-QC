#!/bin/bash
#SBATCH --job-name=WGSQC-metrics
#SBATCH --mail-type=END
#SBATCH --mail-user=edoardo.giacopuzzi@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --output=logs/WGSQC-metrics_%A.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --chdir "/processing_data/shared_datasets/molisani/raw_data/nf-manual-wf"

module load nextflow/22.10.1

#Create a list of new samples to process
#ls /processing_data/shared_datasets/molisani/raw_data/wgs/bam/*/*.bam > filelist_bam
echo "Getting new BAMs"
find /processing_data/shared_datasets/molisani/raw_data/wgs_genomics_facility/data_by_sample -name "*.bam" | while read -r filename; do echo -e "$(date -r $filename '+%Y%m%d')\t$filename"; done > filelist_bam
grep -v --file processed_bam_molisaniGF filelist_bam > bams_to_run

#ls /processing_data/shared_datasets/molisani/raw_data/wgs/vcf/*/*.vcf.gz > filelist_vcf
echo "Getting new VCFs"
#find /processing_data/shared_datasets/molisani/raw_data/wgs/vcf -name "*.vcf.gz" | while read -r filename; do echo -e "$(date -r $filename '+%Y%m%d')\t$filename"; done > filelist_vcf
#grep -v --file processed_vcf filelist_vcf > vcfs_to_run

#run nextflow here
echo "Running nextflow"
nextflow run /ssu/gassu/pipelines/nf-WGS-QC/main.nf -profile ht_cluster,cpuq -c wgs_manual_qc.conf

#when completed succesfully add samples to processed ones
if [ $? == 0 ]
then
	cat bams_to_run >> processed_bam_molisaniGF
	cat vcfs_to_run >> processed_vcf
else
	exit 1
fi
