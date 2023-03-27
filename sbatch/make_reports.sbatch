#!/bin/bash
#SBATCH --job-name=WGSQC-reports
#SBATCH --mail-type=END
#SBATCH --mail-user=edoardo.giacopuzzi@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --output=logs/WGSQC-reports_%A.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --chdir "/processing_data/shared_datasets/molisani/raw_data/nf-report-wf"

QC_folder="/processing_data/shared_datasets/molisani/raw_data/WGS_QC_data"

# if qc_align_summary.tsv exists rename it to include its date
if [ -f ${QC_folder}/qc_align_summary.tsv ]
then
	echo "Renaming qc_align_summary.tsv to qc_align_summary_$(date +%Y%m%d).tsv"
	mv ${QC_folder}/qc_align_summary.tsv ${QC_folder}/qc_align_summary_$(date +%Y%m%d).tsv
fi

# if qc_vars_summary.tsv exists rename it to include its date
if [ -f ${QC_folder}/qc_vars_summary.tsv ]
then
	echo "Renaming qc_vars_summary.tsv to qc_vars_summary_$(date +%Y%m%d).tsv"
	mv ${QC_folder}/qc_vars_summary.tsv ${QC_folder}/qc_vars_summary_$(date +%Y%m%d).tsv
fi

module load nextflow/21.10.6

#run nextflow here
nextflow run /project/alfredo/pipelines/nf-WGS-QC/main.nf -profile ht_cluster -entry MAKE_REPORTS -c make_reports.conf