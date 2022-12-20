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

module load nextflow/21.10.6

#run nextflow here
nextflow run /project/alfredo/pipelines/nf-WGS-QC/main.nf -profile ht_cluster -entry MAKE_REPORTS -c make_reports.conf