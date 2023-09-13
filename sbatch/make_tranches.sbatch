#!/bin/bash
#SBATCH --job-name=WGSQC-maketranche
#SBATCH --mail-type=END
#SBATCH --mail-user=edoardo.giacopuzzi@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --output=logs/WGSQC-maketranche_%A.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --chdir "/processing_data/shared_datasets/molisani/raw_data/WGS_QC_data"

#Activate conda env
echo "Activating conda env"
eval "$(conda shell.bash hook)"
conda activate pythonenv_latest

echo "Running make_tranches_time.py"
/ssu/gassu/pipelines/nf-WGS-QC/bin/make_tranches_time.py --manifest_file bam_manifest.tsv --start_date 20230808 --force_update
