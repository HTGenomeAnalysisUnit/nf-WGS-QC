#!/usr/bin/bash

source_dir="/project/alfredo/pipelines/nf-WGS-QC/sbatch/"

jid1=$(sbatch $source_dir/process_new_files.sbatch)
jid2=$(sbatch --dependency=afterok:${jid1##* } $source_dir/make_tranches.sbatch)
sbatch --dependency=afterok:${jid2##* } $source_dir/make_reports.sbatch