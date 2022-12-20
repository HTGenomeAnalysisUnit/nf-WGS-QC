#!/usr/bin/env python

import os, sys, shutil
from datetime import datetime, timedelta
import pandas as pd
from glob import glob
import re
import argparse

re_exp_align = re.compile('.+/qc_alignment/')
re_exp_vars = re.compile('.+/qc_variants/')

def symlink_files(row, source_dir, dest_dir):
    source_abs_path = os.path.abspath(source_dir)
    files_in_tranche = glob(f"{source_abs_path}/**/{row['external_id']}*", recursive=True)
    for source_file in files_in_tranche:
        dirname = os.path.dirname(source_file)
        subdir = re.sub(re_exp_align,'qc_alignment/',dirname)
        subdir = re.sub(re_exp_vars,'qc_variants/',subdir)
        target_dir = os.path.join(dest_dir,subdir)
        basename = os.path.basename(source_file)
        target_file = os.path.join(target_dir, basename)
        #print(f"file: {source_file}, subdir: {subdir}, target_dir: {target_dir}")
        os.makedirs(target_dir,exist_ok=True)

        os.symlink(source_file, target_file)

# Command line parser
parser = argparse.ArgumentParser()

# Add the arguments
parser.add_argument('--source_dir', default='./', help='The source directory')
parser.add_argument('--dest_dir', default='by_tranche_data', help='The destination directory')
parser.add_argument('--start_date', required=True, help='The start date')
parser.add_argument('--days_interval', default=7, type=int, help='The number of days in each interval')
parser.add_argument('--manifest_file', required=True, help='The manifest file')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
source_dir = args.source_dir
dest_dir = args.dest_dir
start_date = args.start_date
days_interval = args.days_interval
manifest_file = args.manifest_file


df = pd.read_csv(manifest_file, sep="\t", names=['internal_id', 'external_id', 'date', 'file_path'], parse_dates=["date"])
print(f"{len(df)} files found in {manifest_file}")

range_start = datetime.strptime(start_date, "%Y%m%d")
current_date = datetime.now()

# Create a list of date intervals
date_intervals = []
while range_start < current_date:
    range_end = range_start + timedelta(days=7)
    date_intervals.append((range_start, range_end))
    range_start = range_end

for start, end in date_intervals:
    tranche_dir = "tranche_{}_{}".format(start.strftime("%Y%m%d"), end.strftime("%Y%m%d"))
    tranche_dir = os.path.join(dest_dir,tranche_dir)
    print(f"tranche dir: {tranche_dir}")
    
    samples_in_interval = df[(df["date"] >= start.strftime("%Y%m%d")) & (df["date"] < end.strftime("%Y%m%d"))]
    print(f"{len(samples_in_interval)} samples in this interval")

    samples_in_interval.apply(lambda row: symlink_files(row, source_dir, tranche_dir), axis=1)