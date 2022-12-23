#!/usr/bin/env python

import os, shutil
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
parser.add_argument('--force_update', action="store_true", help="Overwrite tranches that already exist")

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

# Symlink files in each interval and create counts report
daily_counts_list = []
for start, end in date_intervals:
    tranche_name = "tranche_{}_{}".format(start.strftime("%Y%m%d"), end.strftime("%Y%m%d"))
    tranche_dir = os.path.join(dest_dir,tranche_name)
    print(f"tranche dir: {tranche_dir}")
    if os.path.exists(tranche_dir):
        if args.force_update:
            shutil.rmtree(tranche_dir)
        else:
            print(f"{tranche_dir} already present. Skip this one")
            continue
    
    
    pd_date_range = pd.date_range(start,end, inclusive='left')
    samples_in_interval = df[df['date'].isin(pd_date_range)]
    print(f"{len(samples_in_interval)} samples in this interval")

    samples_in_interval.apply(lambda row: symlink_files(row, source_dir, tranche_dir), axis=1)
    
    # Count the number of rows for each day in the date range
    counts = df[df['date'].isin(pd_date_range)].groupby('date').size().reset_index(name='count')

    # Merge the counts with the date range to include zero values for days with no rows
    daily_counts = pd.merge(pd_date_range.to_frame(name='date'), counts, on='date', how='left').fillna(0)
    daily_counts['count'] = daily_counts['count'].astype('int')

    # Add tranche id and save
    daily_counts['tranche'] = tranche_name
    daily_counts_list.append(daily_counts)

if (len(daily_counts_list) > 0):
    daily_counts_report = pd.concat(daily_counts_list)
    daily_counts_report.to_csv("daily_counts_report.tsv", sep="\t", index=False)