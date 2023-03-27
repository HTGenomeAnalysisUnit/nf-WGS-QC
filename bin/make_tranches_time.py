#!/usr/bin/env python

print ("### - Make tranches by time - ###")
print("### - Importing libraries - ###")

import os, shutil
from datetime import datetime, timedelta
import pandas as pd
from glob import glob
import re
import argparse

print("### - Program starts - ###")

re_exp_align = re.compile('.+/qc_alignment/')
re_exp_vars = re.compile('.+/qc_variants/')

def symlink_files(row, source_dir, dest_dir, force_update):
    source_abs_path = os.path.abspath(source_dir)
    print(f"### - Sample: {row['external_id']}")
    align_qc_files = glob(f"{source_abs_path}/qc_alignment/**/{row['external_id']}.*", recursive=True)
    vars_qc_file = glob(f"{source_abs_path}/qc_variants/**/{row['external_id']}.*", recursive=True)
    files_in_tranche = align_qc_files + vars_qc_file
    
    for source_file in files_in_tranche:
        dirname = os.path.dirname(source_file)
        subdir = re.sub(re_exp_align,'qc_alignment/',dirname)
        subdir = re.sub(re_exp_vars,'qc_variants/',subdir)
        target_dir = os.path.join(dest_dir,subdir)
        basename = os.path.basename(source_file)
        target_file = os.path.join(target_dir, basename)
        #print(f"file: {source_file}, subdir: {subdir}, target_dir: {target_dir}")
        os.makedirs(target_dir,exist_ok=True)
        #check if file already exists and remove it eventually
        if os.path.exists(target_file):
            if force_update:
                os.remove(target_file)
            else:
                continue
        if not os.path.exists(target_file):
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

print(f"Reading manifest file: {manifest_file}")
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

# Filter df to only include samples with the most recent date
df_most_recent = df.sort_values('date', ascending=False).drop_duplicates('external_id')

# Symlink files in each interval and create counts report
daily_counts_list = []
for start, end in date_intervals:
    tranche_name = "tranche_{}_{}".format(start.strftime("%Y%m%d"), end.strftime("%Y%m%d"))
    tranche_dir = os.path.join(dest_dir,tranche_name)
    print(f"tranche dir: {tranche_dir}")
    if os.path.exists(tranche_dir):
        if args.force_update:
            print(f"{tranche_dir} already present. Force update this one")
            shutil.rmtree(tranche_dir)
        else:
            print(f"{tranche_dir} already present. Skip this one")
            continue
    
    pd_date_range = pd.date_range(start,end, inclusive='left')
    # For this we use df_most_recent to only include the most recent sample for each external_id
    samples_in_interval = df_most_recent[df_most_recent['date'].isin(pd_date_range)]
    print(f"{len(samples_in_interval)} samples in this interval")

    samples_in_interval.apply(lambda row: symlink_files(row, source_dir, tranche_dir, args.force_update), axis=1)
    
    # Count the number of rows for each day in the date range
    # This is done on the original df to include all samples
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