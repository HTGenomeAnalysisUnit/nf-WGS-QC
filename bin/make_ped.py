#!/usr/bin/env python

import sys
import pandas as pd
import argparse

VERSION = 1.0

# Handle --version early so it works without required arguments
version_parser = argparse.ArgumentParser(add_help=False)
version_parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
version_parser.parse_known_args()

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Create PED file from a manifest file")
parser.add_argument("manifest_file", help="Path to the project manifest files. This is a tab-separated file containing sample information")
parser.add_argument("--out", help="Path to the output PED file", required=True)
parser.add_argument("--sample_id_col", help="Name of the column containing sample IDs in the manifest file", default="bam_sample_id")
parser.add_argument("--sex_col", help="Name of the column containing sex information in the manifest file", default="sex")
parser.add_argument("--sex_encode", help="Encoding for sex column.", choices=["01", "12", "MF"], default="12")
args = parser.parse_args()

manifest_file = args.manifest_file
sample_id_col = args.sample_id_col
sex_col = args.sex_col
sex_encode = args.sex_encode
output_file = args.out

# Read the manifest file into a DataFrame
manifest_df = pd.read_csv(manifest_file, sep='\t')

# Check sample ID and sex columns exist in the manifest file, otherwise exit with an error
if sample_id_col not in manifest_df.columns:
    sys.exit(f"Error: Sample ID column '{sample_id_col}' not found in manifest file")
if sex_col not in manifest_df.columns:
    sys.exit(f"Error: Sex column '{sex_col}' not found in manifest file")

# Check that sex column contains valid values based on the specified encoding, otherwise exit with an error
if sex_encode == "01":
    valid_values = {0, 1, -9}
elif sex_encode == "12":
    valid_values = {0, 1, 2, -9}
elif sex_encode == "MF":
    valid_values = {"M", "F"}

invalid_values = set(manifest_df[sex_col].dropna()) - valid_values
if invalid_values:
    print(f"Invalid values found in sex column: {invalid_values}")
    print(f"Valid values for sex column based on encoding '{sex_encode}': {valid_values}")
    sys.exit(f"Error: Sex column contains invalid values: {invalid_values}")

# Encode sex column based on the specified encoding
# Somalier wants sex encoded as 1 for male and 2 for female, with 0 for unknown
if sex_encode == "01":
    manifest_df[sex_col].replace(0, 2, inplace=True)
elif sex_encode == "MF":
    manifest_df[sex_col].replace({"M": 1, "F": 2}, inplace=True)

# This ensures missing values are set to 0 (unknown) as expected by Somalier
manifest_df[sex_col].replace(-9, 0, inplace=True)
manifest_df.fillna({sex_col: 0}, inplace=True)

# Create a dummy ped file with sex information
new_df = pd.DataFrame({
    'FID': manifest_df[sample_id_col],
    'IID': manifest_df[sample_id_col],
    'father_id': 0,
    'mother_id': 0,
    'sex': manifest_df[sex_col],
    'status': 0
})

# Save to PED file
new_df.to_csv(output_file, sep='\t', header=False, index=False)
