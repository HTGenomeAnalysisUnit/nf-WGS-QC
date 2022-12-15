#!/usr/bin/env python

import sys
import pandas as pd

# Read the first command-line argument as the name of the file to read
project_file = sys.argv[1]

# Read the second command-line argument as the name of the file to use for the dictionary
sid_map_file = sys.argv[2]

# Read the project metadata file into a Pandas DataFrame
project_df = pd.read_csv(project_file, sep='\t')

# Read the map_table file into a Pandas DataFrame
sid_map_table_df = pd.read_csv(sid_map_file, sep='\t', names = ['internal_id', 'external_id'])

# Create a dictionary from the map_table DataFrame
# with keys from the "internal_id" column and values from the "external_id" column
sid_map = sid_map_table_df.set_index('internal_id')['external_id'].to_dict()

# If any sample from the dictionary is not found in the metadata raise an error
if not project_df['WGS_id'].isin(sid_map.keys()).all():
    raise ValueError("Some of the sample IDs from BAM/VCF were not found in the WGS_id column in metadata")

# Use the dictionary to convert the values in the "Sample_ID" column of project_df
# to the new external ids from the dictionary
# retain only samples present in the dictionary
project_df = project_df[project_df['WGS_id'].isin(sid_map.keys())]
project_df['WGS_id'] = project_df['WGS_id'].map(sid_map)

# Replace all 0 values in the "sex" column with 2
project_df['Sex'].replace(0, 2, inplace=True)

# Create a dummy ped file with sex information
new_df = pd.DataFrame({
    'FID': project_df['WGS_id'],
    'IID': project_df['WGS_id'],
    'father_id': 0,
    'mother_id': 0,
    'sex': project_df['Sex'],
    'status': 0
})

# Save the new DataFrame to a tab-separated file named "ped_file.ped"
# without header or row names
new_df.to_csv('ped_file.ped', sep='\t', header=False, index=False)
