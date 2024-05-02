# Python script to create DESeq2 metadata file(s) from the pipeline design file

import argparse
import pandas as pd
import numpy as np
import os


# Setup



def read_options():
    #Arguments and options
    parser = argparse.ArgumentParser(
                prog = 'create_metadata.py',
                description = 'Python script to create DESeq2 metadata file(s) from the pipeline design file',
                epilog = 'Steven Wingett 2024, The MRC-LMB, Cambridge, UK'
    )
    parser.add_argument("--design_file", action="store", type=str, metavar='', help="CSV experiment design file", default='analysis_design_file.csv')
    parser.add_argument("--outdir", action="store", type=str, metavar='', help="Output directory", default='metadata_deseq2')

    args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not
    options = args[0]   # Get the 2 arrays of known/unknown arguments from the tuple

    return options


#Read in data 
options = read_options()


# Read in data
print(f'Reading in design file {options.design_file}')
design = pd.read_csv(options.design_file, header=None)   # Importing this way allows duplicate column names!
design.columns = design.iloc[0]
design = design.drop(0) 


# Format - remove uneeded columns
design = pd.concat([design.iloc[:, 1], design.iloc[:, 3:]], axis=1)


# Make the output directory
output_directory = f'{options.outdir}/metadata_deseq2'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Obtain the column values dividing the metadata groups i.e. all NaN columns
divider_columns = design.isna().all().reset_index().iloc[:, 1]
divider_columns = divider_columns[divider_columns].index.to_list()
divider_columns.append(design.shape[1])

metatdata_file_counter = 1

for i in range(len(divider_columns) - 1):
    
    # Split the design by the columns separated by NA values
    start_column = divider_columns[i] + 1
    end_column = divider_columns[i + 1] - 1
    #print(f"Start:{start_column}   End:{end_column}")
    design_subset = design.iloc[:, start_column:end_column + 1]
    design_subset = pd.concat([design.iloc[:, 0], design_subset], axis=1)
    filt = ~ design_subset.isna().any(axis=1)
    design_subset = design_subset[filt].reset_index(drop=True)
    
    # Identify the control vs sample pairwise comparisons
    deseq_reference = design_subset.columns[1].split('/')[1]
    
    for deseq_sample in design_subset.iloc[:, 1].drop_duplicates():
        if deseq_sample == deseq_reference:   # Don't compare with self!
            continue
            
        filt = design_subset.iloc[:, 1] == deseq_reference
        metadata = design_subset[filt]
        
        filt = design_subset.iloc[:, 1] == deseq_sample
        metadata = pd.concat([metadata, design_subset[filt]], axis=0)
        
        metatdata_column_names = metadata.columns.to_list()
        metatdata_column_names[1] = metatdata_column_names[1].split('/')[0]
        metadata.columns = metatdata_column_names
        metadata = metadata.reset_index(drop=True)
        
        # Print out metadata
        outfile = f'{output_directory}/comparison{metatdata_file_counter}_metadata.tsv'
        print(f'Writing to {outfile}')
        metadata.to_csv(outfile, sep='\t', index=False)
        metatdata_file_counter = metatdata_file_counter + 1

print('Done')

