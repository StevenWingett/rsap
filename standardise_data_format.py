# Python script to convert input format into a standardised format for the pipeline

import pandas as pd
import argparse
import os
import numpy as np
import qnorm


ANSI_ESC= {'END':'\033[0m', 'RED':'\033[31m'}


def _color(txt, cname='RED'):
    return f"{ANSI_ESC.get(cname, 'RED')}{txt}{ANSI_ESC['END']}"



# Setup
def read_options():
    #Arguments and options
    parser = argparse.ArgumentParser(
                prog = 'standardise_data_format.py',
                description = 'Python script to parse input data matrix into the standard format for this pipeline ',
                epilog = 'Steven Wingett 2024, The MRC-LMB, Cambridge, UK'
    )
    parser.add_argument("--raw_ex", action="store", type=str, metavar='', default='salmon.merged.gene_counts.tsv',
                        help="Path to the raw expression matrix")
    parser.add_argument("--norm_ex", action="store", type=str, metavar='', default='salmon.merged.gene_tpm.tsv',
                        help="Path to the normalised expression matrix (NOT log-transformed)")
    parser.add_argument("--outdir", action="store", type=str, metavar='', help="Output directory", default='expression_data_pipeline_format')
    parser.add_argument("--format", action="store", type=str, metavar='', help="Input data format: nf_core [default], seqmonk", default='nf_core')


    args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not
    options = args[0]   # Get the 2 arrays of known/unknown arguments from the tuple

    return options

options = read_options()

if options.format in (['nf_core', 'seqmonk']):
    print(f'Input format set to {options.format}')
else:
    print(f'Input format "{options.format}" not valid!\nQuiting')
    exit(1)



# Read in data
input_output_file_lookup = {options.raw_ex : 'raw_expression_data_pipeline_format.tsv.gz',
                            options.norm_ex : 'normalised_expression_data_pipeline_format.tsv.gz'
                            }

# Make output folder
outdir = f'{options.outdir}/expression_data_pipeline_format'
if not os.path.exists(outdir):
    os.makedirs(outdir) 


# Process the 2 input files
#raw_ex_gene_ids = pd.Series()
#norm_ex_gene_ids = pd.Series()

for input_file in (options.raw_ex, options.norm_ex):   # options.norm_ex is used later on to create log2 matrix
    print(f'Reading in data file {input_file}')

    
    # Standardise data format
    if options.format == 'seqmonk':
        expression_data = pd.read_csv(f'{input_file}',         
                            sep='\t', dtype = {'Chromosome': str}    # Import and make sure chromosome are imported as dtypes to avoid warnings
                        )
        expression_data.iloc[:, 0] = expression_data.iloc[:, 1] + '..' + expression_data.iloc[:, 0]   # Add chromosome name
        columns_to_select = [0, 5] + list(range(12, expression_data.shape[1]))
        expression_data = expression_data.iloc[:, columns_to_select ]
        column_names = expression_data.columns.to_list()
        column_names[0] = 'gene_id'
        column_names[1] = 'gene_name'
        expression_data.columns = column_names
    else:
        expression_data = pd.read_csv(f'{input_file}',         
                            sep='\t'
                        )

    outfile = f'{outdir}/{input_output_file_lookup[input_file]}'
    print(f'Writing to {outfile}')
    expression_data.to_csv(outfile, sep='\t', index=False)

    if 'raw_ex_gene_ids' not in globals():   # For checking input
        raw_ex_gene_ids = expression_data['gene_id'].copy()
    else:
        norm_ex_gene_ids = expression_data['gene_id'].copy()  


# Let's check the raw/normalised data correspond
if len(raw_ex_gene_ids) != len(norm_ex_gene_ids):
    print(_color(f'Problem: {len(raw_ex_gene_ids)} genes ids raw data, but {len(norm_ex_gene_ids)} in normalised data!'))
    exit(1)

# Check raw/normalised gene ids match exactly
filt = raw_ex_gene_ids != norm_ex_gene_ids
filt = filt.any()
if filt:
    print(_color(f'Problem: genes ids raw data do not all match normalised data!'))
    exit(1)

# Make the log2-transformed normalised matrix
expression_data.iloc[:, 2:] = np.log2(expression_data.iloc[:, 2:] + 1)
outfile = f'{outdir}/log2_normalised_expression_data__plus_1_pipeline_format.tsv.gz'
print(f'Writing to {outfile}')
expression_data.to_csv(outfile, sep='\t', index=False)

# Create a gene name / gene id lookup table
lookup_table_data = (expression_data
                     .loc[:, ['gene_id', 'gene_name']]
                     .copy()
                     .drop_duplicates()
)

print(f'{lookup_table_data.shape[0]} unique gene_id : gene_name entries in lookup table')
outfile = f'{outdir}/gene_id_gene_name_lookup_table.tsv.gz'
print(f'Writing to {outfile}')
lookup_table_data.to_csv(outfile, sep='\t', index=False)


# Make a quantile-normalised matrix
quantile_normalised_expression_data = expression_data.iloc[:, 2:].copy()
quantile_normalised_expression_data = qnorm.quantile_normalize(quantile_normalised_expression_data, axis=1)
quantile_normalised_expression_data = pd.concat([expression_data.iloc[:, :2], quantile_normalised_expression_data], axis=1)
outfile = f'{outdir}/quantile_normalised_log2_normalised_expression_data__plus_1_pipeline_format.tsv.gz'
print(f'Writing to {outfile}')
quantile_normalised_expression_data.to_csv(outfile, sep='\t', index=False)

print('Done')