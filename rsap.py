# (RNA-seq analysis pipeline takes quantitated RNA-seq data (both raw counts and normalised) and performs QC, DE-seq2, GO-terms analysis

import argparse
import os
import sys
import glob
import pandas as pd


# Setup
RSAP_VERSION = '0.0.1'

ANSI_ESC= {'END':'\033[0m', 'BOLD':'\033[1m',
           'ITALIC':'\033[2m', 'UNDERLINE':'\033[3m',
           'BLACK':'\033[30m', 'RED':'\033[31m',
           'GREEN':'\033[32m', 'YELLOW':'\033[33m',
           'BLUE':'\033[34m', 'MAGENTA':'\033[35m',
           'CYAN':'\033[36m', 'WHITE':'\033[37m',
           'GREY':'\033[90m', 'LT_RED':'\033[91m',
           'LT_GREEN':'\033[92m',  'LT_YELLOW':'\033[93m',
           'LT_BLUE':'\033[94m', 'LT_MAGENTA':'\033[95m',
           'LT_CYAN':'\033[96m', 'LT_WHITE':'\033[97m'}


def _color(txt, cname='BLUE'):
  
    return f"{ANSI_ESC.get(cname, 'BLUE')}{txt}{ANSI_ESC['END']}"



def read_options():
    #Arguments and options
    parser = argparse.ArgumentParser(
                prog = 'rsap.py',
                description = '''\
                    RNA-seq analysis pipeline takes quantitated RNA-seq data (both raw counts and normalised) and performs QC, DE-seq2, GO-terms analysis.
                
                    Note: a container has been made available to run with this script: https://hub.docker.com/repository/docker/swingett/rsap_deseq2/general
                    ''',
                epilog = 'Steven Wingett 2024, The MRC-LMB, Cambridge, UK'
    )

    parser.add_argument("--raw_ex", action="store", type=str, metavar='', default='salmon.merged.gene.tsv',
                        help="Path to the raw count expression matrix")
    parser.add_argument("--norm_ex", action="store", type=str, metavar='', default='salmon.merged.gene_tpm.tsv',
                        help="Path to the normalised expression matrix")  
    parser.add_argument("--design_file", action="store", type=str, metavar='', help="CSV experiment design file", default='analysis_design_file.csv')
    parser.add_argument("--outdir", action="store", type=str, metavar='', help="Output directory", default='results_rsap')
    parser.add_argument("--padj_threshold", action="store", type=float, metavar='', default=0.05, 
                        help="DESeq2 maxmium adjusted p-value threshold"
                        )
    parser.add_argument("--abs_l2fc_threshold", action="store", type=float, metavar='', default=0.584, 
                        help="Minimum DESeq2 absolute log2-fold change threshold"
                        )
    parser.add_argument("--format", action="store", type=str, metavar='', help="Input data format: nf_core [default], seqmonk", default='nf_core')
    parser.add_argument("--skip_go", action='store_true', metavar='', help="Skip the Human GO terms Enrichment Analysis")



    parser.add_argument("--version", action='store_true', help="Print RSAP version and quit")


    args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not
    options = args[0]   # Get the 2 arrays of known/unknown arguments from the tuple

    return options


#Read in data 
options = read_options()

if options.version == True:
    print(f'RSAP v{RSAP_VERSION}')
    exit(0)


def main():
    rsap_folder = os.path.abspath(os.path.dirname(sys.argv[0]))

    # Get options
    options = read_options()

    if options.format in (['nf_core', 'seqmonk']):
        print(f'Input format set to {options.format}')
    else:
        print(f'Input format "{options.format}" not valid!\nQuiting')
        exit(1)


    # Standardise data format
    print(_color('******CONVERT DATA TO PIPELINE FORMAT******'))
    command = f'python3 {rsap_folder}/standardise_data_format.py --raw_ex {options.raw_ex} --norm_ex {options.norm_ex} --outdir {options.outdir} --format {options.format}'
    print(f'Running command:\n{command}')
    os.system(command)

    # Output files produced, (to be used in this pipeline):
    raw_expression_file = f'{options.outdir}/expression_data_pipeline_format/raw_expression_data_pipeline_format.tsv.gz'
    normalised_expression_file = f'{options.outdir}/expression_data_pipeline_format/normalised_expression_data_pipeline_format.tsv.gz'
    log2_normalised_expression_file = f'{options.outdir}/expression_data_pipeline_format/log2_normalised_expression_data__plus_1_pipeline_format.tsv.gz' 


    # Perform QC
    print(_color('******PERFORM QC******'))
    print(f'Running command:\n{command}')
 
    command = f'python3 {rsap_folder}/qc.py --expression_file {normalised_expression_file} --design_file {options.design_file} --outdir {options.outdir}'
    os.system(command)


    # Make DESeq2 Metadata files
    print(_color('******MAKE DESEQ2 METADATA FILES******'))
    command = f'python3 {rsap_folder}/create_metadata.py --design_file {options.design_file} --outdir {options.outdir}'
    print(f'Running command:\n{command}')
    os.system(command)


    # Perform DESeq2 analysis
    metadata_files = glob.glob(f'{options.outdir}/metadata_deseq2/*.tsv')

    for metadata_file in metadata_files:
        #print(metadata_file)
        comparison_description = os.path.basename(metadata_file)
        comparison_description = comparison_description.split('_')[0] + '__'

        # Determine comparison name
        metadata = pd.read_csv(metadata_file, sep='\t')
        comparison_description = comparison_description + '_'.join(metadata.columns[1:].to_list())   # Variable to test and contrast(s)
        comparison_description = comparison_description + '__' + '_vs_'.join(metadata.iloc[:, 1].drop_duplicates().to_list())
        outdir = f'{options.outdir}/deseq2/{comparison_description}'

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        command = f'Rscript {rsap_folder}/run_deseq2.R {raw_expression_file} {metadata_file} {outdir} > {outdir}/deseq.log'
        print(f'Running command:\n{command}')
        os.system(command)


    # Summarise DEseq results
    print(_color('******SUMMARISE DESEQ2 RESULTS******'))

    
     # Use the metadata to detmine the location of the DESeq2 output
    for metadata_file in metadata_files:
        #print(metadata_file)
        metadata = pd.read_csv(metadata_file, sep='\t')

        outdir = f'{options.outdir}/deseq2/'
        outdir = outdir + os.path.basename(metadata_file)[:-13]  # Get comparison id
        outdir = outdir + '__' + ('_'.join(metadata.columns.to_list()[1:]))

        groups = metadata.iloc[:, 1].drop_duplicates().to_list()
        outdir = outdir + '__' + groups[0] + '_vs_' + groups[1]

        deseq_file = glob.glob(f'{outdir}/*/*.deseq2_results.tsv')

        if len(deseq_file) == 1:
            deseq_file = deseq_file[0]
        else:
            print('Mutliple DESeq2 files returned - this should not happen!!!!')
            print(deseq_file)
            exit(1)

        #print(deseq_file)
        outdir = outdir + '/comparison_summary'
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        print(f'Writing DESeq2 summary results to {outdir}')

        command = f'python3 {rsap_folder}/summarising_deseq_results.py --metadata_file {metadata_file} --deseq_file {deseq_file} --log2_norm_ex_file {log2_normalised_expression_file} --padj_threshold {options.padj_threshold} --abs_l2fc_threshold {options.abs_l2fc_threshold} --outdir {outdir}'
        print(f'Running command:\n{command}')
        os.system(command)


    # Perform GO Analysis
    if options.skip_go:
        print(_color('******GO TERM ENRICHMENT ANALYSIS - SKIPPING******'))
    else:
        print(_color('******GO TERM ENRICHMENT ANALYSIS******'))
    
    
        # Use the metadata to detmine the location of the DESeq2 output
        for metadata_file in metadata_files:
            #print(metadata_file)
            metadata = pd.read_csv(metadata_file, sep='\t')

            outdir = f'{options.outdir}/deseq2/'
            outdir = outdir + os.path.basename(metadata_file)[:-13]  # Get comparison id
            outdir = outdir + '__' + ('_'.join(metadata.columns.to_list()[1:]))

            groups = metadata.iloc[:, 1].drop_duplicates().to_list()
            outdir = outdir + '__' + groups[0] + '_vs_' + groups[1]
            outdir = outdir + '/comparison_summary'

            for deg_file in ('down_regulated_genes.tsv', 'up_regulated_genes.tsv'):
                deg_file = f'{outdir}/{deg_file}'   

                if not os.path.exists(deg_file):
                    print(f'Could not find {deg_file}')
                    print('This should not happen!!!!')
                    exit(1)

                command = f'python3 {rsap_folder}/go_terms_panther.py --gene_ids_file {deg_file} --outdir {outdir}/go_enrichment'
                print(f'Running command:\n{command}')
                os.system(command)


    print(_color('******PIPELINE COMPLETE******'))



if __name__ == "__main__":
    main()