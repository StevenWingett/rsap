import os
import json
import argparse
import pandas as pd

# curl -X POST "https://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=Q96PB1&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR" -H "accept: application/json"

# Specify input file

# Setup
def read_options():
    #Arguments and options
    parser = argparse.ArgumentParser(
                prog = 'go_terms_panther.py',
                description = 'Python script to summarise perform GO-enrichment analysis',
                epilog = 'Steven Wingett 2024, The MRC-LMB, Cambridge, UK'
    )
    parser.add_argument("--gene_ids_file", action="store", type=str, metavar='', 
                        help="Path to the file containing gene ids to assess")
    
    parser.add_argument("--outdir", action="store", type=str, metavar='', help="Output directory", default='./go_enrichment')

    args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not
    options = args[0]   # Get the 2 arrays of known/unknown arguments from the tuple

    return options

options = read_options()


fdr_threshold = 0.05

# Input
#geneInputList = 'Q96PB1'
organism = '9606'   # Only works with human at the moment
refInputList  = ''
refOrganism = ''  


print('GO enrichment currently only works for HUMAN!!!!')

# Read in data
#Read in data 
print('Reading in gene list: ' + options.gene_ids_file)
geneInputList = pd.read_csv(options.gene_ids_file, sep='\t')

print(f'{geneInputList.shape[0]} genes to assess')
geneInputList = geneInputList.iloc[:, 0].to_list()
geneInputList = ','.join(geneInputList)


# Make the output directory
if not os.path.exists(options.outdir):
    os.makedirs(options.outdir) 

json_outfile = f'{options.outdir}/{os.path.basename(options.gene_ids_file)}.pantherdb.json'
tsv_outfile = f'{options.outdir}/{os.path.basename(options.gene_ids_file)}.pantherdb.tsv'

# Make the command
curl_command = 'curl -X POST "https://pantherdb.org/services/oai/pantherdb/enrich/overrep?'

geneInputList = f'geneInputList={geneInputList}'
organism = f'organism={organism}'

if refInputList != '':
    refInputList = f'refInputList={refInputList}'
    refOrganism = f'refOrganism={refOrganism}'    # Set so target and reference gene lists are from the same organism
    command_parameters = '&'.join((geneInputList, organism, refInputList, refOrganism, ''))  # Add a final & at the end
else:
    command_parameters = '&'.join((geneInputList, organism, ''))  # Add a final & at the end

# Sets other parameters, including search for Molecular Function GO terms (GO:0008150)
remaining_command = 'annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR" -H "accept: application/json"'
redirect = f' > {json_outfile}'

#print(command_parameters)

command = curl_command + command_parameters + remaining_command + redirect

print(command)

#exit()
os.system(command)


# JSON file
with open(json_outfile, 'r') as f, open(tsv_outfile, 'w') as f_out:

    # Reading from file
    data = json.loads(f.read())

    # Iterating through the json list
    header = '\t'.join(['id', 'label', 'number_in_list', 'fold_enrichment', 'fdr'])
    f_out.writelines(header + '\n')

    for record in data['results']['result']:
        if record['fdr'] <= fdr_threshold:

            if 'id' in record['term'].keys():     # Unclassified entries will be missing an id
                output_line = '\t'.join([
                    str(record['term']['id']),
                    str(record['term']['label']),
                    str(record['number_in_list']), 
                    str(record['fold_enrichment']),
                    str(record['fdr'])]
                )
            
                #print(output_line)
                f_out.writelines(output_line + '\n')

f.close()
f_out.close()


print('Done')