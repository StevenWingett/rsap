import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os


# Setup
def read_options():
    #Arguments and options
    parser = argparse.ArgumentParser(
                prog = 'summarising_deseq_results.py',
                description = 'Python script to summarise the results of the DESeq2 output',
                epilog = 'Steven Wingett 2024, The MRC-LMB, Cambridge, UK'
    )
    parser.add_argument("--metadata_file", action="store", type=str, metavar='', 
                        help="Path to the DESeq2 metadata file - generated previously by the pipeline i.e. NOT the design file")
    
    parser.add_argument("--deseq_file", action="store", type=str, metavar='',
                        help="Path to the DESeq output file")
    
    parser.add_argument("--log2_norm_ex_file", action="store", type=str, metavar='', 
                        default='results_rsap/expression_data_pipeline_format/log2_normalised_expression_data__plus_1_pipeline_format.tsv.gz',
                        help="Path to log2 normalised expression file - produced previously in the pipeline"
                        )

    parser.add_argument("--padj_threshold", action="store", type=float, metavar='', default=0.05, 
                        help="DESeq2 maxmium adjusted p-value threshold"
                        )
    
    parser.add_argument("--abs_l2fc_threshold", action="store", type=float, metavar='', default=1, 
                        help="Minimum DESeq2 absolute log2-fold change threshold"
                        )

    parser.add_argument("--outdir", action="store", type=str, metavar='', help="Output directory", default='./deseq2_summarising')

    args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not
    options = args[0]   # Get the 2 arrays of known/unknown arguments from the tuple

    return options

options = read_options()


#metadata_file = 'results_rsap/metadata_deseq2/comparison1_metadata.tsv'
#deseq_file = 'results_rsap/deseq2/comparison1__DAY_Zatch__Day_0_vs_Day_6/Comparison_Day_6_vs_Day_0/Comparison_Day_6_vs_Day_0.deseq2_results.tsv'
#log2_normalised_expression_file = 'results_rsap/expression_data_pipeline_format/log2_normalised_expression_data__plus_1_pipeline_format.tsv.gz'
#padj_threshold = 0.05
#abs_l2fc_threshold = 1
#outdir = './deseq2_summarising'

image_formats = ('png', 'svg', 'eps')


# Import data
print(f'Reading in expression file: {options.log2_norm_ex_file}')
log2_normalised_expression_data = pd.read_csv(options.log2_norm_ex_file, sep='\t')
print(f'{log2_normalised_expression_data.shape[0]} genes and {log2_normalised_expression_data.shape[1] - 2} samples in count matrix')
print()

print(f'Reading in metadata file: {options.metadata_file}')
metadata = pd.read_csv(options.metadata_file, sep='\t')
print(f'{metadata.shape[0]} samples present')
print()

print(f'Reading in DESeq2 data file: {options.deseq_file}')
deseq_data= pd.read_csv(options.deseq_file, sep='\t')
print(f'{deseq_data.shape[0]} genes reported')


# Create and identifier for the type of comparison being performed
sample_types = metadata.iloc[:, 1].drop_duplicates().reset_index(drop=True)
contrasts = ''

if metadata.shape[1] > 2:
    contrasts = metadata.columns.to_list()[2:]
    contrasts = '_'.join(contrasts)
    contrasts = f'__({contrasts})'

comparison = f'{sample_types[1]}_vs_{sample_types[0]}{contrasts}'


# Filter log2 normalised data to only contain subsets of interest
filt = pd.Series(log2_normalised_expression_data.columns.isin(metadata['Pipeline_Name'])).to_list()
filt[0:2] = (True, True)  # ID columns
log2_normalised_expression_data = log2_normalised_expression_data.iloc[:, filt]
print(f'{log2_normalised_expression_data.shape[1] - 2} samples after filtering using metadata')


# Make output directory
if not (os.path.exists(options.outdir)):
    os.mkdir(options.outdir)


# List DE genes
up_regulated_genes = deseq_data.query('padj <= @options.padj_threshold')
down_regulated_genes = up_regulated_genes.query('log2FoldChange <= -@options.abs_l2fc_threshold').loc[:, 'region'] 
up_regulated_genes = up_regulated_genes.query('log2FoldChange >= @options.abs_l2fc_threshold').loc[:, 'region']

outfile = f'{options.outdir}/up_regulated_genes.tsv'
print(f'Writing to {outfile}')
up_regulated_genes.to_csv(outfile, sep='\t', index=False)

outfile = f'{options.outdir}/down_regulated_genes.tsv'
print(f'Writing to {outfile}')
down_regulated_genes.to_csv(outfile, sep='\t', index=False)


# Volcano plot
#deseq_data['minus_log10(padj)'] = -np.log10(deseq_data['padj'])

#sns.scatterplot(data=deseq_data, 
#                x="log2FoldChange", 
#                y="minus_log10(padj)", 
#                s=2
#               )

#plt.title(comparison)
#plt.axhline(y = -np.log10(options.padj_threshold), color = 'r', linestyle = '--', lw=0.5) 
#plt.axvline(x = options.abs_l2fc_threshold, color = 'r', linestyle = '--', lw=0.5) 
#plt.axvline(x = -options.abs_l2fc_threshold, color = 'r', linestyle = '--', lw=0.5) 

#outfile = f'{options.outdir}/{comparison}.volcano_plot'
#for image_format in image_formats:
#    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
#plt.clf()
#plt.show()


# Identify differentially expressed genes in log2 data file
log2_normalised_expression_data['DEG'] = 'NO'

filt = log2_normalised_expression_data['gene_id'].isin(down_regulated_genes)
log2_normalised_expression_data.loc[filt, 'DEG'] = 'DOWN'

filt = log2_normalised_expression_data['gene_id'].isin(up_regulated_genes)
log2_normalised_expression_data.loc[filt, 'DEG'] = 'UP'


# Volcano Plots
volcano_data = deseq_data.copy()
volcano_data['minus_log10(padj)'] = -np.log10(volcano_data['padj'])
volcano_data = volcano_data.loc[:, ['region', 'log2FoldChange', 'minus_log10(padj)']]


# Remove entries containing NA
print(f'{volcano_data.shape[0]} genes BEFORE filter for NA')
filt = ~ volcano_data.isna().any(axis=1)
volcano_data = volcano_data[filt]
print(f'{volcano_data.shape[0]} genes AFTER filtering for NA')

volcano_data['Significant'] = 'NO'

filt = volcano_data['region'].isin(down_regulated_genes)
volcano_data.loc[filt, 'Significant'] = 'DOWN'

filt = volcano_data['region'].isin(up_regulated_genes)
volcano_data.loc[filt, 'Significant'] = 'UP'

log2FoldChange_off_scale = 10
minus_log10_padj_off_scale = 10

volcano_data['Off_scale'] = False

filt = volcano_data['log2FoldChange'] > log2FoldChange_off_scale
volcano_data.loc[filt, 'log2FoldChange'] = log2FoldChange_off_scale
volcano_data.loc[filt, 'Off_scale'] = True

filt = volcano_data['log2FoldChange'] < -log2FoldChange_off_scale
volcano_data.loc[filt, 'log2FoldChange'] = -log2FoldChange_off_scale
volcano_data.loc[filt, 'Off_scale'] = True

filt = volcano_data['minus_log10(padj)'] > minus_log10_padj_off_scale
volcano_data.loc[filt, 'minus_log10(padj)'] = minus_log10_padj_off_scale
volcano_data.loc[filt, 'Off_scale'] = True

filt = volcano_data['minus_log10(padj)'] < -minus_log10_padj_off_scale
volcano_data.loc[filt, 'minus_log10(padj)'] = -minus_log10_padj_off_scale
volcano_data.loc[filt, 'Off_scale'] = True

volcano_data = volcano_data.reset_index(drop=True) #The needs doing

# Make volcano plot with no annotations
if volcano_data['Off_scale'].sum() > 0:   # Prevents error
    markers = ['o', '*']
else:
    markers = ['o']

colors = ["grey", "red", 'blue']
sns.set_palette(sns.color_palette(colors))

sns.scatterplot(data=volcano_data, 
                x="log2FoldChange", 
                y="minus_log10(padj)", 
                hue="Significant",
                hue_order=['NO', 'UP', 'DOWN'],
                style="Off_scale",
                markers=markers,
                s=7,
                edgecolor = None
               )
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

outfile = f'{options.outdir}/{comparison}.volcano_plot'
#outfile = f'{options.outdir}/{os.path.basename(deseq_data_file)}.padj{padj_threshold}_abs_l2fc_{abs_l2fc_threshold}_volcano_plot'
for image_format in image_formats:
    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
plt.clf()










# Cluster heatmap of de genes
cluster_heatmap_data = log2_normalised_expression_data.query('DEG != "NO"')
cluster_heatmap_data = cluster_heatmap_data.drop(['gene_id', 'DEG'], axis=1)
cluster_heatmap_data = cluster_heatmap_data.set_index('gene_name')


try:    # clustermap will fail if not many DE genes
    cg = sns.clustermap(data=cluster_heatmap_data, 
                        z_score=0,
                        col_cluster=False,
                        xticklabels=True, 
                        yticklabels=True,
                        figsize=(10, 20),
                        #cbar_pos=(0.01, 0.95, 0.01, 0.01)
                    )

    plt.xticks(fontsize=30, rotation=0)
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize = 3)
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize = 1)

    outfile = f'{options.outdir}/{comparison}.deg_cluster_heatmap'
    for image_format in image_formats:
        plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
    plt.clf()
    #plt.show()
except:
    print(f'Unable to create clustermap for {comparison} - perhaps there is limited data returned?')



# Prepare scatter plot data
scatterplot_data = (log2_normalised_expression_data
                    .copy()
                    .drop('gene_name', axis=1)
                   )

# Use metadata to classify samples into groups
scatterplot_data = pd.melt(scatterplot_data, id_vars=['gene_id', 'DEG'], var_name='Pipeline_Name', value_name='expression')
scatterplot_data = pd.merge(scatterplot_data, metadata.iloc[:, 0:2], how='left', on='Pipeline_Name')

# Rename columns and remove Pipeline_Name column which is no longer needed
column_names = scatterplot_data.columns.tolist()
column_names[-1] = 'sample_group'
scatterplot_data.columns = column_names
scatterplot_data = scatterplot_data.drop('Pipeline_Name', axis=1)


# Calculate the mean expression for gene for each group
scatterplot_data_grouped = scatterplot_data.groupby(by=['gene_id', 'sample_group', 'DEG'])
scatterplot_data = scatterplot_data_grouped.mean().reset_index()

# Wrangle the data to scatterplot format
scatterplot_data = scatterplot_data.pivot(index=['gene_id', 'DEG'], columns='sample_group', values='expression')
scatterplot_data = scatterplot_data.reset_index()


# Plot scatterplot - plot mulitple plots so DEGs are not covered by grey spots
sns.set_theme(style="whitegrid")
colors = ["grey", "blue", 'red']
point_size=5
sns.set_palette(sns.color_palette(colors))
plt.figure(figsize=(10 ,10))


sns.scatterplot(data=scatterplot_data[scatterplot_data['DEG'] == 'NO'], 
                x=sample_types[0], 
                y=sample_types[1], 
                legend=False,
                linewidth=0,
                s=point_size,
                alpha=1,
                color = 'grey'
               )

sns.scatterplot(data=scatterplot_data[scatterplot_data['DEG'] == 'UP'], 
                x=sample_types[0], 
                y=sample_types[1], 
                legend=False,
                linewidth=0,
                s=point_size,
                alpha=1,
                color = 'red'
               )

sns.scatterplot(data=scatterplot_data[scatterplot_data['DEG'] == 'DOWN'], 
                x=sample_types[0], 
                y=sample_types[1], 
                legend=False,
                linewidth=0,
                s=point_size,
                alpha=1,
                color = 'blue'
               )

plt.title(comparison)
#plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

outfile = f'{options.outdir}/{comparison}.scatterplot'
for image_format in image_formats:
    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
plt.clf()
#plt.show()


print('Done')
