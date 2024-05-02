# Python script to perform QC on the expression dataset

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import math
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import plotly.express as px
from sklearn.decomposition import PCA
import argparse


# Setup
#expression_data_file = 'expression_data_pipeline_format.tsv.gz'
image_formats = ('png', 'svg', 'eps')
#outdir = 'QC'
#analysis_design_file = 'analysis_design_file.csv'
palette_to_choose = 'colorblind'
principal_components_to_select = 4


def read_options():
    #Arguments and options
    parser = argparse.ArgumentParser(
                prog = 'qc.py',
                description = 'Python script to perform basic QC on mapped and quantitated RNA-seq data matrix files',
                epilog = 'Steven Wingett 2024, The MRC-LMB, Cambridge, UK'
    )
    parser.add_argument("--expression_file", action="store", type=str, metavar='', default='expression_data_pipeline_format.tsv.gz',
                        help="Path to expression matrix (in pipeline format)")
    parser.add_argument("--design_file", action="store", type=str, metavar='', help="CSV experiment design file", default=None)
    parser.add_argument("--outdir", action="store", type=str, metavar='', help="Output directory", default='QC')

    args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not
    options = args[0]   # Get the 2 arrays of known/unknown arguments from the tuple

    return options


options = read_options()
outdir = f'{options.outdir}/QC'


#Read in data 
print('Reading in expression data: ' + options.expression_file)
expression_data = pd.read_csv(options.expression_file, sep='\t')
print(f'{expression_data.shape[1] - 2} samples with {expression_data.shape[0]} genes')


colorDict = {}

if options.design_file is not None:
    analysis_design = pd.read_csv(options.design_file)

    # Check no duplicate values
    if analysis_design.loc[:, ['Sample', 'Pipeline_Name']].duplicated().any():
        print('Design file contains duplicated Sample/Pipeline_Name combinations')
    
    # Check all values in expression data
    print('Identifying samples of interest')
    not_found = analysis_design['Sample'].isin(expression_data.iloc[:, 2:])
    not_found = ~ not_found.any()
    if not_found:
        print('Warning: sample name(s) found in the design file that are not in the expression data file')
    
    analysis_design['Sample']
        
    filt = (expression_data
            .iloc[:, 2:]
            .columns.isin(analysis_design['Sample'])
            .tolist()
                         )
    filt = [True, True] + filt   # Select first 2 columns
    expression_data = expression_data.iloc[:, filt]
    
    print(f'{expression_data.shape[1] - 2} samples remaining') 
    
    # Rename the column headers
    sample_lookup_dict = analysis_design.set_index('Sample')['Pipeline_Name'].to_dict()
    expression_data.rename(columns=sample_lookup_dict)
    
    # Set up the colour scheme
    palette_to_choose = sns.color_palette(palette_to_choose, analysis_design.shape[0]).as_hex()
    palette_ids = pd.Categorical(analysis_design['Group']).codes

    for i in range(analysis_design.shape[0]):
        colorDict[analysis_design['Sample'][i]] = palette_to_choose[palette_ids[i]]
else:
    sample_names = expression_data.columns.tolist()[2:]
    palette_to_choose = sns.color_palette(palette_to_choose, len(sample_names))
    # Set up the colour scheme - every sample a different colour
    for i in range(len(sample_names)):
        colorDict[sample_names[i]] = palette_to_choose[i]


# Remove rows with no data
filt = expression_data.iloc[:,2:].sum(axis=1) != 0
expression_data = expression_data[filt]

print(f'{expression_data.shape[0]} genes after removing genes with no expression values')

# Convert to log2
print('Convert to log2(expression + 1)')
expression_data.iloc[:, 2:] = np.log2(expression_data.iloc[:, 2:] + 1)

# Make the ouput directories
output_directories = ['histograms', 'cumulative_distribution', 'pca',
                      'correlation_heatmap', 'similarity_tree']

for output_directory in output_directories:
    output_directory = f'{outdir}/{output_directory}'
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        

# Histograms
#Process
fig_width = 20
ncol = 4
ngraph = expression_data.shape[1] - 3
nrow = int(math.ceil(ngraph / ncol))
fig_height = fig_width * (nrow / ncol)
fig = plt.figure(figsize=(fig_width, fig_height))
xmin = expression_data.iloc[:, 3:].min().min()
xmax= expression_data.iloc[:, 3:].max().max()

counter = 0

for i in range(3, expression_data.shape[1]):
    
    # Graph code1
    sample_name = expression_data.iloc[:, i].name
    counter += 1
    ax = fig.add_subplot(nrow, ncol, counter)
    ax.hist(x=expression_data.iloc[:, i], bins=15, color=colorDict[sample_name])    
    ax.set_title(sample_name)
    ax.set_xlim(xmin, xmax)
 
outfile = f'{outdir}/histograms/expression_histogram'
for image_format in image_formats:
    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
plt.clf()


# Combined data histogram
sns.histplot(data=pd.Series(expression_data.iloc[:, 2:].to_numpy().flatten()), bins=30, color='blue')

outfile = f'{outdir}/histograms/combined_expression_histogram'
for image_format in image_formats:
    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
plt.clf()


# Cumulative distribution plot
custom_palette = []
for sample in expression_data.iloc[:, 2:].columns.tolist():
    custom_palette.append(colorDict[sample])

plt.figure(figsize=(7,7))
plot=sns.ecdfplot(data=expression_data.iloc[:, 2:], legend=False, palette=custom_palette)
plt.xlabel('Log2(Normalised expression + 1)')

outfile = f'{outdir}/cumulative_distribution/cumulative_distribution'
for image_format in image_formats:
    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
plt.clf()


# Correlation heatmap
pearson_matrix = expression_data.iloc[:, 2:].corr()
plt.figure(figsize=(0.5 * expression_data.shape[1], 0.4 * expression_data.shape[1]))
sns.heatmap(pearson_matrix, annot=True);
plt.title('Correlation heatmap of datasets')

#plt.savefig('correlation_heatmap_before_filtering.svg', bbox_inches='tight')

outfile = f'{outdir}/correlation_heatmap/expression_correlation_heatmap'
for image_format in image_formats:
    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
plt.clf()


# Make dendrogram
pearson_matrix.index = pearson_matrix.columns
#matrix
dissimilarity = 1 - abs(pearson_matrix)

Z = linkage(squareform(dissimilarity), 'average')

dendrogram(Z, labels=pearson_matrix.index, orientation='left', color_threshold=0, above_threshold_color='black')

if options.design_file is not None:
    ax = plt.gca()
    y_labels = ax.get_ymajorticklabels()
    for y in y_labels:
        color_to_lookup = y.get_text()
        y.set_color(colorDict[color_to_lookup])

outfile = f'{outdir}/similarity_tree/expression_correlation_heatmap'
for image_format in image_formats:
    plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
plt.clf()


# Format PCA data
pca_data = expression_data.iloc[:, 2:].transpose()

sample_labels = (pca_data
          .index
          .to_series()
          .reset_index(drop=True)
         )

pca = PCA()
components = pca.fit_transform(pca_data)
labels = {
    str(i): f"PC {i+1} ({var:.1f}%)"
    for i, var in enumerate(pca.explained_variance_ratio_ * 100)
}

scatter_plot_data = pd.DataFrame(components)

column_names = []
for i in range(1, scatter_plot_data.shape[1] + 1):
    column_names.append(f'PC{i}')
scatter_plot_data.columns = column_names    

scatter_plot_data['Sample'] = pca_data.index
scatter_plot_data = pd.concat([scatter_plot_data.iloc[:, -1:], scatter_plot_data.iloc[:, :-1]], axis=1)    # Re-order columns 

# make PCA plots
for i in range(1, principal_components_to_select):
    for j in range(i + 1, principal_components_to_select + 1):
        pcx = i
        pcy = j
        marker_size=300

        # Format data for graph
        pc_variance_x = round(pca.explained_variance_ratio_[pcx - 1] * 100, 1)
        pc_variance_y = round(pca.explained_variance_ratio_[pcy - 1] * 100, 1)

        pcx = f'PC{pcx}'
        pcy = f'PC{pcy}'

        x_axis_label = f'{pcx} ({pc_variance_x}%)'
        y_axis_label = f'{pcy} ({pc_variance_y}%)'

        # Plot graph
        custom_palette = []
        for sample in scatter_plot_data.loc[:, 'Sample'].drop_duplicates():
            custom_palette.append(colorDict[sample])

        sns.scatterplot(data=scatter_plot_data, x=pcx, y=pcy, hue='Sample', s=marker_size, palette=custom_palette)
        plt.xlabel(x_axis_label)
        plt.ylabel(y_axis_label)
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        
        outfile = f'{outdir}/pca/pc{i}_pc{j}_plot'
        for image_format in image_formats:
            plt.savefig(fname=f'{outfile}.{image_format}', bbox_inches='tight', pad_inches=0.5)
        plt.clf()

print('Done')