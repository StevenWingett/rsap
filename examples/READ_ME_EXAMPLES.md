# Examples

The file "rsap_example_design_file.csv" is an example design file for rsap.

The first 3 column need to be named:

Sample	- Sample name as it appears in the input matrix data files
Pipeline_Name - How the sample should be named in the pipeline - TODO - at present column names won't be changed so for now make sure Pipeline_Name == Sample for all samples
Group - How to group samples plots (e.g. PCA plots), but not for DEseq2 analysis

Then, if required, we list the DESeq comparisons:
Comparison name in the column header
DEseq baseline factor in the column name (after forward slash)
Then list the groups that need to be compared to the baseline i.e. allocate samples to groups
Add any contrasts in an adjacent column

Separate different comparisons with a blank column.