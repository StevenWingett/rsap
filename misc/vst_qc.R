# DEseq2 comparing high / low secretors.
# Use same antibody on same day for comparisons

# # # # # # # # # # # # # # # # # #
# Real data2
library("DESeq2")
library("vsn")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
#library("apeglm")

rm(list=ls())

# Arguments
#args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
#if (length(args)!=3) {
#  stop("3 arguments need to be supplied", call.=FALSE)
#}

#count_matrix_file = args[1]
#metadata_file = args[2]
#outdir_base = args[3]

count_matrix_file <- "results_stefan_rna_seq/expression_data_pipeline_format/raw_expression_data_pipeline_format.tsv.gz"
metadata_file <- "Stefan_RNA_seq_design_file.csv"
outdir_base <- "vst_qc"

# Setup
#count_matrix_file <- "../../Quantitated_Results/salmon.merged.gene_counts.tsv.gz"
#count_matrix_file <- "./results_rsap/expression_data_pipeline_format/raw_expression_data_pipeline_format.tsv.gz"
#metadata_file <- "./results_rsap/metadata_deseq2/comparison1_metadata.tsv"


# Read in data
cts  <- read.csv(count_matrix_file, sep="\t", row.names="gene_id", check.names = FALSE)
cts$gene_name <- NULL
cts <- as.matrix(cts)
cts <- round(cts)   # Convert to integers

#head(cts)
#class(cts)

coldata <- read.csv(metadata_file, sep=",", check.names = FALSE)
#head(coldata)
#class(coldata)

print(paste(ncol(cts), 'datasets imported'))
coldata[,1] <- NULL   # Remove first column
coldata[,3:length(coldata)] <- NULL  # Remove trailing columns
rownames(coldata) <- coldata[,1]
coldata[,1] <- NULL

# Select data of interest
filt = colnames(cts) %in% rownames(coldata)
cts = cts[, filt ] 
print(paste(ncol(cts), 'datasets of interest selected'))

if(ncol(cts) == nrow(coldata)){
  print("All samples in metadata detected in count data matrix")
} else{
  stop("Sample(s) in metadata and not found in count data matrix!!!!!!")
}

# Filter out genes with low numbers of reads
print('Filtering out genes with a low number of reads')
print(paste(dim(cts)[1], 'genes BEFORE filtering'))
filt = rowSums(cts >= 10) >= 1   # 10 reads considered ok cutoff - see DEseq2 manual
cts <- cts[filt, ]
print(paste(dim(cts)[1], 'genes AFTER filtering'))



# Re-order the metadata so it matches the samples
coldata$dummy_row = rownames(coldata)   # Use a dummy row to prevent rownames disappearing when we have 1 column
coldata <- coldata[colnames(cts), ]
coldata$dummy_row <- NULL

# Check the data and metadata order are the same
checking = colnames(cts) != rownames(coldata)
if(any(checking)){
  stop("Samples in data and metadata do not match!!!!!!!!!!!!!")
} else{
  print("Data/metadata match correctly")
}

# Set up experimental design
reference_levels <- unique(coldata[, 1]) 
baseline_reference <- reference_levels[1]    # Reference appears first in metadata
reference_levels <- unique(reference_levels)
coldata[, 1] <- factor(coldata[, 1], levels=reference_levels)
colnames(coldata)[1] <- "Comparison"
print(colnames(coldata))

# Create output directory
if(!dir.exists(outdir_base)) {
  dir.create(outdir_base)
}

# # # # # # # # # # # # # # # # # #
# Perform analysis
comparison_column <- colnames(coldata)[1]
contrasts_columns <- ""

if(ncol(coldata) > 1){
  contrasts_columns <- paste(colnames(coldata)[-1], collapse = " + ")
  contrasts_columns <- paste(contrasts_columns, '+')
}

my_design_formula <- paste(contrasts_columns, comparison_column)
print(paste("DESeq2 design:", my_design_formula))

# The formula code below is needed to embed dynamic code in the experimental design
dds <- DESeqDataSetFromMatrix(countData = cts,
                            colData = coldata,
                            design = formula(paste("~", my_design_formula))
                            )
dds

print(colData(dds))

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

dds <- DESeq(dds)

#resultsNames(dds)

#comparison <- resultsNames(dds)[length(resultsNames(dds))]  # Last value
#print(paste("--------------- Performing comparison:", comparison, "---------------"))

#res <- results(dds, name=comparison)
#res

#library("apeglm")
#resLFC <- lfcShrink(dds, coef=comparison, type="apeglm")   # Shrunk l2fc
#resLFC <- lfcShrink(dds, coef=comparison, type="normal")
#resLFC


# # # # # # # # # # # # # # # # # #
# Generate output files
# Check res and resLFC correspond
#to_check_res <- data.frame(res[, c(-2, -3)])    # For normal
#to_check_res <- data.frame(res[, c(-2, -3, -4)])     # For apeglm, as 'stat' column not returned
#to_check_LFC <- data.frame(resLFC[, c(-2, -3)])
#to_check_res[is.na(to_check_res)] <- 'replaced'
#to_check_LFC[is.na(to_check_LFC)] <- 'replaced'

#if(identical(to_check_res, to_check_LFC)){
#  print("Comparison ok")
#} else {
#  print("DEseq and shrunkL2FC output do not correspond!!!!")
#  stop()
#}

# Combine oringal and shruk results togehter
#combined_results <- cbind(res, resLFC$log2FoldChange)
#combined_results <- cbind(combined_results, resLFC$lfcSE)
#combined_results$region <- rownames(combined_results)
#combined_results <- combined_results[, c(9, 1, 2, 3, 7, 8, 4, 5, 6)]
#colnames(combined_results)[5] = 'shrunk_Log2FoldChange'
#colnames(combined_results)[6] = 'shrunk_lfcSE'
#rownames(combined_results) <- NULL
#combined_results


# Print out results
#outdir <- comparison, day, antibody, sep="_" )
#outdir <- paste(outdir_base, comparison, sep="/")
#if(!dir.exists(outdir)) {
#  dir.create(outdir)
#}


# Write out to a file all the genes that were evaluated by DESeq2 (i.e. after filtering for low read counts)

#outfile <- paste0(outdir, '/', comparison, '.background_list.tsv')
#write.table(rownames(cts), 
#            file=outfile, 
#            row.names=FALSE,
#            quote = FALSE,
#            sep="\t")


#outfile <- paste0(outdir, '/', comparison, '.MA_plot_l2fc.png')
#print(paste('Creating MA plot', outfile))
#png(outfile)
#plotMA(res, ylim=c(-3,3))
#dev.off()

#outfile <- paste0(outdir, '/', comparison, '.MA_plot_shrunk_l2fc.png')
#print(paste('Creating shrunk lfc MA plot', outfile))
#png(outfile)
#plotMA(resLFC, ylim=c(-3,3))
#dev.off()

#outfile <- paste0(outdir, '/', comparison, '.deseq2_results.tsv')
#write.table(as.data.frame(combined_results), 
#            file=outfile, 
#            row.names=FALSE,
#            quote = FALSE,
#            sep="\t")

#####
# Further metrics

# this gives log2(n + 1)
#Effects of transformations on the variance
#The figure below plots the standard deviation of the transformed data, 
#across samples, against the mean, using the shifted logarithm transformation, 
#the regularized log transformation and the variance stabilizing 
#transformation. The shifted logarithm has elevated standard deviation in the 
#lower count range, and the regularized log to a lesser extent, while for the 
#variance stabilized data the standard deviation is roughly constant along the 
#whole dynamic range.

#Note that the vertical axis in such plots is the square root of the variance 
#over all samples, so including the variance due to the experimental 
#conditions. While a flat curve of the square root of variance over the mean 
#may seem like the goal of such transformations, this may be unreasonable in 
#the case of datasets with many true differences due to the experimental 
#conditions.
#ntd <- normTransform(dds)
vsd <- vst(dds, blind=TRUE)  #It only looks at the design to know the _global_ distribution of dispersion values. It doesn't use the design in the transformation itself.
#head(assay(vsd), 3)
#rld <- rlog(dds, blind=FALSE)


outfile <- paste0(outdir_base, '/', 'vst_normalised_data.tsv')

vsd_for_file <- as.data.frame(assay(vsd))   # Prevent R shifting column names on output
gene_id <- rownames(vsd_for_file)
vsd_for_file <- cbind(gene_id, vsd_for_file)

write.table(vsd_for_file, 
            file=outfile, 
            sep="\t", 
            quote=FALSE, 
            col.names = TRUE,
            row.names = FALSE)

#outfile <- paste0(outdir, '/', comparison, '.rlog_normalised_data.tsv')
#write.table(assay(rld), file=outfile, sep="\t", quote=FALSE)

#outfile <- paste0(outdir, '/', comparison, '.meanSdPlot_ntd.svg')
#print(paste('Creating meanSdPlot plot', outfile))
#svg(outfile)
#meanSdPlot(assay(ntd))
#dev.off()

outfile <- paste0(outdir_base, 'meanSdPlot_vsd.svg')
print(paste('Creating meanSdPlot VSD plot', outfile))
svg(outfile)
meanSdPlot(assay(vsd))
dev.off()

#outfile <- paste0(outdir, '/', comparison, '.meanSdPlot_rld.svg')
#print(paste('Creating meanSdPlot RLD plot', outfile))
#svg(outfile)
#meanSdPlot(assay(vsd))
#dev.off()

# PCA
outfile <- paste0(outdir_base, '/', 'pca.png')
print(paste('Creating PCA plot', outfile))
#svg(outfile)

pcaData <- plotPCA(vsd, intgroup=c("Comparison"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Comparison)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggsave(outfile)
#dev.off()

# Cluster heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

outfile <- paste0(outdir_base, '/', 'distanceMatrixHeatmap.png')
print(paste('Creating clusterheatmap', outfile))
png(outfile)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
      
print("Done")
