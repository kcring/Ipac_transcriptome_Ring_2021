#########################################################################################################################
# purpose: generate DESeq2 results objects of DE genes
# input: path to directory which contains (at any level) quant.sf files from the same experiment (and only the same experiment), path to tx2gene, list of sample_rep, sample_conditions, path to functions script
# output: results tables, MA plots
# written by: chase mateusiak, chase.mateusiak@gmail.com
# modified by: Anne Sapiro, annesapiro@gmail.com and Kacie Ring, kacie.ring@gmail.com
# credit: https://bioconductor.github.io/BiocWorkshops/rna-seq-data-analysis-with-deseq2.html#exploratory-data-analysis
#
#########################################################################################################################
#
# experiment:Host blood meal experiment- Ixodes pacificus mapped to ISE6 (Ixodes scapularis)
#
################################################## user input #############################################################
getwd()
# path to working directrory, which should include a directory `counts_for_deseq` that includes raw counts files
mappingDir <- '/Users/kaciering/Desktop/'
# point to directory with raw read counts, in this case called `counts_for_deseq`
quantDir <- file.path(mappingDir, 'Bloodmeal_quants')


# add an output directory called analysis
analysisDir <- file.path(mappingDir, 'analysis')
dir.create(analysisDir)

############################################### end user input ##############################################################
#############################################################################################################################
############################################### install packages ###########################################################

# description of packages
# DESeq2 -- differential gene expression, used throughough
# ggplot2 -- used specifically for PCA plot, but generally has nice plotting functionality

# libraries will be installed and loaded via the librarian package (https://cran.r-project.org/web/packages/librarian/vignettes/intro-to-librarian.html)
# you need the librarian package and Biobase installed for this to work. if needed, uncomment lines below and install

#install.packages("librarian") # uncomment to use

# biobase, for installing packages from Bioconductor

#if (!requireNamespace("BiocManager", quietly = TRUE))  #uncomment these 3 lines to install
#  install.packages("BiocManager")
#BiocManager::install("Biobase")

librarian::shelf(ggplot2, ggfortify, tidyverse, DESeq2, pheatmap, RColorBrewer, cluster, broom, ggforce, EnhancedVolcano)


################################################### start main script #######################################################

# set working directory
setwd(mappingDir)

#import counts files (finds all files ending in counts.cut.txt in count directory)
sampleFiles <- grep('counts.cut.txt', list.files(quantDir), value = TRUE)
#print out the files found
sampleFiles

#for each file (in alphabetical order, see prited out list), write what group it belongs to
sampleCondition<-c("UF_MF","UF_MF","UF_MF","UF_LF","UF_LF","UF_LF", "posBb_MF","posBb_MF","posBb_MF","posBb_LF","posBb_LF","posBb_LF","negBb_MF","negBb_MF","negBb_MF","negBb_LF","negBb_LF","negBb_LF")
#create sample table
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

#double check to make sure it all looks good
sampleTable

# pull all of the info into deseq
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=quantDir, design=~condition)


#### look at dds object
colSums(assay(ddsHTSeq))
colData(ddsHTSeq)
rowData(ddsHTSeq)

############################### use heatmap, pca and dispersion to visualize data ########################################
setwd(analysisDir)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
#### check the percent of genes without zero counts.

# count number of genes
geneCounts <- counts(ddsHTSeq)

# make logical maxtrix of genes with (true) and without (false) zero count
indexNotZero <- apply(geneCounts, 1, function(x) {all (x > 0)})

# display percent of genes with zero counts
sum(indexNotZero == TRUE)/nrow(assay(ddsHTSeq))

#### create heatmap to examine similiarities btwn samples

#normalize data - log transformation
rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
#str(rld)

#calculate sample-to-sample distances  
distsRL <- dist(t(assay(rld)))
distsRL

#create a matrix of sample-tosample distances 
mat <- as.matrix(distsRL)
rownames(mat) <- colData(rld)$condition
colnames(mat) <- colData(rld)$condition
#view(mat)

#generate heatmap based on sample-to-sample distances 
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
pheatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
heatmap <- pheatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
pdf(file="heatmap.pdf")
heatmap
dev.off() 

# PCA
# from the vignette re: vst : DESeq2 offers transformations for count data that stabilize the variance across the mean: the regularize logarithm (rlog) and the variance stabilizing transformation (VST). 
# These have slightly different implementations, discussed a bit in the DESeq2 paper and in the vignette, but a similar goal of stablizing the variance across the range of values. 
# Both produce log2-like values for high counts

#normalize data using variance stabilizing transformation 
vsd <- vst(ddsHTSeq)


#apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
#create a matrix - columns=variables and rows=samples
distsVS <- dist(t(assay(vsd)))
pca_matrix <- as.matrix(distsVS)
rownames(pca_matrix) <- colData(vsd)$condition
colnames(pca_matrix) <- colData(vsd)$condition

# Transpose the matrix so that rows = samples and columns = variables
pca_matrix <- t(pca_matrix)
view(pca_matrix)
sample_pca <- prcomp(pca_matrix)

#look at the first 16 rows and columns 
pca_matrix[1:16, 1:16]

#convert matrix into a tibble 
as_tibble(pca_matrix)
as_tibble(pca_matrix, rownames = "sample")

#tibble matrix must have unique column names. Use name_repair function to indicate replicates of the condition
tibble::as_tibble(diag(3), .name_repair = "unique")

#Eigenvalues - these represent the variance explained by each PC. 
#We can use these to calculate the proportion of variance in the original data that each axis explains.
tidy(sample_pca, matrix = "eigenvalues")

#This table can easily be used to produce a Scree Plot, which shows the fraction of total variance explained by each principal component. 
tidy(sample_pca, matrix = "eigenvalues") %>% 
  ggplot(aes(x = factor(PC))) +
  geom_col(aes(y = percent)) +
  geom_line(aes(y = cumulative, group = 1)) + 
  geom_point(aes(y = cumulative)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# The PC scores are stored in the "x" value of the prcomp object, which is a matrix
pc_scores <- sample_pca$x 

#calculate percent PCA 
percentage <- round(sample_pca$sdev / sum(sample_pca$sdev) * 100, 2)
percentage <- paste( colnames(sample_pca), "(", paste( as.character(percentage), "%", ")", sep="") )

# Make PCA plot
pca_plot <- pc_scores %>% 
  # convert it to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # create the plot
  ggplot(aes(PC1,PC2, color = sample)) +
  geom_point() +
  expand_limits(x = -250, y = -210) +
  xlab(percentage[1]) + 
  ylab(percentage[2]) +
  #use ggforce when < 4 replicates per sample condition
  ggforce::geom_mark_ellipse(aes(fill = sample,
                                 color = sample))

pdf(file="PCA.pdf")
pca_plot
dev.off() 

### plot dispersion
ddsHTSeq_disp <- estimateSizeFactors(ddsHTSeq)
ddsHTSeq_disp <- estimateDispersions(ddsHTSeq_disp)
plotDispEsts(ddsHTSeq_disp, main = "Dispersion")

############################################# begin DE analysis #####################################################

# create DESeq object
dds <-  DESeq(ddsHTSeq, betaPrior = TRUE)

########################### create results objects. There are two ways to do this: ##################################
# 
# pairwise comparison for differential expression (things that increase in PScoc will be + numbers)
#Comparison below is for Unfed nymphs with either a lizard(LF) of mouse(MF) larval bloodmeal

res_UFlizard_UFmouse <- results (dds, contrast = c('condition', 'UF_LF','UF_MF'))
res_UFlizard_UFmouse <- res_UFlizard_UFmouse[order(res_UFlizard_UFmouse $padj),]
head(res_UFlizard_UFmouse)
# make an MA plot to visualize data

plotMA(res_UFlizard_UFmouse, ylim=c(-8,8), main="Lizard-fed larvae vs. Mouse-fed larvae")

pdf(file="res_UFlizard_UFmouse_maplot.pdf")
plotMA(res_UFlizard_UFmouse, ylim=c(-5,5), main="Lizard-fed larvae vs. Mouse-fed larvae")
dev.off() 

# print out the results file as a csv
write.csv(as.data.frame(res_UFlizard_UFmouse),file="res_UFlizard_UFmouse.csv", quote=F)
#####################################################################
# you can copy and paste ^ code for any other pairwise comparisons, just change the groups/names

#for unfed vs. engorged, we changed sample conditions (in Line 58) to 
#sampleCondition<-c("Unfed","Unfed","Unfed","Unfed","Unfed","Unfed", "Engorged","Engorged","Engorged","Engorged","Engorged","Engorged","Engorged","Engorged","Engorged","Engorged","Engorged","Engorged")
#our deseq comparison was as follows 
#res_Engorged_Unfed <- results (dds, contrast = c('condition','Engorged', 'Unfed'))



# make a file of normalized counts, which can help look at genes over time course
normCounts <- counts(dds, normalized=TRUE)
write.csv(normCounts, file='Host_Legacy_normCounts.csv', quote=F)
#####################################################################
# you can copy and paste ^ code for any other pairwise comparisons, just change the groups/names

# these snippet make a graph of 1 gene over all samples
#plotCounts(dds, gene='[Accession number of gene of interest]')
# don't forget to save graphs as files, save the output csv files, and your sessionInfo() so you don't have to remember which version of software you were using


###Volcano plot script####
#Volcano plots visualize differentially expressed genes in the above pairwise comparision

# below pulls out only the significant values and writes them to a .csv
res_UFlizard_UFmouse_sig<- subset(res_UFlizard_UFmouse, padj <= 0.1) 
write.csv(res_UFlizard_UFmouse_sig, file ='res_UFlizard_UFmouse_sig.csv') 

dataframe <- res_UFlizard_UFmouse

head(dataframe)

#Enhanced volcano plot script

EnhancedVolcano (dataframe,
                 lab = rownames(dataframe),
                 x = 'log2FoldChange',
                 y = 'padj',
                 drawConnectors =FALSE,
                 colConnectors = '#d55e00',
                 transcriptPointSize = 3.0, ## changes size of points
                 title = 'results_of_deseq_comparison', 
                 xlim = c(-5, 5), ## x axis limit
                 ylim = c(0,10),## y axis limit 
                 xlab = bquote(~Log[2]~ 'fold change'), ## label of x axis
                 ylab = bquote(~-Log[10]~adjusted~italic(P)), ## label of y axis
                 pCutoff = 0.1, ## p values that are being considered significant
                 col = c('grey30', '#56B4E9', '#0072B2', '#009E73'), ## uses a color blind friendly pallete from here http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
                 colAlpha = 0.7, ## from 0 to 1 this changes how opaque the points are
                 cutoffLineType = 'dotted',
                 gridlines.major = FALSE, ## gets rid of grid marks from the background
                 gridlines.minor = FALSE)

ggsave('results_of_deseq_comparison_plot_logandPadj_labeled2.pdf', width = 15, height = 12, units = c("in"),
       dpi = 300)


