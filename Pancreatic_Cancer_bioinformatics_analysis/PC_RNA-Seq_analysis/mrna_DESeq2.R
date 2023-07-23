## Reference URL: https://hbctraining.github.io/DGE_workshop/schedule/1.5-day.html
setwd()
#Uploading metadata and counts data files into the environment
meta <- read.csv('meta.csv', header = TRUE, sep = ",", quote = "\"", fill = TRUE, comment.char = "")
rownames(meta) <- meta[ , 1]
meta <- meta[ , -1]
data <- read.csv('symbol.csv', header = FALSE, sep = ",", quote = "\"", dec = "-", fill = TRUE, comment.char = "")
rownames(data) = make.names(data[,1], unique = TRUE)
data <- data[ , -1]
names(data) <- lapply(data[1, ], as.character)
data <- data[-1,]
d <- as.matrix(data)
storage.mode(d) = "integer" #Storing the counts data matrix in integer format
#Inspecting the match between row names of the metadata file and column names of the counts data matrix
all(colnames(d) %in% rownames(meta))
all(colnames(d) == rownames(meta))
#Since they don't match, rearrangement of the mismatched rows/columns are done using match function
idx <- match(d,meta)
idx <- match(colnames(d), rownames(meta))
metam <- meta[idx, ] #Rearranging the meta variable based on the new index
#Checking the match again, since the values are TRUE, we move on to the next step
all(colnames(d) == rownames(metam))
all(colnames(d) %in% rownames(metam))
#Loading the required libraries
library(tidyverse)
library(DESeq2)
library(DEGreport)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(tibble)
#Creating a DESeq matrix from the counts data
dds <- DESeqDataSetFromMatrix(countData = d, colData = metam, design = ~cases.0.samples.0.sample_type, ignoreRank = FALSE)
dds <- estimateSizeFactors(dds)
#Obtaining normalized counts data
normalized_counts <- counts(dds, normalized=TRUE)
#Transforming counts data using VST (Variance Stabilizing Transformations)
rld <- vst(dds, blind=TRUE)
#Plotting PCA plot based on sampletype
plotPCA(rld, intgroup="cases.0.samples.0.sample_type")
rld_mat <- assay(rld)
#Getting correlation matrix using the cor function
rld_cor <- cor(rld_mat)
#Generating a heatmap of the correlation matrix
pheatmap(rld_cor)
dds <- DESeq(dds)
#Plotting dispersion estimate plots
plotDispEsts(dds, legend = FALSE, cex = 0.5)
contrast_oe <- c("cases.0.samples.0.sample_type", "Primary Tumor", "Solid Tissue Normal")
res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
#Shrinking the log2foldchange estimates
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken, type = “normal”)
#Plotting MA plots for the unshrunken and shrunk data
plotMA(res_tableOE_unshrunken, ylim=c(-1,1))
plotMA(res_tableOE, ylim=c(-1,1))
#Creating tibble objects of the metadata, normalized counts and shrunk data
met <- metam %>% 
             rownames_to_column(var="samplename") %>% 
             as_tibble()
normalized_counts <- normalized_counts %>% 
             data.frame() %>%
             rownames_to_column(var="gene") %>% 
             as_tibble()
#Setting cutoff criteria for differential expression analysis
padj.cutoff <- 0.05
lfc.cutoff <- 1.00
res_tableOE_tb <- res_tableOE %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()
#Filtering differentially expressed genes from the shrunk data
sigOE <- res_tableOE_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
top20_sigOE_genes <- res_tableOE_tb %>% 
        arrange(padj) %>% 	#Arrange rows by padj values
        pull(gene) %>% 		#Extract character vector of ordered genes
        head(n=20) 
#Collecting normalized counts data of the top 20 significant genes
top20_sigOE_norm <- normalized_counts %>%
        filter(gene %in% top20_sigOE_genes)
#Collating data based on sample name key
gathered_top20_sigOE <- top20_sigOE_norm %>%
        gather(colnames(top20_sigOE_norm)[2:9], key = "samplename", value = "normalized_counts")
#Writing the above collated data into a .csv file
write.csv(gathered_top20_sigOE, "change.csv")
gathered_top20_sigOE <- read.csv('change.csv', header = TRUE, sep = ",", quote = "\"", fill = TRUE, comment.char = "")
#Doing an inner join to add the normalized counts data to the metadata tibble variable
gathered_top20_sigOE <- inner_join(met, gathered_top20_sigOE)
#Extracting normalized expression for significant genes
norm_sig <- normalized_counts[,c(1,4:9)] %>% 
             filter(gene %in% gathered_top20_sigOE$gene) %>% 
             data.frame() %>%
             column_to_rownames(var = "gene")
#Getting the corresponding annotation
annotation <- met %>% 
             select(samplename, cases.0.samples.0.sample_type) %>% 
             data.frame(row.names = "samplename")
heat_colors <- brewer.pal(6, "YlOrRd") #Setting the colour palette for the heatmap
#Plotting heatmap of the gene expression
pheatmap(norm_sig, color = heat_colors, border_color=NA, fontsize = 10, 
                    fontsize_row = 10, height=20, annotation = annotation)