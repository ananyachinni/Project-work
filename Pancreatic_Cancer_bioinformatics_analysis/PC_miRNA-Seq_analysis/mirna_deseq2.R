setwd()
meta <- read.csv('metadatafiles.csv', header = TRUE, sep = ",", quote = "\"", fill = TRUE, comment.char = "")
rownames(meta) <- meta[ , 1]
meta <- meta[ , -1]
data <- read.csv('mirmatrix.csv', header = FALSE, sep = ",", quote = "\"", dec = "-", fill = TRUE, comment.char = "")
rownames(data) <- data[ , 1]
data <- data[ , -1]
names(data) <- lapply(data[1, ], as.character)
data <- data[-1,]
d <- as.matrix(data)
storage.mode(d) = "integer"
all(colnames(d) %in% rownames(meta))
all(colnames(d) == rownames(meta))
idx <- match(colnames(d), rownames(meta))
metam <- meta[idx, ]
all(colnames(d) == rownames(metam))
all(colnames(d) %in% rownames(metam))
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tidyverse)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = d, colData = metam, design = ~ cases.0.samples.0.sample_type, ignoreRank = FALSE)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
rld <- varianceStabilizingTransformation(dds, blind = TRUE)
plotPCA(rld, intgroup="cases.0.samples.0.sample_type")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor)
dds <- DESeq(dds)
plotDispEsts(dds, legend = FALSE, cex = 0.5)
contrast_oe <- c("cases.0.samples.0.sample_type", "Primary Tumor", "Solid Tissue Normal")
res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken)
#using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).

#Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#See?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#Reference: https://doi.org/10.1093/bioinformatics/bty895
plotMA(res_tableOE_unshrunken, ylim=c(-1,1))
plotMA(res_tableOE, ylim=c(-1,1))
library(ggrepel)
library(ggplot2)
met <- metam %>% 
             rownames_to_column(var="samplename") %>% 
             as_tibble()
normalized_counts <- normalized_counts %>% 
             data.frame() %>%
             rownames_to_column(var="gene") %>% 
             as_tibble()
padj.cutoff <- 0.05
lfc.cutoff <- 1.00
res_tableOE_tb <- res_tableOE %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()
sigOE <- res_tableOE_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
top20_sigOE_genes <- res_tableOE_tb %>% 
        arrange(padj) %>% 	#Arrange rows by padj values
        pull(gene) %>% 		#Extract character vector of ordered genes
        head(n=20)
top20_sigOE_norm <- normalized_counts %>%
       filter(gene %in% top20_sigOE_genes)
gathered_top20_sigOE <- top20_sigOE_norm %>%
        gather(colnames(top20_sigOE_norm)[2:9], key = "samplename", value = "normalized_counts")
write.csv(gathered_top20_sigOE, "change.csv")
gathered_top20_sigOE <- read.csv('change.csv', header = TRUE, sep = ",", quote = "\"", fill = TRUE, comment.char = "")
gathered_top20_sigOE <- inner_join(met, gathered_top20_sigOE)
#Joining, by = "samplename"
norm_sig <- normalized_counts[,c(1,4:9)] %>% 
             filter(gene %in% gathered_top20_sigOE$gene) %>% 
             data.frame() %>%
             column_to_rownames(var = "gene")
annotation <- met %>% 
             select(samplename, cases.0.samples.0.sample_type) %>% 
             data.frame(row.names = "samplename")
View(annotation)
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(norm_sig, color = heat_colors, border_color=NA, fontsize = 10, 
                                fontsize_row = 10, height=20)
res_tableOE_tb <- res_tableOE_tb %>% 
        mutate(threshold_OE = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
