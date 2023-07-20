setwd()
vol1 <- read.csv("galaxyplot.csv", header = TRUE)
rownames(vol1) <- vol1[,1]
vol1 <- vol1[,-1]
library(EnhancedVolcano)
EnhancedVolcano(vol1, lab= rownames(vol1), x= "log2FoldChange", y= "padj", pCutoff = 0.05, FCcutoff = 1, ylim = c(0, -log10(10e-4)), xlim= c (-4, 4), pointSize = 2, labSize = 3.0, selectLab = c("hsa-mir-126", "hsa-mir-144", "hsa-mir-142", "hsa-mir-486-1", "hsa-mir-424", "hsa-mir-196a-2", "hsa-mir-196b", "hsa-mir-184", "hsa-mir-196a-1", "hsa-mir-203a"), xlab = bquote(~Log[2]~ 'fold change'), boxedLabels = TRUE, colAlpha = 4/5, cutoffLineType = 'twodash', cutoffLineWidth = 0.8, legendLabels=c('Not sig.','Log2FC','adj.p-value', 'adj.p-value & Log2FC'), legendPosition = 'top', legendLabSize = 10, legendIconSize = 5.0, drawConnectors = TRUE, widthConnectors = 1.0, colConnectors = 'black')