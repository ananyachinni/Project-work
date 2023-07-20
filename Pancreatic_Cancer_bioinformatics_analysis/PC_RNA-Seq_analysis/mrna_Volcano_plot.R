setwd()
vol <- read.csv("Volcano_galaxy.csv", header = TRUE)
rownames(vol) <- vol[,1]
vol <- vol[,-1]
library(EnhancedVolcano)
EnhancedVolcano(vol, lab= rownames(vol), x= "log2FoldChange", y= "padj", pCutoff = 0.05, FCcutoff = 1, ylim = c(0, -log10(10e-4)), xlim= c (-4, 4), pointSize = 2, labSize = 3.0, selectLab = c("MMP12", "CD68", "LILRB1", "MATK", "CD36", "AGR2", "FXYD3", "MUC13", "ERN2", "SPDEF"), xlab = bquote(~Log[2]~ 'fold change'), boxedLabels = TRUE, colAlpha = 4/5, cutoffLineType = 'twodash', cutoffLineWidth = 0.8, legendLabels=c('Not sig.','Log2FC','adj.p-value', 'adj.p-value & Log2FC'), legendPosition = 'top', legendLabSize = 10, legendIconSize = 5.0, drawConnectors = TRUE, widthConnectors = 1.0, colConnectors = 'black')