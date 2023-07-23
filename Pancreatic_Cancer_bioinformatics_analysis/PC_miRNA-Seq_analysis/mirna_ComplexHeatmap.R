setwd()
library(readxl)
#Reading in the data from an .xlsx file
heat <- read_excel("heatmap_kegg_input2_5.xlsx")
heat <- as.data.frame(heat)
rownames(heat) <- heat$Pathway
heat <- heat[,-1]
heatM <- as.matrix(heat) #Creating the input matrix
#Reading in the annotation data
met <- read_excel("met.xlsx")
met <- data.frame(met)
#Loading the packages into the environment
library(ComplexHeatmap)
library(RColorBrewer)
heat_colors <- brewer.pal(9, "Pastel1") #Setting the colour palette
#Calling the pheatmap function in the ComplexHeatmap package
pheatmap(heatM, color = heat_colors, cluster_rows = FALSE, cluster_cols = FALSE, angle_col = "90", legend_breaks = c(0.01,0.02, 0.03, 0.04, 0.05), annotation_col = met)
