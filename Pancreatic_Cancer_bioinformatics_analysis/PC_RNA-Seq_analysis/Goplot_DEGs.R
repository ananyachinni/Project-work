setwd()
library(readxl)
david <- read_excel("GO_upreg.xlsx")
david <- as.data.frame(david)
genelist <- read.csv("degs.csv")
gene <- read.csv("Chord_Go2.csv")
process <- c("inflammatory response",
             "B cell receptor signaling pathway",
             "signal transduction",
             "receptor internalization",
             "cell adhesion",
             "cytoskeleton",
             "myosin complex",
             "lamellipodium",
             "mast cell granule", "actin binding", "receptor activity", "low-density lipoprotein particle binding", "guanyl-nucleotide exchange factor activity")
list <- list("david" = david, "genelist" = genelist, "genes" = gene, "process" = process)
library(GOplot)
circ <- circle_dat(list$david, list$genelist)
chord <- chord_dat(data = circ, genes = list$genes, process = list$process)
cho <- GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 6, ribbon.col = mycolors)
cho
clust <- GOCluster(chord, list$process, clust.by = 'term', term.width = 2, term.col = mycolors)
clust
nb.cols <- 7
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
ggsave("GOChord_red_pre5.png", cho, height = 12, width = 4 * 3)
