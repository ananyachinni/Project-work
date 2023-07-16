#Microarray Data Analysis of GSE28735
setwd = (“C:/gse28735”)
celpath = "C:/gse28735"
library(affy)
data = ReadAffy(celfile.path = celpath) #Reading in the .CEL files; throws an error asking to use oligo or xps packages for these array data types
ph = data@phenoData
ph
ph@data [ ,1] = c("mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", 
                  "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", 
                  "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", 
                  "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", 
                  "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", 
                  "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", 
                  "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control")
name = "boxplot.jpg"      #Boxplot before normalization
jpeg(name)
boxplot(data,which='pm',col='red',names=ph@data$sample,xlab='Samplename',ylab='logInt')
dev.off()
for (i in 1:90)               #MA plot for the 90 samples before normalization
   {
         name = paste("MAplot",i,".jpg",sep="")
         jpeg(name)
         MAplot(data,which=i)
         dev.off()
         }
data.rma = rma(data)  #RMA normalization
data.matrix = exprs(data.rma)
library(oligo)
name = "boxplotnorm.jpg"    #Boxplot after normalization
jpeg(name)
boxplot(data.matrix,col='red',names=ph@data$sample,xlab='Samplename',ylab='logInt')
dev.off()
for (i in 1:90)             #MA plot for the 90 samples after normalization
   {
         name = paste("MAplotnorm",i,".jpg",sep="")
         jpeg(name)
         MAplot(data.rma,which=i)
         dev.off()
         }
ph@data [ ,2] = c("mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", 
                  "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", 
                  "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", 
                  "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", 
                  "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", 
                  "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", 
                  "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control", "mutant", "control")
colnames(ph@data)[2]="source"
ph@data
groups = ph@data$source
f = factor(groups,levels=c("mutant","control"))
design = model.matrix(~ 0 + f)
colnames(design) = c("mutant","control")
library(limma)
data.fit = lmFit(data.matrix,design)
contrast.matrix = makeContrasts(mutant-control,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)
volcanoplot(data.fit.eb, coef = 1L, style = "p-value", highlight = 10, names = ph@data$source, main = "Volcano plot", xlab = "Log2 Fold Change", ylab = "-log10Pvalue", pch=16, cex=0.35)
lp <- -log10(data.fit.eb$p.value[,1])
ord <- order(lp, decreasing = TRUE)[1:500]
points(data.fit.eb$coef[ord,1], lp[ord], pch = 16, cex = 0.45, col = 2)
options(digits = 2)
tab = topTable(data.fit.eb,coef=1,number=500,adjust.method="BH") #Using topTable for extracting differentially expressed probes
topgenes = tab[tab[, "adj.P.Val"] < 0.05, ]
topups = topgenes[topgenes[, "logFC"] > 1, ]
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
IDs.up = rownames(topups)
IDs.down = rownames(topdowns)
write.table(IDs.up,row.names=FALSE,col.names=FALSE,quote=FALSE,file="C:/gse28735/upIDs.txt")
write.table(IDs.down,row.names=FALSE,col.names=FALSE,quote=FALSE,file="C:/gse28735/downIDs.txt")
library(biomaRt)
listMarts()
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
affyids = c("7908072","7924029","8049487","7901175","7944164","8045835","7981514","8093950","7920258","8169504","8081818","7996819","8087547",
            "7915472","7961215","7936144","7970441","7909164","8129082","8152703","7993588","8050160","7979710","8028924","8105267","8020551",
            "7950933","8015349","8012953","8150889","8009951","7915184","8029086","7962579","7918064","8021584","8139488","7985317","8027748",
            "7938485","7930980","8138381","8002303","8083941","7952290","8120043","7969493","7989985","8059350","8112045","7993606","8064904",
            "7899265","8067125","7971077","8029098","7909332","7920291","8008237","8072735","8065948","7976443","7979179","7919669","8053417",
            "7920128","8058765","8015016","8071758","7954065","7978544","8106727","7965541","8086517","8113433","7956271","7990151","8135601",
            "7910611","8105300","7969438","8111387","8102232","7942135","7934898","7950391","8029773","7902127","7972713","7945245","8084794",
            "8042942","8066493","7909708","8037205","8056184","8098246","7935058","7950810","8132318","8149774","8042439","8126820","8140967",
            "8018305","7954527","7973336","8048717","8062927","8146863","8092541","7919637","8102950","8148572","8118061","7928429","8162744",
            "7958884","8112342","7964907","8015387","8138613","8065412","8054479","8169263","7920297","8104788","8179221","8101429","8082928",
            "8139207","7962183","8065416","8085984","8078014","8090162","8049349","8094028","7898413","7927367","8146000","8056113","8047738",
            "7930498","8097692","8059279","8061579","8056151","8180414","8179617","7953532","8100298","8004184","8014974","8120719","8070579",
            "8051322","8048864","8129880","8138566","8033674","8179861","7923086","8127563","8130867")
getBM(attributes=c('affy_hugene_1_0_st_v1', 'entrezgene', 'hgnc_symbol'), filters = 'affy_hugene_1_0_st_v1', values = affyids, mart = ensembl)
affyids = c("8047300","8054281","7933194","7961580","8017098","8145532","7922598","8057056","8149927","8111677","8160168","8113666","8144786",
            "7910111","7925342","8146687","8131927","8057377","8166784","7908917","7904361","8098654","7976073","7976350","8060134","7917304",
            "7924309","8106827","7909681","7954377","8135378","8156134","7930454","8161755","7934936","8113234","7986195")
getBM(attributes=c('affy_hugene_1_0_st_v1', 'entrezgene', 'hgnc_symbol'), filters = 'affy_hugene_1_0_st_v1', values = affyids, mart = ensembl)