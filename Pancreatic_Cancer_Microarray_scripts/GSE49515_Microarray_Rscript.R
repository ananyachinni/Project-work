#Microarray Data analysis of GSE49515
setwd("C:/gse49515")
celpath = "C:/gse49515"
library(affy)
data = ReadAffy(celfile.path = celpath)   #Reading in the .CEL files
ph = data@phenoData
#ph
ph@data [ ,1] = c("mutant", "mutant", "mutant", "mutant", "mutant", "mutant", "mutant", "mutant", "mutant", 
                  "mutant", "control", "control", "control", "control", "control", "control", "control", 
                  "control", "control", "control")
ph@data
name = "boxplot.jpg"         #Boxplot before normalization
jpeg(name)
boxplot(data,which='pm',col='red',names=ph@data$sample,xlab='Samplename',ylab='logInt')
dev.off()
for (i in 1:20)              #MA plot for the 20 samples before normalization
   {
         name = paste("MAplot",i,".jpg",sep="")
         jpeg(name)
         MAplot(data,which=i)
         dev.off()
         }
data.rma = rma(data) #RMA normalization
data.matrix = exprs(data.rma)
name = "boxplotnorm.jpg"     #Boxplot after normalization
jpeg(name)
boxplot(data.matrix,col='red',names=ph@data$sample,xlab='Samplename',ylab='logInt')
dev.off()
library(oligo)
for (i in 1:20)              #MA plot for the 20 samples after normalization
   {
         name = paste("MAplotnorm",i,".jpg",sep="")
         jpeg(name)
         MAplot(data.rma,which=i)
         dev.off()
         }
ph@data [ ,2] = c("mutant", "mutant", "mutant", "mutant", "mutant", "mutant", "mutant", "mutant", "mutant", 
                  "mutant", "control", "control", "control", "control", "control", "control", "control", 
                  "control", "control", "control")
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
library(biomaRt)
listMarts()
options(digits = 2)
tab = topTable(data.fit.eb,coef=1,number=4000,adjust.method="BH")  #Using topTable for extracting differentially expressed probes
topgenes = tab[tab[, "adj.P.Val"] < 0.05, ]
topups = topgenes[topgenes[, "logFC"] > 1, ]
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
IDs.up = rownames(topups)
IDs.down = rownames(topdowns)
write.table(IDs.up,row.names=FALSE,col.names=FALSE,quote=FALSE,file="C:/gse49515/upIDs.txt")
write.table(IDs.down,row.names=FALSE,col.names=FALSE,quote=FALSE,file="C:/gse49515/downIDs.txt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
affyids = c("228282_at" ,"222061_at" ,"213853_at" ,"218138_at" ,"228925_at" ,"227351_at" ,"225509_at" ,"239698_at" ,"242818_x_at" ,"209585_s_at" ,"235803_at" ,"227335_at" ,"235151_at" ,"227361_at" ,"235728_at"
            ,"219074_at" ,"226230_at" ,"227693_at" ,"239944_at" ,"228619_x_at" ,"232020_at" ,"212905_at" ,"227288_at" ,"228469_at" ,"222500_at" ,"238695_s_at" ,"216392_s_at" ,"226868_at" ,"226826_at" ,"236401_at"
            ,"226040_at" ,"208891_at" ,"235085_at" ,"51146_at" ,"238653_at" ,"223268_at" ,"229366_at" ,"212714_at" ,"236314_at" ,"222867_s_at" ,"220235_s_at" ,"1558250_s_at" ,"203102_s_at" ,"204859_s_at" ,"223412_at"
            ,"242590_at" ,"227247_at" ,"242981_at" ,"243407_at" ,"224959_at" ,"212367_at" ,"208893_s_at" ,"201089_at" ,"225558_at" ,"230235_at" ,"240824_at" ,"206978_at" ,"219777_at" ,"226423_at" ,"1554108_at" ,"227029_at"
            ,"1552318_at" ,"228157_at" ,"218882_s_at" ,"202925_s_at" ,"224881_at" ,"226977_at" ,"230012_at" ,"229794_at" ,"224748_at" ,"218361_at" ,"244786_at" ,"230300_at" ,"209288_s_at" ,"223062_s_at" ,"230836_at"
            ,"218783_at" ,"215985_at" ,"219353_at" ,"227164_at" ,"203143_s_at" ,"227385_at" ,"240673_at" ,"239138_at" ,"206314_at" ,"202666_s_at" ,"229560_at" ,"228650_at" ,"215318_at" ,"212215_at" ,"236583_at" ,"221786_at"
            ,"225290_at" ,"231853_at" ,"229367_s_at" ,"208892_s_at" ,"228453_at" ,"212731_at" ,"227426_at" ,"200799_at" ,"206765_at" ,"234983_at" ,"228336_at" ,"209175_at" ,"203144_s_at" ,"227792_at" ,"235304_at"
            ,"235940_at" ,"218732_at" ,"229285_at" ,"222691_at" ,"1552789_at" ,"223583_at" ,"200703_at" ,"224789_at" ,"229167_at" ,"223351_at" ,"205202_at" ,"222566_at" ,"213626_at" ,"201735_s_at" ,"242317_at" ,"231919_at"
            ,"225256_at" ,"224722_at" ,"225537_at" ,"222589_at" ,"228920_at" ,"243303_at" ,"230860_at" ,"202321_at" ,"225669_at" ,"227040_at" ,"227184_at" ,"201098_at" ,"231784_s_at" ,"204700_x_at" ,"228306_at" ,"227980_at"
            ,"204523_at" ,"203765_at" ,"213552_at" ,"1552330_at" ,"228373_at" ,"220926_s_at" ,"1562608_at" ,"229295_at" ,"225255_at" ,"218303_x_at" ,"229315_at" ,"201824_at" ,"228234_at" ,"1558700_s_at" ,"201193_at"
            ,"238122_at" ,"225086_at" ,"226041_at" ,"228694_at" ,"228243_at" ,"239024_at" ,"218999_at" ,"218805_at" ,"202451_at" ,"213365_at" ,"228071_at" ,"228869_at" ,"242140_at" ,"207125_at" ,"202710_at" ,"230261_at"
            ,"220704_at" ,"228963_at" ,"223809_at" ,"1555037_a_at" ,"204807_at" ,"242648_at" ,"240948_at" ,"230292_at" ,"233329_s_at" ,"1554455_at" ,"235215_at" ,"224827_at" ,"225989_at" ,"235306_at" ,"203614_at"
            ,"1556090_at" ,"213596_at" ,"235478_at" ,"243009_at" ,"205898_at" ,"212565_at" ,"225348_at" ,"220391_at" ,"219003_s_at" ,"226483_at" ,"219487_at" ,"226109_at" ,"238992_at" ,"228087_at" ,"238649_at" ,"242759_at"
            ,"239376_at" ,"204075_s_at" ,"217216_x_at" ,"243405_at" ,"212378_at" ,"222843_at" ,"1563509_at" ,"204172_at" ,"226479_at" ,"1557238_s_at" ,"201503_at" ,"223403_s_at" ,"201964_at" ,"202540_s_at" ,"222980_at"
            ,"212872_s_at" ,"218948_at" ,"240064_at" ,"218979_at" ,"228280_at" ,"236194_at" ,"1556821_x_at" ,"206082_at" ,"213233_s_at" ,"204838_s_at" ,"226458_at" ,"243378_at" ,"229614_at" ,"224953_at" ,"1557487_at"
            ,"223044_at" ,"235675_at" ,"237759_at" ,"1562612_at" ,"230098_at" ,"212407_at" ,"234986_at" ,"236128_at" ,"242953_at" ,"52285_f_at" ,"230252_at" ,"244450_at" ,"212199_at" ,"1552790_a_at" ,"224414_s_at"
            ,"225994_at" ,"226159_at" ,"242943_at" ,"225987_at" ,"205930_at" ,"230009_at" ,"204204_at" ,"206188_at" ,"232180_at" ,"212405_s_at" ,"221211_s_at" ,"227916_x_at" ,"235424_at" ,"222754_at" ,"218643_s_at"
            ,"224840_at" ,"243166_at" ,"219421_at" ,"223433_at" ,"64064_at" ,"222077_s_at" ,"219146_at" ,"218689_at" ,"231323_at" ,"222805_at" ,"232555_at" ,"220146_at" ,"222714_s_at" ,"235651_at" ,"1558801_at" ,"227618_at"
            ,"207394_at" ,"214525_x_at" ,"230399_at" ,"225847_at" ,"209662_at" ,"1554007_at" ,"232068_s_at" ,"230821_at" ,"213817_at" ,"206110_at" ,"229298_at" ,"204085_s_at" ,"207338_s_at" ,"219128_at" ,"227636_at"
            ,"218513_at" ,"242827_x_at" ,"222286_at" ,"218195_at" ,"209301_at" ,"239083_at" ,"201703_s_at" ,"202706_s_at" ,"224341_x_at" ,"204055_s_at" ,"227689_at" ,"219283_at" ,"204256_at" ,"213005_s_at" ,"207008_at"
            ,"227312_at" ,"207513_s_at" ,"222985_at" ,"217542_at" ,"214193_s_at" ,"222585_x_at" ,"206734_at" ,"218519_at" ,"223490_s_at" ,"227187_at" ,"238004_at" ,"204544_at" ,"235286_at" ,"235220_at" ,"223060_at"
            ,"225765_at" ,"209828_s_at" ,"226371_at" ,"236520_at" ,"223434_at" ,"1554456_a_at" ,"1552316_a_at" ,"225352_at" ,"233461_x_at" ,"232725_s_at" ,"203427_at" ,"230174_at" ,"219243_at" ,"225774_at" ,"225769_at"
            ,"215165_x_at" ,"201660_at" ,"243417_at" ,"219017_at" ,"238462_at" ,"218701_at" ,"222811_at" ,"225095_at" ,"234148_at" ,"243003_at" ,"221556_at" ,"229333_at" ,"204125_at" ,"229455_at" ,"203159_at" ,"225661_at"
            ,"217892_s_at" ,"225734_at" ,"236313_at" ,"238937_at" ,"223404_s_at" ,"228167_at" ,"219130_at" ,"223507_at" ,"225191_at" ,"222402_at" ,"213506_at" ,"202872_at" ,"230211_at" ,"223275_at" ,"234295_at" ,"230192_at"
            ,"219123_at" ,"1569189_at" ,"235625_at" ,"227027_at" ,"1563259_at" ,"214440_at" ,"223386_at" ,"1560800_at" ,"207643_s_at" ,"235767_x_at" ,"1556185_a_at" ,"217403_s_at" ,"200800_s_at" ,"237201_at" ,"232883_at"
            ,"227813_at" ,"240057_at" ,"204084_s_at" ,"1552553_a_at" ,"227536_at" ,"210772_at" ,"236249_at" ,"228190_at" ,"228561_at" ,"230383_x_at" ,"1566608_at" ,"233208_x_at" ,"206983_at" ,"225580_at" ,"235885_at"
            ,"214155_s_at" ,"224856_at" ,"1554345_a_at" ,"210176_at" ,"238902_at" ,"1562250_at" ,"235117_at" ,"222869_s_at" ,"224365_s_at" ,"219294_at" ,"220005_at" ,"205098_at" ,"230226_s_at" ,"204959_at" ,"219913_s_at"
            ,"213625_at" ,"220969_s_at" ,"227600_at" ,"222872_x_at" ,"204103_at" ,"205140_at" ,"223588_at" ,"240481_at" ,"232504_at" ,"219467_at" ,"235443_at" ,"229934_at" ,"227559_at" ,"209906_at" ,"219383_at" ,"216069_at"
            ,"214414_x_at" ,"214329_x_at" ,"213218_at" ,"223396_at" ,"227626_at" ,"1552315_at" ,"215071_s_at" ,"232623_at" ,"229970_at" ,"1556820_a_at" ,"202949_s_at" ,"206026_s_at" ,"206698_at" ,"219684_at" ,"241721_at"
            ,"226481_at" ,"221027_s_at" ,"233085_s_at" ,"214130_s_at" ,"232628_at" ,"206991_s_at" ,"228387_at" ,"217832_at")
results1 <- getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene', 'hgnc_symbol'), filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)      #Mapping upregulated DE probes to HGNC symbols
write.csv(results1, file="ups.csv")
affyids = c("234055_s_at" ,"225262_at" ,"230802_at" ,"218708_at" ,"204141_at" ,"200989_at" ,"218631_at" ,"201465_s_at" ,"237643_at" ,"225884_s_at" ,"214211_at" ,"1553338_at" ,"202499_s_at" ,"38037_at",
            "226140_s_at" ,"237718_at" ,"211458_s_at" ,"233899_x_at" ,"203821_at" ,"211924_s_at" ,"210512_s_at" ,"1553133_at" ,"230170_at" ,"238736_at" ,"234907_x_at" ,"225699_at" ,"206877_at" ,"211199_s_at" ,"209545_s_at",
            "201416_at" ,"209305_s_at" ,"201939_at" ,"208536_s_at" ,"230348_at" ,"231889_at" ,"209383_at" ,"222088_s_at" ,"219312_s_at" ,"237485_at" ,"217738_at" ,"235739_at" ,"208078_s_at" ,"211998_at" ,"216260_at",
            "1564274_at" ,"208868_s_at" ,"202340_x_at" ,"208869_s_at" ,"1564430_at" ,"201466_s_at" ,"222343_at" ,"237510_at" ,"238893_at" ,"1568768_s_at" ,"218880_at" ,"210845_s_at" ,"230511_at" ,"211324_s_at",
            "201464_x_at" ,"210634_at" ,"1558691_a_at" ,"212171_x_at" ,"203574_at" ,"228562_at" ,"241824_at" ,"1554291_at" ,"232045_at" ,"235102_x_at" ,"213844_at" ,"204621_s_at" ,"209398_at" ,"217739_s_at" ,"243918_at",
            "219624_at" ,"223169_s_at" ,"1556806_at" ,"216248_s_at" ,"237107_at" ,"205330_at" ,"200731_s_at" ,"204622_x_at" ,"243538_at" ,"1564972_x_at" ,"222044_at" ,"243605_at" ,"213281_at" ,"1554929_at" ,"1569263_at",
            "243798_at" ,"204285_s_at" ,"209020_at" ,"244257_at" ,"229106_at" ,"239124_at" ,"224453_s_at" ,"214290_s_at" ,"223584_s_at" ,"1557166_at" ,"1554309_at" ,"210232_at" ,"218280_x_at" ,"1560058_at" ,"230134_s_at",
            "1556607_at" ,"231182_at" ,"236213_at" ,"224836_at" ,"224164_at" ,"1553861_at" ,"214349_at" ,"1554229_at" ,"239845_at" ,"213668_s_at" ,"231165_at" ,"241985_at" ,"1560846_at" ,"241938_at" ,"1564970_at",
            "1556608_a_at" ,"230048_at" ,"1553134_s_at" ,"1557459_at" ,"219622_at" ,"242975_s_at" ,"214230_at" ,"207574_s_at" ,"224454_at" ,"237009_at" ,"231630_at" ,"205548_s_at" ,"227501_at" ,"230133_at" ,"222045_s_at",
            "232632_at" ,"243296_at" ,"244840_x_at" ,"1558331_at" ,"226342_at" ,"1559060_a_at" ,"239404_at" ,"1554676_at" ,"1559582_at" ,"1554089_s_at" ,"213537_at" ,"227757_at" ,"211527_x_at" ,"241027_at" ,"1569864_at",
            "1554906_a_at" ,"216027_at" ,"238389_s_at" ,"219228_at" ,"244753_at" ,"226608_at" ,"237464_at" ,"244047_at" ,"202464_s_at" ,"214696_at" ,"238509_at" ,"216236_s_at" ,"1563113_at" ,"212374_at" ,"204286_s_at",
            "230156_x_at" ,"237246_at" ,"222309_at" ,"208553_at" ,"200730_s_at" ,"228097_at" ,"244546_at" ,"1555167_s_at" ,"222142_at" , "213524_s_at")
results2 <- getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene', 'hgnc_symbol'), filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)     #Mapping downregulated DE probes to HGNC symbols
write.csv(results2, file="downs.csv")