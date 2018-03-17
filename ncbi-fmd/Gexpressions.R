library(tximport)
library(rjson)
library(readr)
library("edgeR")
library(biomaRt)
library(org.Hs.eg.db) 
library(gplots)


samples <- read.table(file.path("/data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/ncbi-fmd", "samples.txt"), header = FALSE)

files <- file.path("/data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz", "ncbi-fmd", samples$V1, "quant.sf")
#print ('Files are :') 
#print (files)
names(files) <- paste0("sample", 1:20)
tx2gene <- read.csv(file.path("/data/WHRI-EndocrinePituitaryGroup/sherine/rnaseq/rnaseq-kaz/ncbi-fmd", "NCBItx2gene.csv"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,ignoreTxVersion = TRUE)
write.csv(txi$counts, "ControlTumorcounts.csv")

cts <- txi$counts
group = factor(c(rep("Control", 10), rep("Tumor",10)) )
dge = DGEList(counts=cts, genes= rownames(cts), group=group)
countsPerMillion <- cpm(dge) #, prior.count=2, log=TRUE)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1 
summary(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dge <- dge[keep,]
summary(cpm(dge))
dge <- calcNormFactors(dge, method="TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
dge <- estimateTrendedDisp(dge)
et <- exactTest(dge, pair=c("Control", "Tumor"))
etp <- topTags(et, n= 100000, adjust.method="BH", sort.by="PValue", p.value = 1) 
#--------------------------------------------Some Plots-------------------------------------
#labels=c("C1", "C2", "C3", "C4", "C5", "C6","C7", "C8", "C9", "C10", "T1","T2","T3","T4","T5","T6","T7","T8","T9","T10")
labels = c("c08","c10","c14","c15","c16","c17","c19","c20","c21","c22","t70","t66","t52","t76","t77","t79","t81","t83","t85","t87")
pdf("ControlTumorMA.pdf")
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR < .2, "black", "red" ) )
dev.off()
pdf("ControlTumorMDS.pdf")
plotMDS(dge, labels = labels)
dev.off()
pdf("ControlTumorbcv.pdf") 
plotBCV(dge)
dev.off()
#------------------------------------------------------Heatmap-----------------------------------------
pdf("ControlTumorheatmap.pdf")
logCPM = countsPerMillion
o = rownames(etp$table[abs(etp$table$logFC)>2 & etp$table$PValue<0.01, ])
logCPM <- logCPM[o[1:100],]
colnames(logCPM) = labels
logCPM <- t(scale(t(logCPM)))
write.csv(logCPM, "ControlTumorCPM.csv")
require("RColorBrewer") 
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-3, 3, length.out=101)
#heatmap.2(logCPM, col=myCol, breaks=myBreaks, Rowv=TRUE,Colv=TRUE,  main="Controls vs Tumors Heatmap", key=T, keysize=0.6,scale="none",trace="none", dendrogram="both", cexRow=0.2, cexCol=0.9, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),  reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2"))

heatmap.2(logCPM, col=myCol, breaks=myBreaks, Rowv=TRUE,Colv=TRUE, main="Controls vs Tumors Heatmap", key=T, keysize=0.7,scale="none",trace="none", dendrogram="both", cexRow=0.2, cexCol=0.9, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),  distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()
#---------------------------------------------------------------------------------
write.csv(etp$table, "dgeControlTumor.csv")
