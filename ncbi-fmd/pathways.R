library(tximport)
library(rjson)
library(readr)
library("edgeR")
library(biomaRt)
library(dplyr) 
library("AnnotationDbi")
library("org.Hs.eg.db")
library(magrittr)
library(pathview)
library(gage)
library(gageData)

samples <- read.table(file.path("/home/sherine/work/rnaseq/kaz", "samples.txt"), header = FALSE)
files <- file.path("/home/sherine/work/rnaseq", "kaz", samples$V1, "quant.sf")
print ('Files are :') 
print (files)
names(files) <- paste0("sample", 1:20)
tx2gene <- read.csv(file.path("/home/sherine/work/rnaseq/kaz", "NCBItx2gene.csv"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,ignoreTxVersion = TRUE)
write.csv(txi$counts, "counts.csv")

cts <- txi$counts
group = factor(c(rep("Control", 10), rep("Tumor",10)) )
dge = DGEList(counts=cts, genes= rownames(data), group=group)
countsPerMillion <- cpm(dge)
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
#-----------------------------------Pathway-Analysis---------------------------------------
write.csv(etp$table, "dge-NCBI.csv")
res = etp$table 
write.csv(res, "res.csv")
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ALIAS",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ALIAS",
                     multiVals="first")
write.csv(res, "res.csv")

#data(kegg.sets.hs)
#data(sigmet.idx.hs)

#kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
#head(kegg.sets.hs, 3)
foldchanges = res$logFC
names(foldchanges) = res$entrez
write.csv((names(foldchanges)), "namefolds.csv")
#---------------------------------------------------Use Kegg and gage to get upregulated and downregulated pathways
data(kegg.gs)
keggres = gage(foldchanges, gsets =kegg.gs, same.dir =TRUE, compare="unpaired")
#keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE,compare="unpaired")
lapply(keggres, head)
write.csv(keggres,"keggres.csv")
keggrespathwaysup = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=5) %>%
  .$id %>%
  as.character()
keggrespathwaysdn = data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=5) %>%
  .$id %>%
  as.character()
write.csv(keggrespathwaysup, "keggspathsup.csv")
write.csv(keggrespathwaysdn, "keggspathsdn.csv")
#--------------------------------------------------------------------------------------------------------------------------------
keggresidsup = substr(keggrespathwaysup, start=1, stop=8)
keggresidsup
keggresidsdn = substr(keggrespathwaysdn, start=1, stop=8)
#data(go.sets.hs)
#data(go.subs.hs)
#gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=kegg.gs, same.dir =TRUE, compare ="unpaired")
#gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE, compare ="unpaired")
lapply(gobpres, head)
#---------------------------------------------------------Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa", new.signature=FALSE)

#---------------------------------------------------------plot multiple pathways ups and downs 
tmpup = sapply(keggresidsup, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa"))
tmpdn = sapply(keggresidsdn, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa"))
#-----------------------------------------------------------Gene Ontology----------------------------------------------------------------
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE, compare = "unpaired", cutoff = 0.01)
lapply(gobpres, head)
#gog = head(gobpres$greater[,1:5], 10)
#gol = head(gobpres$less[,1:5],50)  
#write.csv(gog, "gog.csv")
#write.csv(gol, "gol.csv")

sig<-sigGeneSet(gobpres, outname="gse16873.kegg")
str(sig, strict.width='wrap')
write.csv(sig$greater, "sig.greater.csv") 
write.csv(sig$less, "sig.less.csv") 
write.table(rbind(sig$greater, sig$less), file = "sig.txt", sep ="\t") 
#-----------------------------------------------------------Print WArnings---------------------------------------------------------------
warnings()
