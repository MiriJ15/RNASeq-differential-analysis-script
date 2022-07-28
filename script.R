#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
BiocManager::install("ComplexHeatmap")

library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library("org.Hs.eg.db")
setwd("/home/yeshua/valm/Miri/RNAseq/results_matru")
getwd()

Count<-read.delim("count_matru.csv", header= TRUE, row.names=1, sep = ",")
#filter in order to minimize noise (0's)
Counts<- Count[which(rowSums(Counts)>50),]
condition<- factor(c("AloneMedia","AloneMedia","AloneMedia","AloneMedia","AloneMedia","AloneMedia", "AloneSaliva","AloneSaliva","AloneSaliva","AloneSaliva","AloneSaliva","AloneSaliva", "Coagg","Coagg","Coagg","Coagg","Coagg","Coagg"))
coldata<- data.frame(row.names = colnames(Counts), condition)
dds<- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)
dds<- DESeq(dds)
vsdata<- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="condition")
plotDispEsts(dds)
res1 <- results(dds, contrast = c("condition", "Coagg", "AloneMedia"))
res2 <- results(dds, contrast = c("condition", "AloneSaliva", "AloneMedia"))
res3 <- results(dds, contrast = c("condition", "Coagg", "AloneSaliva"))
sigs1<-na.omit(res1)
sigs1<- sigs1[sigs1$padj <0.05,]
sigs2<-na.omit(res2)
sigs2<- sigs2[sigs2$padj <0.05,]
sigs3<-na.omit(res3)
sigs3<- sigs3[sigs3$padj <0.05,]
sigs1
sigs2
sigs3
write.csv(sigs1, file ="CoaggVsAloneMedia.csv")
write.csv(sigs2, file ="AloneSalivaVsAloneMedia.csv")
write.csv(sigs3, file ="CoaggVsAloneSaliva.csv")

sigs1.df<- as.data.frame(sigs1)
sigs2.df<- as.data.frame(sigs2)
sigs3.df<- as.data.frame(sigs3)
sigs1.df$symbol<-mapIds(org.Hs.eg.db, keys= rownames(sigs1.df), keytype = "ENTREZID", column = "SYMBOL")
sigs2.df$symbol<-mapIds(org.Hs.eg.db, keys= rownames(sigs2.df), keytype = "ENSEMBL", column = "SYMBOL")
sigs3.df$symbol<-mapIds(org.Hs.eg.db, keys= rownames(sigs3.df), keytype = "ENSEMBL", column = "SYMBOL")
sigs1.df
sigs2.df
sigs3.df