# BiocManager::install("DESeq2")

# START
rm(list=ls())
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(clusterProfiler)
library(dplyr)

# read original counts from featurecounts
getwd()
setwd("/Users/lubeifang/Desktop/")

counts <- read.csv("/Users/lubeifang/Desktop/sRNA/PA/counts2023.txt",
                   row.names=1,sep="\t", header = T, skip = 1)
# counts <- read.csv("/Users/lubeifang/Desktop/sRNA/PA/srnacounts.txt",
#                    row.names=1,sep="\t", header = T, skip = 1)
# counts$WTEV.2.sorted.bam <- counts$WTEV.1.sorted.bam
# load annotation 
PA_anno = read.csv("/Users/lubeifang/Desktop/BIOTOOLS/ref/PAO1_ref.csv", header = T)

# GO term
PAgo = read.csv("/Users/lubeifang/Desktop/BIOTOOLS/ref/PA_GOterm.csv",header = T)
PAgo_term2gene <- as.data.frame(PAgo[, c(5,1)])
PAgo_term2name <- as.data.frame(PAgo[, c(5,6)])

PAkegg = read.csv("/Users/lubeifang/Desktop/BIOTOOLS/ref/PA_KEGG.csv",header = T)
PAkegg_term2gene <- as.data.frame(PAkegg[, c(4,1)])
PAkegg_term2name <- as.data.frame(PAkegg[, c(4,3)])

# Delete the first 5 cols
exprSet=counts[,6:ncol(counts)]
colnames(exprSet)

OE_name.df <- exprSet[, 1:(ncol(exprSet) - 2)]
OE_name = substr(colnames(OE_name.df),1,8)
OE_name.list = unique(OE_name)
print(OE_name.list)
WT1 <- "WTEV.1.1.sorted.bam"
WT2 <- "WTEV.1.2.sorted.bam"

# generated annotaed DEG files in batches
for (OE_name in OE_name.list) {
  print(OE_name)
  col_name1 <- paste0(OE_name,".1.sorted.bam")
  col_name2 <- paste0(OE_name,".2.sorted.bam")
  OE <- data.frame()
  OE <- cbind(as.data.frame(as.integer(exprSet[[WT1]])),
              as.data.frame(as.integer(exprSet[[WT2]])),
              as.data.frame(as.integer(exprSet[[col_name1]])),
              as.data.frame(as.integer(exprSet[[col_name2]])))
  colnames(OE) <- c("WT1","WT2",paste0(OE_name,".1"),paste0(OE_name,".2"))
  rownames(OE) <- rownames(counts)
  condition <- factor(c(rep("control",2), rep("treat",2)))
  coldata <- data.frame(row.names = colnames(OE), condition)
  head(coldata)
  
  dds <- DESeqDataSetFromMatrix(countData=OE, colData=coldata,
                                design=~condition)
  dim(dds)
  dds <- DESeq(dds)
  vsdata <- vst(dds, blind=FALSE)
  res <- results(dds,contrast = c('condition', 'treat', 'control'))
  res <- res[order(res$padj),]
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,
                                                            normalized=TRUE)),by="row.names",sort=FALSE)
  annotation <- merge(resdata, PA_anno, by.x="Row.names", by.y="Locus.Tag", all.x= T)
  annotation <- annotation[order(annotation$padj),]
  annotation <- subset(annotation, padj < 0.05 & pvalue < 0.01)
  annotation_DEG <- subset(annotation, log2FoldChange > 2 | log2FoldChange < -2)
  annotation_DEG <- annotation_DEG[order(annotation_DEG$log2FoldChange, decreasing = TRUE),]
  write.csv(annotation_DEG,file = paste0("/Users/lubeifang/Desktop/sRNA/PA/srnadeg/",OE_name,"_DEG.csv"),row.names = FALSE)
}
  
