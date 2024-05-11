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

  
  # DEG volcano
  print(OE_name)
  # Filter the columns of log2FC and value2 to see if they are not numbers
  deg.data <- annotation
  deg.data$log10FDR <- -log10(deg.data$pvalue)
  
  # add a new row named Group
  deg.data$Group = "normal"
  # Genes with padj less than 0.05 and logFC less than 2 were set as significantly down-regulated genes
  # Adjust the value according to your own data
  deg.data$Group[which( (deg.data$pvalue < 0.05) & (deg.data$log2FoldChange > 2) )] = "up"
  deg.data$Group[which( (deg.data$pvalue < 0.05) & (deg.data$log2FoldChange < -2) )] = "down"
  # View the number of up-regulated and down-regulated genes
  table(deg.data$Group)
  
  # add a new row named Label
  deg.data$Label = ""
  # Sort the FDR values of differentially expressed genes from small to large
  deg.data <- deg.data[order(deg.data$padj), ]
  # Among the highly expressed genes, select the 10 with the smallest adj.P.Val
  # up.genes <- head(deg.data$Row.names[which(deg.data$Group == "up")], 10)
  up.genes <- head(deg.data$Gene.Name[which(deg.data$Group == "up")], 10)
  up.genes <- up.genes[nzchar(up.genes)]
  # Low expressed adj.P.Val minimum 10
  # down.genes <- head(deg.data$Row.names[which(deg.data$Group == "down")], 10)
  down.genes <- head(deg.data$Gene.Name[which(deg.data$Group == "down")], 10)
  down.genes <- down.genes[nzchar(down.genes)]
  # Merge up.genes and down.genes and add them to Label
  deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
  deg.data$Label[match(deg.top10.genes, deg.data$Gene.Name)] <- deg.top10.genes
  # Adding the top ten significantly differentially expressed genes to the volcano plot
  p <- ggscatter(deg.data, x = "log2FoldChange", y = "log10FDR",
            color = "Group", 
            palette = c("#00AFBB", "#999999", "#FC4E07"),
            size = 2,
            label = deg.data$Label, 
            font.label = 8, 
            repel = T) + 
    theme_base() + 
    xlab("log2(FC)") + ylab("-log10(adujsted_P)") + 
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    geom_vline(xintercept = c(-1.5,1.5), linetype="dashed") +
    labs(title = OE_name)
  
  
  # Output the picture 
  pdf(file=paste0("/Users/lubeifang/Desktop/sRNA/Figure/Volcano/",OE_name,"_vol.pdf"), bg="transparent")
  print(p)
  dev.off()

  
  # GO
  updeg = subset(annotation, log2FoldChange > 1)
  GO <- enricher(updeg$Row.names,
                TERM2GENE = PAgo_term2gene,
                TERM2NAME = PAgo_term2name,
                pvalueCutoff = 0.05,)
  pdf(file=paste0("/Users/lubeifang/Desktop/sRNA/Figure/GO/",OE_name,"_GOup.pdf"), bg="transparent")
  p <- dotplot(GO) + labs(title = paste0(OE_name,"up"))
  print(p)
  dev.off()
  
  KEGG <- enricher(updeg$Row.names,
                 TERM2GENE = PAkegg_term2gene,
                 TERM2NAME = PAkegg_term2name,
                 pvalueCutoff = 0.05,)
  pdf(file=paste0("/Users/lubeifang/Desktop/sRNA/Figure/KEGG/",OE_name,"_KEGGup.pdf"), bg="transparent")
  p <- barplot(KEGG) + labs(title = paste0(OE_name,"up"))
  print(p)
  dev.off()

  downdeg = subset(annotation, log2FoldChange < -1)
  GO <- enricher(downdeg$Row.names,
                 TERM2GENE = PAgo_term2gene,
                 TERM2NAME = PAgo_term2name,
                 pvalueCutoff = 0.05,)
  pdf(file=paste0("/Users/lubeifang/Desktop/sRNA/Figure/GO/",OE_name,"_GOdown.pdf"), bg="transparent")
  p <- dotplot(GO) + labs(title = paste0(OE_name,"down"))
  print(p)
  dev.off()
  
  KEGG <- enricher(downdeg$Row.names,
                   TERM2GENE = PAkegg_term2gene,
                   TERM2NAME = PAkegg_term2name,
                   pvalueCutoff = 0.05,)
  pdf(file=paste0("/Users/lubeifang/Desktop/sRNA/Figure/KEGG/",OE_name,"_KEGGdown.pdf"), bg="transparent")
  p <- barplot(KEGG) + labs(title = paste0(OE_name,"down"))
  print(p)
  dev.off()
}


# build edge file of DEG
edge <- data.frame()
print(OE_name.list)

for (OE_name in OE_name.list) {
  print(OE_name)
  deg=read.csv(paste0("/Users/lubeifang/Desktop/sRNA/PA/srnadeg/",OE_name,"_DEG.csv"))
  deg<-subset(deg, padj < 0.05)
  deg<-subset(deg,log2FoldChange > 1 | log2FoldChange < -1)
  sRNA=rep(OE_name, times = nrow(deg))
  temp_edge <- data.frame(data.frame(fromNode = sRNA), toNode = deg$Row.names, weight = deg$log2FoldChange)
  edge <- rbind(edge,temp_edge)
}
a=unique(edge$toNode)
edge$relation <- ifelse(edge$weight>0,1,0)

# remove sRNA in whole edge figure
edge <- edge[!grepl("sRNA", edge$toNode),]
write.csv(edge,"/Users/lubeifang/Desktop/sRNA/PA/sRNAedge.csv",row.names = F)

know44rna <- read.csv('/Users/lubeifang/Desktop/sRNA/sRNAprediction/ncRNA-PAO1.csv', sep = ';', header = F)
know44rna$V3 <- substr(know44rna$V3,7,14)
srna.df <- merge(know44rna,edge,by.x = "V3", by.y = "toNode")
srna.df <- srna.df[,c(6,1,7,8)]
na.omit(srna.df)

remain = data.frame(name=setdiff(know44rna$V3,srna.df$V3))
write.csv(srna.df,"/Users/lubeifang/Desktop/sRNA/PA/srna_inter2.txt", quote = F, row.names = F)

# 计算每个值在列中出现的频次
value_counts <- table(edge$toNode)
plot(value_counts)
# 筛选值
occurrence_values <- names(value_counts[value_counts > 2])

# GO and KEGG of most regular genes ( more than 2 times)
GO <- enricher(occurrence_values,
               TERM2GENE = PAgo_term2gene,
               TERM2NAME = PAgo_term2name,
               pvalueCutoff = 0.05,)
pdf(file=paste0("/Users/lubeifang/Desktop/sRNA/Figure/commongenes_GO.pdf"), bg="transparent")
p <- dotplot(GO) + labs(title = "GO analysis of enriched DEGs")
print(p)
dev.off()
KEGG <- enricher(occurrence_values,
                 TERM2GENE = PAkegg_term2gene,
                 TERM2NAME = PAkegg_term2name,
                 pvalueCutoff = 0.05,)
pdf(file=paste0("/Users/lubeifang/Desktop/sRNA/Figure/commongenes_KEGG.pdf"), bg="transparent")
p <- barplot(KEGG) + labs(title = "KEGG analysis of enriched DEGs")
print(p)
dev.off()

# 使用逻辑条件筛选需要保留的行
# !（取反）和 %in%（判断是否包含）
filtered_df <- edge[!(edge$toNode %in% occurrence_values), ]
filtered_df <- filtered_df[!grepl("sRNA", filtered_df$toNode), ]
unique(filtered_df$fromNode)
a=table(filtered_df$toNode)
occurrence_values <- names(a[a > 5])
# 打印筛选后的表格
print(filtered_df)
write.csv(filtered_df,"sRNAfiltered.csv",row.names = F)


#PCA
head(exprSet)
colnames(exprSet)
pattern <- "(\\d{4}).*" # 正则表达式模式，提取四位数字
extracted <- sub(pattern, "\\1", colnames(exprSet))
print(extracted)

condition <- factor(extracted)
coldata <- data.frame(row.names = colnames(OE_name), condition)
head(coldata)

dds <- DESeqDataSetFromMatrix(countData=OE_name, colData=coldata,
                              design=~condition)
dim(dds)
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata)



## test the intersection between Combat, Combat-seq and raw counts.
library(VennDiagram)
setwd("/Users/lubeifang/Desktop/")
a = read.csv("/Users/lubeifang/Desktop/sRNA/PA/PA_oldWT/sPAF0807_DEG.csv", header = T)
b = read.csv("/Users/lubeifang/Desktop/sRNA/PA/batchremoved_combatseq/sPAF0807_DEG.csv", header = T)
c = read.csv("/Users/lubeifang/Desktop/sRNA/PA/removebatch_combat/sPAF0807_DEG.csv", header = T)
venn.diagram(
  x = list(a$Row.names, b$Row.names, c$Row.names),
  category.names = c("Raw" , "Combat-seq" , "Combat"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)

edge <- read.csv("/Users/lubeifang/Desktop/sRNA/PA/edge.txt",sep = "\t")
frequency_table <- as.data.frame(table(edge[, 1], edge[, 4]))

ggplot(frequency_table[frequency_table$Var2 == "1",]) + 
  geom_hline(yintercept=mean(frequency_table[frequency_table$Var2 == "1",]$Freq),linetype=2,size=.25,colour="grey") + 
  geom_bar(aes(x=Var1,y=Freq),stat="identity",fill="#FDCF9E",colour=NA) + 
#  ylim(0,150) + 
  coord_flip() +
  theme_classic() +
  geom_text(aes(x=Var1, y=Freq + 7, label=Freq)) +
  xlab("") +
  ggtitle("up")
ggsave("/Users/lubeifang/Desktop/sRNA/Figure/Figure2/updegcounts.pdf")

ggplot(frequency_table[frequency_table$Var2 == "0",]) + 
  geom_hline(yintercept=mean(frequency_table[frequency_table$Var2 == "0",]$Freq),linetype=2,size=.25,colour="grey") + 
  geom_bar(aes(x=Var1,y=Freq),stat="identity",fill="#B8A8CF",colour=NA) + 
#  ylim(0,200) + 
  coord_flip() +
  theme_classic() +
  geom_text(aes(x=Var1, y=Freq + 7, label=Freq)) +
  xlab("") +
  scale_y_reverse() +
  ggtitle("down") +
  geom_hline(yintercept = -8) 
ggsave("/Users/lubeifang/Desktop/sRNA/Figure/Figure2/downdegcounts.pdf")


library(gridExtra)
grid.arrange(up, down, nrow = 1)



  