#!/usr/bin/env Rscript

rm(list=objects())
graphics.off()

setwd("/home/ubuntu/Documents/AMI2B/Hackathon/Projet/R")

counts <- read.table("counts.txt",header = T, row.names = 1)
colnames(counts) <- sub("\\.bam", "", colnames(counts))
countData <- counts[rowSums(counts[,6:13]) > 0,][,6:13]

metadata <- read.table("SraRunTable.txt",header = T, row.names = 1, sep=",")

id <- rownames(metadata)

tumor <- metadata$tumor_class
tumor <- replace(tumor, tumor=="Class 1", "Class1")
tumor <- replace(tumor, tumor=="", "Class2")

mutation <- metadata$sf3b1_mutation_status
mutation <- replace(mutation, mutation=="Class 1 R625C", "Yes")
mutation <- replace(mutation, mutation=="R625C", "Yes")
mutation <- replace(mutation, mutation=="R625H", "Yes")
mutation <- replace(mutation, mutation=="WT", "No")
mutation <- replace(mutation, mutation=="Class 1 WT", "No")

metaData <- data.frame(id, tumor, mutation)

metaData$tumor <- factor(metaData$tumor)
metaData$mutation <- factor(metaData$mutation)

library(DESeq2)

# des <- ~ tumor + mutation + mutation:tumor
des <- ~ mutation

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData, 
                              design = des, tidy = FALSE)

dds <- DESeq(dds)

# estimateSizeFactors => computes the relative library depth of each sample 
# estimateDispersions => estimates the dispersion of counts for each gene
# nbinomWaldTest => computes the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs

res <- results(dds, tidy=FALSE)
summary(res)

### P-values

par(mfrow=c(1,2))

library(ggplot2)
library(hrbrthemes)

pvalue.df <- data.frame(type = c(rep("pvalue", nrow(res)), rep("pvalue adjusted", nrow(res))),
                        pvalue = c(res$pvalue, res$padj))

ggplot(pvalue.df, aes(x=pvalue, fill=type)) + 
  geom_histogram(aes(y=..density..), binwidth=0.05, color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic(base_size = 11) +
  labs(fill="")

plot1 <- ggplot(data.frame(res), aes(x=pvalue)) +
  geom_histogram(aes(y=..density..), binwidth=0.07, fill="#69b3a2", color="#e9ecef", alpha=0.6, position = 'identity') +
  theme_classic(base_size = 11)

plot2 <- ggplot(data.frame(res), aes(x=padj)) +
  geom_histogram(aes(y=..density..), binwidth=0.07, fill="#404080", color="#e9ecef", alpha=0.6, position = 'identity') +
  theme_classic(base_size = 11)

require(gridExtra)
grid.arrange(plot1, plot2, ncol=2)

table(res$pvalue < 0.05)
table(res$padj < 0.05)

### MA plot ###

par(mfrow=c(1,1))

plotMA(res, ylim=c(-20,20))

### Count plot ###

res <- res[order(res$padj),] # Order by adjusted p-value
head(res)

par(mfrow=c(1,3))

plotCounts(dds, gene="ENSG00000245694", intgroup="mutation") # CRNDE
plotCounts(dds, gene="ENSG00000101019", intgroup="mutation") # UQCC
plotCounts(dds, gene="ENSG00000114770", intgroup="mutation") # ABCC5

par(mfrow=c(2,3))

for (i in 1:6){
  plotCounts(dds, gene=rownames(res)[i], intgroup="mutation")
}

# Next steps in exploring these data...BLAST to database to find associated gene function

### Volcano plot ###

library(ggrepel)

res <- results(dds, tidy=TRUE)
res.condition1 <- subset(res, padj<.1)
res.condition2 <- subset(res, (padj<.005 & -log10(pvalue) > 6))

ggplot(data = res.condition2, aes(x = log2FoldChange, y = -log10(pvalue), label = row)) +
  geom_point(data = res, aes(x = log2FoldChange, y = -log10(pvalue)), color = "black") +
  geom_point(data = res.condition1, aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue") +
  geom_point(color = "red") +
  geom_text_repel(size=3) +
  xlim(-7.3, 7.3) +
  theme_classic(base_size = 11)

### PCA ###

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="mutation")
