#!/usr/bin/env Rscript

rm(list=objects())
graphics.off()

library(DESeq2)
library(tidyr)
library(ggplot2)
library(ggrepel)

path_to_figures <- "/home/ubuntu/Documents/AMI2B/Hackathon/Projet/Figures"
path_to_counts <- "/home/ubuntu/Documents/AMI2B/Hackathon/Projet/Data/Counts"
path_to_metadata <- "/home/ubuntu/Documents/AMI2B/Hackathon/Projet/R"

counts <- paste(path_to_counts, "counts.txt", sep="/") %>% read.table(header = T, row.names = 1)
colnames(counts) <- sub("\\.bam", "", colnames(counts))
countData <- counts[rowSums(counts[,6:13]) > 0,][,6:13]

metadata <- paste(path_to_metadata, "SraRunTable.txt", sep="/") %>% read.table(header = T, row.names = 1, sep=",")

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

pvalue.df <- data.frame(type = c(rep("pvalue", nrow(res)), rep("pvalue adjusted", nrow(res))),
                        pvalue = c(res$pvalue, res$padj))

plot.pvalue.bind <- ggplot(pvalue.df, aes(x=pvalue, fill=type)) + 
  geom_histogram(aes(y=..density..), binwidth=0.05, color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic(base_size = 11) +
  labs(fill="")

paste(path_to_figures, "histogram_pvalue_bind.png", sep="/") %>%
  ggsave(plot = plot.pvalue.bind)
rm(plot.pvalue.bind)

plot.pvalue <- ggplot(data.frame(res), aes(x=pvalue)) +
  geom_histogram(aes(y=..density..), binwidth=0.07, fill="#69b3a2", color="#e9ecef", alpha=0.6, position = 'identity') +
  theme_classic(base_size = 11)

plot.padj <- ggplot(data.frame(res), aes(x=padj)) +
  geom_histogram(aes(y=..density..), binwidth=0.07, fill="#404080", color="#e9ecef", alpha=0.6, position = 'identity') +
  theme_classic(base_size = 11)

require(gridExtra)
plot.pvalue.unbind <- grid.arrange(plot.pvalue, plot.padj, ncol=2)

paste(path_to_figures, "histogram_pvalue_unbind.png", sep="/") %>%
  ggsave(plot = plot.pvalue.unbind)
rm(plot.pvalue, plot.padj, plot.pvalue.unbind)

table(res$pvalue < 0.05)
table(res$padj < 0.05)

### MA plot ###

par(mfrow=c(1,1))

paste(path_to_figures, "MA_plot.png", sep="/") %>% png(width = 800, height = 800)
plotMA(res, ylim=c(-7,7))
dev.off()

### Count plot ###

res <- res[order(res$padj),] # Order by adjusted p-value
head(res)

par(mfrow=c(1,3))

plotCounts(dds, gene="ENSG00000245694", intgroup="mutation") # CRNDE
plotCounts(dds, gene="ENSG00000101019", intgroup="mutation") # UQCC
plotCounts(dds, gene="ENSG00000114770", intgroup="mutation") # ABCC5

paste(path_to_figures, "count_plot.png", sep="/") %>% png(width = 800, height = 800)
par(mfrow=c(2,3))
for (i in 1:6){
  plotCounts(dds, gene=rownames(res)[i], intgroup="mutation")
}
dev.off()

### Volcano plot ###

res <- results(dds, tidy=TRUE)
res.condition1 <- subset(res, padj<.1)
res.condition2 <- subset(res, padj<.001)

plot.volcano <- ggplot(data = res.condition2, aes(x = log2FoldChange, y = -log10(pvalue), label = row)) +
  geom_point(data = res, aes(x = log2FoldChange, y = -log10(pvalue)), color = "black") +
  geom_point(data = res.condition1, aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue") +
  geom_point(color = "red") +
  geom_text_repel(size=3) +
  xlim(-7.3, 7.3) +
  theme_classic(base_size = 11)

paste(path_to_figures, "volcano_plot.png", sep="/") %>%
  ggsave(plot = plot.volcano)
rm(plot.volcano)

### PCA ###

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="mutation")
