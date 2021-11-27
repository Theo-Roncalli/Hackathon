#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

install.packages("pacman", repos="https://cran.irsn.fr/")
library(pacman)

p_load(tidyr, ggplot2, ggrepel, BiocManager)

BiocManager::install("DESeq2")

library(DESeq2)

path_to_counts <- args[1]
path_to_metadata <- args[2]

counts <- read.table(path_to_counts, header = T, row.names = 1)
colnames(counts) <- sub("\\.bam", "", colnames(counts))
countData <- counts[rowSums(counts[,6:13]) > 0,][,6:13]

metadata <- read.table(path_to_metadata, header = T, row.names = 1, sep=",")

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

des <- ~ mutation

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData, 
                              design = des, tidy = FALSE)

dds <- DESeq(dds)

res <- results(dds, tidy=FALSE)

### P-values

pvalue.df <- data.frame(type = c(rep("pvalue", nrow(res)), rep("pvalue adjusted", nrow(res))),
                        pvalue = c(res$pvalue, res$padj))

plot.pvalue.bind <- ggplot(pvalue.df, aes(x=pvalue, fill=type)) + 
  geom_histogram(aes(y=..density..), binwidth=0.05, color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic(base_size = 11) +
  labs(fill="")

ggsave("Figures/histogram_pvalue_bind.png", plot = plot.pvalue.bind)
rm(plot.pvalue.bind)

plot.pvalue <- ggplot(data.frame(res), aes(x=pvalue)) +
  geom_histogram(aes(y=..density..), binwidth=0.07, fill="#69b3a2", color="#e9ecef", alpha=0.6, position = 'identity') +
  theme_classic(base_size = 11)

plot.padj <- ggplot(data.frame(res), aes(x=padj)) +
  geom_histogram(aes(y=..density..), binwidth=0.07, fill="#404080", color="#e9ecef", alpha=0.6, position = 'identity') +
  theme_classic(base_size = 11)

require(gridExtra)
plot.pvalue.unbind <- grid.arrange(plot.pvalue, plot.padj, ncol=2)

ggsave("Figures/histogram_pvalue_unbind.png", plot = plot.pvalue.unbind)
rm(plot.pvalue, plot.padj, plot.pvalue.unbind)

### MA plot ###

par(mfrow=c(1,1))

png("Figures/MA_plot.png", width = 800, height = 800)
plotMA(res, ylim=c(-7,7))
dev.off()

### Count plot ###

res <- res[order(res$padj),] # Order by adjusted p-value

png("Figures/count_plot.png", width = 800, height = 800)
par(mfrow=c(2,3))
for (i in 1:6){
  plotCounts(dds, gene=rownames(res)[i], intgroup="mutation")
}
dev.off()

### Volcano plot ###

res <- results(dds, tidy=TRUE)
res.condition1 <- subset(res, padj<.1)
res.condition2 <- subset(res, (padj<.005 & -log10(pvalue) > 6))

plot.volcano <- ggplot(data = res.condition2, aes(x = log2FoldChange, y = -log10(pvalue), label = row)) +
  geom_point(data = res, aes(x = log2FoldChange, y = -log10(pvalue)), color = "black") +
  geom_point(data = res.condition1, aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue") +
  geom_point(color = "red") +
  geom_text_repel(size=3) +
  xlim(-7.3, 7.3) +
  theme_classic(base_size = 11)

ggsave("Figures/volcano_plot.png", plot = plot.volcano)
rm(plot.volcano)

