rm(list=objects())
graphics.off()

setwd("/home/ubuntu/Documents/AMI2B/Hackathon/Projet/R")

counts <- read.table("counts.txt",header = T, row.names = 1)
colnames(counts) <- sub("\\.bam", "", colnames(counts))
countData <- counts[rowSums(counts[,6:13]) > 0,][,6:13]

metadata <- read.table("SraRunTable.txt",header = T, row.names = 1, sep=",")

id <- rownames(metadata)

tumor <- metadata$tumor_class
tumor <- replace(tumor, tumor=="", "Class 2")

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

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~mutation, tidy = FALSE)

dds <- DESeq(dds)

# estimateSizeFactors => computes the relative library depth of each sample 
# estimateDispersions => estimates the dispersion of counts for each gene
# nbinomWaldTest => computes the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs

res <- results(dds, tidy=FALSE)
summary(res)

res <- res[order(res$padj),] # Order by p-value
head(res)

### Count plot ###

par(mfrow=c(2,3))

for (i in 1:6){
  plotCounts(dds, gene=rownames(res)[i], intgroup="mutation")
}

# Next steps in exploring these data...BLAST to database to find associated gene function

### Volcano plot ###

par(mfrow=c(1,1))

plot(res$log2FoldChange, -log10(res$pvalue), pch=20, main="Volcano plot", xlim=c(-7.3,7.3), xlab="log2FoldChange", ylab="-log(pvalue)")
res.condition1 <- subset(res, padj<.01)
points(res.condition1$log2FoldChange, -log10(res.condition1$pvalue), pch=20, col="blue")
res.condition2 <- subset(res, padj<.01 & abs(log2FoldChange)>3)
points(res.condition2$log2FoldChange, -log10(res.condition2$pvalue), pch=20, col="red")

text(res.condition1$log2FoldChange, -log10(res.condition1$pvalue), labels=rownames(res.condition1), cex=0.6, font=2)

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
