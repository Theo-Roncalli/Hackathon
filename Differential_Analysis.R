#Desq2 accepte une matrice comptage qui contient en lignes les génes et en colonnes les échantillons (dans notre cas 8).

#Premiérement, créer un objet DEseq Dataset 
#Importation du fichier

countsTable <- read.delim("~/Desktop/Master/Master2/hackathon/counts.txt",header = T)

rownames(countsTable) <- countsTable[,1]
countsTable <- countsTable[,c(-1,-2,-3,-4,-5,-6)]

colnames(countsTable)<-c("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589")


dim(countsTable)
#nous avons 60676 génes pour 8 échantillons


#on prépare le vecteur des conditions pour les échantillons
conds<-factor(c("Mut","Mut","Mut","WT","WT","WT","WT","Mut"))

#création de colData (utile pour la suite)
colData=data.frame(condition=conds)


#on importe la package.
library(DESeq2)
require(MASS)
require(calibrate)



#Création de l'objet DESeq 
dds<- DESeqDataSetFromMatrix(countsTable,colData,design= ~condition)

#cette présicion est utile pour la modélisation des données RNA-seq. 
#On précise que les reads de chaque géne dépendent de la variable condition (i.e si muté ou non )
#c'est un choix de modélisation , on aurait pu faire autre chose. 


dds<-DESeq(dds)
#Tout est inclu dans la pipeline DESeq , il fait tout ça : 

#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

#Toutes ces étapes peuvent être éffectués individuellement 



#On plot la dispersion des Reads.

png("~/dispersions_reads.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()


#On peut choisir de faire une transformation log2 pour normaliser nos données et minimiser les différences de Reads entre échantillons.

rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(conds))])


sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots) 

png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20) 
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[conds], RowSideColors=mycols[conds],
          margin=c(10, 10), main="Sample Distance Matrix")

dev.off()


####Resultats#####
#on tire les résultats de la pipeline DESeq
res <- results(dds)

table(res$padj<0.05)

summary(res)

#ordonner selon les p-values (pour le volcano plot)
res <- res[order(res$padj),]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

#vérification des pvalues

#hist(res$pvalue, breaks=50, col="grey")


######Plot CPA#####



png("PCArld.png", w=1000, h=1000, pointsize=20)

DESeq2::plotPCA(rld, intgroup="condition") 
dev.off()


vsdata <- vst(dds, blind=FALSE) #on applique une stabilisation de la varaiance

png("PCA.png", w=1000, h=1000, pointsize=20)

plotPCA(vsdata, intgroup="condition")

dev.off()



##########Volcano plot#####



volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topleft", labelsig=FALSE, textcx=1, Stock =TRUE ,...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) { #FALSE donc pas de nom de génes dans le plot
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  if (Stock){
    volcano_Gene_threshold <- subset(res, padj<thresh)[,1] # je ne garde que la premiére colonne , donc les noms de génes.
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}


png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-4.3, 2))

dev.off()



#####MA plot ####


maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, Stock = TRUE ,...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2)) #textxy permet d'ajuster les libéllées dans une figure
  }
  if (Stock){
    MA_Gene_threshold <- subset(res, padj<thresh)[,1] # je ne garde que la premiére colonne , donc les noms de génes.
  }
}


png("diffexpr-MAplot.png", 1200, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()






