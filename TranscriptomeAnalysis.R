# preparation of data
samples <- read.table(file.path("/RNA_Trial/salmon/pangenome", "samples.txt"), header=TRUE)
files.NoDup <- file.path("/RNA_Trial/salmon/pangenome", "quants_NoDup", samples$file, "quant.sf")
names(files.NoDup) <- samples$name
all(file.exists(files.NoDup))

# adding gene annotation
library(readr)
tx2gene <- read_csv(file.path("/RNA_Trial/ref", "Bd30_txdb.csv"))
head(tx2gene)
library(tximport)
txi.genes.NoDup <- tximport(files.NoDup, type = "salmon", tx2gene = tx2gene)
head(txi.genes.NoDup$counts)
txi.trans.NoDup <- tximport(files.NoDup, type = "salmon", txOut = TRUE)
head(txi.trans.NoDup$counts)
txi.sum.NoDup <- summarizeToGene(txi.trans.NoDup, tx2gene)
all.equal(txi.genes.NoDup$counts, txi.sum.NoDup$counts)
head(txi.sum.NoDup$counts)
samples_DESeq <- samples
rownames(samples_DESeq) <- samples$name
samples_DESeq

library(DESeq2)

# creating DESeq2 DataSet
#dds.trans.noDup <- DESeqDataSetFromTximport(txi.trans.NoDup, colData = samples_DESeq, design = ~ stand)
dds.trans.noDup <- DESeqDataSetFromTximport(txi.trans.NoDup, colData = samples_DESeq, design = ~ type)

# filtering out genes with lowest number of reads
keep <- rowSums(counts(dds.trans.noDup)) >= 10
dds.trans.noDup.10 <- dds.trans.noDup[keep,]
nrow(dds.trans.noDup.10)

# transformation and data visualisations
library(vsn)
vsd.trans <- vst(dds.trans.noDup.10, blind=FALSE)
rld.trans <- rlog(dds.trans.noDup.10, blind=FALSE)

library("dplyr")
library("ggplot2")


dds.trans.noDup.10.eSF <- estimateSizeFactors(dds.trans.noDup.10)

###########################################################################

df <- bind_rows(as_data_frame(log2(counts(dds.trans.noDup.10.eSF,  normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd.trans)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld.trans)[, 1:2]) %>% mutate(transformation = "rlog"))
  
colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
ggsave("comparison.pdf", plot=last_plot(),device="pdf")

#####################################################

sampleDists.vsd.trans <- dist(t(assay(vsd.trans)))
sampleDists.rld.trans <- dist(t(assay(rld.trans)))

##########################################

library("pheatmap")
library("RColorBrewer")


######################################################
pdf('matrix_vsd.trans.pdf')
sampleDistMatrix.vsd.trans <- as.matrix( sampleDists.vsd.trans )
rownames(sampleDistMatrix.vsd.trans) <- NULL
colnames(sampleDistMatrix.vsd.trans) <- vsd.trans$name
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.vsd.trans,
         clustering_distance_rows = sampleDists.vsd.trans,
         clustering_distance_cols = sampleDists.vsd.trans,
         col = colors)
dev.off()

pdf('matrix_rld.trans.pdf')
sampleDistMatrix.rld.trans <- as.matrix( sampleDists.rld.trans )
rownames(sampleDistMatrix.rld.trans) <- NULL
colnames(sampleDistMatrix.rld.trans) <- vsd.trans$name
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.rld.trans,
         clustering_distance_rows = sampleDists.rld.trans,
         clustering_distance_cols = sampleDists.rld.trans,
         col = colors)
dev.off()

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds.trans.noDup.10)))

pdf('poissondistance.pdf')
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds.trans.noDup.10$stand, dds.trans.noDup.10$name, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col = colors)
dev.off()

pdf('PCAplot_vsd.trans.pdf')
plotPCA(vsd.trans, intgroup = c("region", "stand"))
dev.off()

pdf('PCAplot_rld.trans.pdf')
plotPCA(rld.trans, intgroup = c("region", "stand"))
dev.off()

pdf('PCAplot_vsd.trans2.pdf')
pcaData <- plotPCA(vsd.trans, intgroup = c( "region", "stand"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = stand, shape = region)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with VST data")
dev.off()

pdf('PCAplot_rld.trans2.pdf')
pcaData <- plotPCA(rld.trans, intgroup = c( "region", "stand"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = stand, shape = region)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with rld.trans data")
dev.off()

library("glmpca")

gpca <- glmpca(counts(dds.trans.noDup.10), L=2)
gpca.dat <- gpca$factors
gpca.dat$region <- dds.trans.noDup.10$region
gpca.dat$stand <- dds.trans.noDup.10$stand
pdf('GLMPCAplot_vsd.trans2.pdf')
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = stand, shape = region)) + geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
dev.off()

mds <- as.data.frame(colData(vsd.trans))  %>%
         cbind(cmdscale(sampleDistMatrix.vsd.trans))
pdf('MDSwithVSTplot_vsd.trans2.pdf')
ggplot(mds, aes(x = `1`, y = `2`, color = stand, shape = region)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
dev.off()

mds <- as.data.frame(colData(rld.trans))  %>%
         cbind(cmdscale(sampleDistMatrix.rld.trans))
pdf('MDSwithrld.transplot_vsd.trans2.pdf')
ggplot(mds, aes(x = `1`, y = `2`, color = stand, shape = region)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with rld.trans data")
dev.off()

mdsPois <- as.data.frame(colData(dds.trans.noDup.10)) %>%
   cbind(cmdscale(samplePoisDistMatrix))
pdf('MDSwithPoissonDistplot_vsd.trans2.pdf')
ggplot(mdsPois, aes(x = `1`, y = `2`, color = stand, shape = region)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
dev.off()

# DE analysis

dds.trans.noDup.10.deseq <- DESeq(dds.trans.noDup.10)

#region

dds.trans.noDup.region <- DESeqDataSetFromTximport(txi.trans.NoDup, colData = samples_DESeq, design = ~ region)

keep <- rowSums(counts(dds.trans.noDup.region)) >= 10

dds.trans.noDup.region.10 <- dds.trans.noDup.region[keep,]

dds.trans.noDup.region.10.deseq <- DESeq(dds.trans.noDup.region.10)


#type

dds.trans.noDup.type <- DESeqDataSetFromTximport(txi.trans.NoDup, colData = samples_DESeq, design = ~ type)

keep <- rowSums(counts(dds.trans.noDup.type)) >= 10

dds.trans.noDup.type.10 <- dds.trans.noDup.type[keep,]

dds.trans.noDup.type.10.deseq <- DESeq(dds.trans.noDup.type.10)


###############################################################################################################

###STAND

res.trans.stand.s26_s34 <- results(dds.trans.noDup.10.deseq, contrast=c("stand","s26","s34"))

summary(res.trans.stand.s26_s34)

res.trans.stand.s26_s34 <- results(dds.trans.noDup.10.deseq, contrast=c("stand","s26","s34"), lfcThreshold=1, alpha=0.01)

write.table(res.trans.stand.s26_s34,file="res.trans.stand.s26_s34.table",quote=FALSE,sep="\t")

summary(res.trans.stand.s26_s34)



topGene <- rownames(res.trans.stand.s26_s34)[which.min(res.trans.stand.s26_s34$padj)]

#
pdf('topvargenesVSDTRANSSTands26s34PLOT.pdf')
plotCounts(dds.trans.noDup.10.deseq, gene = topGene, intgroup=c("stand"))
dev.off()
#

library("ggbeeswarm")

geneCounts <- plotCounts(dds.trans.noDup.10.deseq, gene = topGene, intgroup = c("stand","name"), returnData=TRUE)

pdf('topvargenesVSDTRANSSTands26s34.pdf')
ggplot(geneCounts, aes(x = stand, y = count, color = name)) + scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()


library("genefilter")

topVarGenes <- head(order(rowVars(assay(vsd.trans)), decreasing = TRUE), 20)

mat  <- assay(vsd.trans)[ topVarGenes, ]

mat  <- mat - rowMeans(mat)

anno <- as.data.frame(colData(vsd.trans)[, c("name","stand")])

pdf('topvargenesVSDTRANSStand.pdf')
pheatmap(mat, annotation_col = anno)
dev.off()



###REGION
#R4vsR3

res.trans.region.r4_r3 <- results(dds.trans.noDup.region.10.deseq, contrast=c("region","r4","r3"))

summary(res.trans.region.r4_r3)

res.trans.region.r4_r3 <- results(dds.trans.noDup.region.10.deseq, contrast=c("region","r4","r3"), lfcThreshold=1, alpha=0.01)

summary(res.trans.region.r4_r3)

write.table(res.trans.region.r4_r3,file="res.trans.region.r4_r3.table",quote=FALSE,sep="\t")



topGene <- rownames(res.trans.region.r4_r3)[which.min(res.trans.region.r4_r3$padj)]


#
pdf('topvargenesVSDTRANSRegionR3R4PLOT.pdf')
plotCounts(dds.trans.noDup.region.10.deseq, gene = topGene, intgroup=c("region"))
dev.off()
#

library("ggbeeswarm")

geneCounts <- plotCounts(dds.trans.noDup.region.10.deseq, gene = topGene, intgroup = c("region","name"), returnData=TRUE)

pdf('topvargenesVSDTRANSRegionR3R4.pdf')
ggplot(geneCounts, aes(x = region, y = count, color = name)) + scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()




#R1c vs R1a

res.trans.region.r1c_r1a <- results(dds.trans.noDup.region.10.deseq, contrast=c("region","r1c","r1a"))

summary(res.trans.region.r1c_r1a)

res.trans.region.r1c_r1a <- results(dds.trans.noDup.region.10.deseq, contrast=c("region","r1c","r1a"), lfcThreshold=1, alpha=0.01)

summary(res.trans.region.r1c_r1a)

write.table(res.trans.region.r1c_r1a,file="res.trans.region.r1c_r1a.table",quote=FALSE,sep="\t")



topGene <- rownames(res.trans.region.r1c_r1a)[which.min(res.trans.region.r1c_r1a$padj)]

#
pdf('topvargenesVSDTRANSRegion1c1aPLOT.pdf')
plotCounts(dds.trans.noDup.region.10.deseq, gene = topGene, intgroup=c("region"))
dev.off()
#

library("ggbeeswarm")

geneCounts <- plotCounts(dds.trans.noDup.region.10.deseq, gene = topGene, intgroup = c("region","name"), returnData=TRUE)

pdf('topvargenesVSDTRANSRegion1c1a.pdf')
ggplot(geneCounts, aes(x = region, y = count, color = name)) + scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()



vsd.trans.region <- vst(dds.trans.noDup.region.10, blind=FALSE)

rld.trans.region <- rlog(dds.trans.noDup.region.10, blind=FALSE)

library("genefilter")

topVarGenes <- head(order(rowVars(assay(vsd.trans.region)), decreasing = TRUE), 20)

mat  <- assay(vsd.trans.region)[ topVarGenes, ]

mat  <- mat - rowMeans(mat)

anno <- as.data.frame(colData(vsd.trans.region)[, c("name","region")])

pdf('topvargenesVSDTRANSRegion.pdf')
pheatmap(mat, annotation_col = anno)
dev.off()




###TYPE
#central vs coastal

res.trans.type.tcen_tcoa <- results(dds.trans.noDup.type.10.deseq, contrast=c("type","central","coastal"))

summary(res.trans.type.tcen_tcoa)

res.trans.type.tcen_tcoa <- results(dds.trans.noDup.type.10.deseq, contrast=c("type","central","coastal"), lfcThreshold=1, alpha=0.01)

summary(res.trans.type.tcen_tcoa)

write.table(res.trans.type.tcen_tcoa,file="res.trans.type.tcen_tcoa.table",quote=FALSE,sep="\t")



topGene <- rownames(res.trans.type.tcen_tcoa)[which.min(res.trans.type.tcen_tcoa$padj)]

#
pdf('topgeneDDSTRANSTypePLOT.pdf')
plotCounts(dds.trans.noDup.type.10.deseq, gene = topGene, intgroup=c("type"))
dev.off()
#

library("ggbeeswarm")

geneCounts <- plotCounts(dds.trans.noDup.type.10.deseq, gene = topGene, intgroup = c("type","name"), returnData=TRUE)

pdf('topgeneDDSTRANSType.pdf')
ggplot(geneCounts, aes(x = type, y = count, color = name)) + scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

vsd.trans.type <- vst(dds.trans.noDup.type.10, blind=FALSE)

rld.trans.type <- rlog(dds.trans.noDup.type.10, blind=FALSE)

library("genefilter")

topVarGenes <- head(order(rowVars(assay(vsd.trans.type)), decreasing = TRUE), 20)

mat  <- assay(vsd.trans.type)[ topVarGenes, ]

mat  <- mat - rowMeans(mat)

anno <- as.data.frame(colData(vsd.trans.type)[, c("name","type")])

pdf('topvargenesVSDTRANSType.pdf')
pheatmap(mat, annotation_col = anno)
dev.off()



#########################################################################################

DATA QUALITY ASSESMENT

library("pheatmap")

dds.trans.noDup.10.eSF <- estimateSizeFactors(dds.trans.noDup.10)

select <- order(rowMeans(counts(dds.trans.noDup.10.eSF,normalized=TRUE)),
                decreasing=TRUE)[1:50]

df <- as.data.frame(colData(dds.trans.noDup.10.eSF)[,c("region","stand")])

pheatmap(assay(dds.trans.noDup.10.eSF)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd.trans)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld.trans)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds.trans.noDup)[["cooks"]]), range=0, las=2)


plotDispEsts(dds.trans.noDup.10.eSF)





