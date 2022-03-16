samples <- read.table(file.path("", "samplesMajaSalmon.txt"), header=TRUE)

files.NoDup <- file.path("", "quants_TrinityFR_lISF", samples$file, "quant.sf")

names(files.NoDup) <- samples$name

all(file.exists(files.NoDup))

library(readr)
library(tximport)

tx2gene <- read_csv(file.path("", "Trinity.fasta_txdb.csv"))

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

#######################################################################################

dds.genes.noDup <- DESeqDataSetFromTximport(txi.genes.NoDup, colData = samples, design = ~ type)

keep <- rowSums(counts(dds.genes.noDup)) >= 10

dds.genes.noDup.10 <- dds.genes.noDup[keep,]

nrow(dds.genes.noDup.10)

library(vsn)

vsd.genes <- vst(dds.genes.noDup.10, blind=FALSE)

rld.genes <- rlog(dds.genes.noDup.10, blind=FALSE)

library("dplyr")
library("ggplot2")

#################################################################################

dds.genes.noDup.10.eSF <- estimateSizeFactors(dds.genes.noDup.10)

df <- bind_rows(as_data_frame(log2(counts(dds.genes.noDup.10.eSF,  normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"), as_data_frame(assay(vsd.genes)[, 1:2]) %>% mutate(transformation = "vst"), as_data_frame(assay(rld.genes)[, 1:2]) %>% mutate(transformation = "rlog"))
  
colnames(df)[1:2] <- c("x", "y") 

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)
ggsave("comparison.pdf", plot=last_plot(),device="pdf")

#############################################################################

sampleDists.vsd.genes <- dist(t(assay(vsd.genes)))
sampleDists.rld.genes <- dist(t(assay(rld.genes)))

library("pheatmap")
library("RColorBrewer")

pdf('matrix_vsd.genes.pdf')
sampleDistMatrix.vsd.genes <- as.matrix( sampleDists.vsd.genes )
rownames(sampleDistMatrix.vsd.genes) <- paste( vsd.genes$type, vsd.genes$name, sep = " - " )
colnames(sampleDistMatrix.vsd.genes) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.vsd.genes, clustering_distance_rows = sampleDists.vsd.genes, clustering_distance_cols = sampleDists.vsd.genes, col = colors)
dev.off()

pdf('matrix_rld.genes.pdf')
sampleDistMatrix.rld.genes <- as.matrix( sampleDists.rld.genes )
rownames(sampleDistMatrix.rld.genes) <- paste( vsd.genes$type, vsd.genes$name, sep = " - " )
colnames(sampleDistMatrix.rld.genes) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.rld.genes, clustering_distance_rows = sampleDists.rld.genes, clustering_distance_cols = sampleDists.rld.genes, col = colors)
dev.off()

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds.genes.noDup.10)))

pdf('poissondistance.pdf')
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds.genes.noDup.10$type, dds.genes.noDup.10$name, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col = colors)
dev.off()

pdf('PCAplot_vsd.genes.pdf')
plotPCA(vsd.genes, intgroup = "type")
dev.off()

pdf('PCAplot_rld.genes.pdf')
plotPCA(rld.genes, intgroup = "type")
dev.off()

pdf('PCAplot_vsd.genes2.pdf')
pcaData <- plotPCA(vsd.genes, intgroup = "type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with VST data")
dev.off()

pdf('PCAplot_rld.genes2.pdf')
pcaData <- plotPCA(rld.genes, intgroup = "type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with rld.genes data")
dev.off()

library("glmpca")

gpca <- glmpca(counts(dds.genes.noDup.10), L=2)
gpca.dat <- gpca$factors
gpca.dat$type <- dds.genes.noDup.10$type
pdf('GLMPCAplot_vsd.genes2.pdf')
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = type)) + geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
dev.off()

mds <- as.data.frame(colData(vsd.genes)) %>% cbind(cmdscale(sampleDistMatrix.vsd.genes))
pdf('MDSwithVSTplot_vsd.genes2.pdf')
ggplot(mds, aes(x = `1`, y = `2`, color = type)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
dev.off()

mds <- as.data.frame(colData(rld.genes)) %>% cbind(cmdscale(sampleDistMatrix.rld.genes))
pdf('MDSwithrld.genesplot_vsd.genes2.pdf')
ggplot(mds, aes(x = `1`, y = `2`, color = type)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with rld.genes data")
dev.off()

mdsPois <- as.data.frame(colData(dds.genes.noDup.10)) %>% cbind(cmdscale(samplePoisDistMatrix))
pdf('MDSwithPoissonDistplot_vsd.genes2.pdf')
ggplot(mdsPois, aes(x = `1`, y = `2`, color = type)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
dev.off()

#################################################################

dds.genes.noDup.10.deseq <- DESeq(dds.genes.noDup.10)

####################################################################################

#137 vs 140

res.genes.type.137_140_all <- results(dds.genes.noDup.10.deseq, contrast=c("type","s137","s140"))
summary(res.genes.type.137_140_all)

res.genes.type.137_140 <- results(dds.genes.noDup.10.deseq, contrast=c("type","s137","s140"), lfcThreshold=1, alpha=0.1)
summary(res.genes.type.137_140)

write.table(res.genes.type.137_140_all,file="res.genes.type.137_140_all.table",quote=FALSE,sep="\t")
topGene <- rownames(res.genes.type.137_140_all)[which.min(res.genes.type.137_140$padj)]

#
pdf('topgeneDDSGeneTypePLOT137_140.pdf')
plotCounts(dds.genes.noDup.10.deseq, gene = topGene, intgroup=c("type"))
dev.off()
#

library("ggbeeswarm")
geneCounts <- plotCounts(dds.genes.noDup.10.deseq, gene = topGene, intgroup = c("type","name"), returnData=TRUE)

pdf('topgeneDDSGeneType137_140.pdf')
ggplot(geneCounts, aes(x = type, y = count, color = name)) + scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

vsd.genes.type <- vst(dds.genes.noDup.10, blind=FALSE)
rld.genes.type <- rlog(dds.genes.noDup.10, blind=FALSE)

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd.genes.type)), decreasing = TRUE), 20)
mat  <- assay(vsd.genes.type)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd.genes.type)[, c("name","type")])

pdf('topvargenesVSDGenesType137_140.pdf')
pheatmap(mat, annotation_col = anno)
dev.off()
#########################################################################################

#137 vs 181

res.genes.type.137_181_all <- results(dds.genes.noDup.10.deseq, contrast=c("type","s137","s181"))
summary(res.genes.type.137_140_all)

res.genes.type.137_181 <- results(dds.genes.noDup.10.deseq, contrast=c("type","s137","s181"), lfcThreshold=1, alpha=0.1)
summary(res.genes.type.137_181)

write.table(res.genes.type.137_181_all,file="res.genes.type.137_181_all.table",quote=FALSE,sep="\t")
topGene <- rownames(res.genes.type.137_181_all)[which.min(res.genes.type.137_181$padj)]

#
pdf('topgeneDDSGeneTypePLOT137_181.pdf')
plotCounts(dds.genes.noDup.10.deseq, gene = topGene, intgroup=c("type"))
dev.off()
#

library("ggbeeswarm")
geneCounts <- plotCounts(dds.genes.noDup.10.deseq, gene = topGene, intgroup = c("type","name"), returnData=TRUE)

pdf('topgeneDDSGeneType137_181.pdf')
ggplot(geneCounts, aes(x = type, y = count, color = name)) + scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

vsd.genes.type <- vst(dds.genes.noDup.10, blind=FALSE)
rld.genes.type <- rlog(dds.genes.noDup.10, blind=FALSE)

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd.genes.type)), decreasing = TRUE), 20)
mat  <- assay(vsd.genes.type)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd.genes.type)[, c("name","type")])

pdf('topvargenesVSDGenesType137_181.pdf')
pheatmap(mat, annotation_col = anno)
dev.off()
#########################################################################################


#140 vs 181


res.genes.type.140_181_all <- results(dds.genes.noDup.10.deseq, contrast=c("type","s140","s181"))
summary(res.genes.type.140_181_all)

res.genes.type.140_181 <- results(dds.genes.noDup.10.deseq, contrast=c("type","s140","s181"), lfcThreshold=1, alpha=0.1)
summary(res.genes.type.140_181)

write.table(res.genes.type.140_181_all,file="res.genes.type.140_181_all.table",quote=FALSE,sep="\t")
topGene <- rownames(res.genes.type.140_181_all)[which.min(res.genes.type.140_181$padj)]

#
pdf('topgeneDDSGeneTypePLOT140_181.pdf')
plotCounts(dds.genes.noDup.10.deseq, gene = topGene, intgroup=c("type"))
dev.off()
#

library("ggbeeswarm")
geneCounts <- plotCounts(dds.genes.noDup.10.deseq, gene = topGene, intgroup = c("type","name"), returnData=TRUE)

pdf('topgeneDDSGeneType140_181.pdf')
ggplot(geneCounts, aes(x = type, y = count, color = name)) + scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

vsd.genes.type <- vst(dds.genes.noDup.10, blind=FALSE)
rld.genes.type <- rlog(dds.genes.noDup.10, blind=FALSE)

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd.genes.type)), decreasing = TRUE), 20)
mat  <- assay(vsd.genes.type)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd.genes.type)[, c("name","type")])

pdf('topvargenesVSDGenesType140_181.pdf')
pheatmap(mat, annotation_col = anno)
dev.off()
#########################################################################################

