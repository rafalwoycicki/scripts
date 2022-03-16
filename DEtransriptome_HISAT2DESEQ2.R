### R script to analyze DE of genes using DESeq2 on the basis of transcriptome mappings done by the HiSAT2. 

library(DESeq2)

# loading the data
countData <- as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))
colData <- read.csv("samplesHDt1.txt", sep="\t", row.names=1)
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

#creating DESeq DataSet
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ type)

# filtering out genes with lowest number of reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
nrow(dds)

# VST transformations and plots for data visualisation
library(vsn)
vsd <- vst(dds, blind=FALSE)
vsd.blind <- vst(dds, blind=TRUE)

library("dplyr")
library("ggplot2")
sampleDists.vsd <- dist(t(assay(vsd)))
sampleDists.vsd.blind <- dist(t(assay(vsd.blind)))

library("pheatmap")
library("RColorBrewer")

pdf('PCAplot_vsd.noblind_TYPE.pdf')
pcaData <- plotPCA(vsd, intgroup = c("type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color=type)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with VST data")
dev.off()

pdf('PCAplot_vsd.Blind_TYPE.pdf')
pcaData <- plotPCA(vsd.blind, intgroup = c("type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color=type)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with VST data")
dev.off()

pdf('matrix_vsd.genes.NoBlind.pdf')
sampleDistMatrix.vsd <- as.matrix( sampleDists.vsd )
rownames(sampleDistMatrix.vsd) <- NULL
colnames(sampleDistMatrix.vsd) <- vsd$type
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.vsd, clustering_distance_rows = sampleDists.vsd, clustering_distance_cols = sampleDists.vsd, col = colors)
dev.off()

pdf('matrix_vsd.genes.Blind.pdf')
sampleDistMatrix.vsd <- as.matrix( sampleDists.vsd.blind )
rownames(sampleDistMatrix.vsd) <- NULL
colnames(sampleDistMatrix.vsd) <- vsd.blind$type
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.vsd, clustering_distance_rows = sampleDists.vsd.blind, clustering_distance_cols = sampleDists.vsd.blind, col = colors)
dev.off()

library("PoiClaClu")

poisd <- PoissonDistance(t(counts(dds)))

pdf('poissondistane_Blind.pdf')
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- NULL
colnames(samplePoisDistMatrix) <- vsd.blind$type
pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col = colors)
dev.off()

library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$type <- dds$type

pdf('GLMPCAplot.genes.pdf')
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = type)) + geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
dev.off()

pdf('MDSwithVSTplot_vsd.blind.genes.pdf')
sampleDistMatrix.vsd <- as.matrix( sampleDists.vsd.blind )
mds <- as.data.frame(colData(vsd.blind))  %>% cbind(cmdscale(sampleDistMatrix.vsd))
ggplot(mds, aes(x = `1`, y = `2`, color = type)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
dev.off()

pdf('MDSwithPoissonDistplot.genes.pdf')
mdsPois <- as.data.frame(colData(dds)) %>% cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = type)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
dev.off()

# rlog transformation and data visualisations

rld <- rlog(dds, blind=FALSE, fitType='local')
sampleDists.rld <- dist(t(assay(rld)))

rld.blind <- rlog(dds, blind=TRUE, fitType='local')
sampleDists.rld.blind <- dist(t(assay(rld)))

pdf('matrix_rld.NoBlind.genes.pdf')
sampleDistMatrix.rld <- as.matrix( sampleDists.rld )
rownames(sampleDistMatrix.rld) <- NULL
colnames(sampleDistMatrix.rld) <- rld$type
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.rld, clustering_distance_rows = sampleDists.rld, clustering_distance_cols = sampleDists.rld, col = colors)
dev.off()

pdf('matrix_rld.Blind.genes.pdf')
sampleDistMatrix.rld.blind <- as.matrix( sampleDists.rld.blind )
rownames(sampleDistMatrix.rld.blind) <- NULL
colnames(sampleDistMatrix.rld.blind) <- rld$type
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.rld.blind, clustering_distance_rows = sampleDists.rld.blind, clustering_distance_cols = sampleDists.rld.blind, col = colors)
dev.off()

pdf('PCAplot_rld.genes.pdf')
pcaData <- plotPCA(rld, intgroup = c("type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with rld.genes data")
dev.off()

pdf('PCAplot_rld.blind.genes.pdf')
pcaData <- plotPCA(rld.blind, intgroup = c("type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + ggtitle("PCA with rld.blind.genes data")
dev.off()

pdf('MDSwithrld.genesplot.genes.pdf')
sampleDistMatrix.rld <- as.matrix( sampleDists.rld )
mds <- as.data.frame(colData(rld))  %>% cbind(cmdscale(sampleDistMatrix.rld))
ggplot(mds, aes(x = `1`, y = `2`, color = type)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with rld.genes data")
dev.off()

pdf('MDSwithrld.blind.genesplot.genes.pdf')
sampleDistMatrix.rld.blind <- as.matrix( sampleDists.rld.blind )
mds <- as.data.frame(colData(rld.blind))  %>% cbind(cmdscale(sampleDistMatrix.rld.blind))
ggplot(mds, aes(x = `1`, y = `2`, color = type)) + geom_point(size = 3) + coord_fixed() + ggtitle("MDS with rld.genes data")
dev.off()

# DE analysis
#dds.deseq <- DESeq(dds)
dds.deseq <- DESeq(dds, fitType='local')

###############################################################

#137 vs 140

res.transcripts.type.137_140_all <- results(dds.deseq, contrast=c("type","s137","s140"))
summary(res.transcripts.type.137_140_all)
write.table(res.transcripts.type.137_140_all,file="results.transcripts.type.137_140_all.table",quote=FALSE,sep="\t")

#######################################################################

#137 vs 181

res.transcripts.type.137_181_all <- results(dds.deseq, contrast=c("type","s137","s181"))
summary(res.transcripts.type.137_181_all)
write.table(res.transcripts.type.137_181_all,file="results.transcripts.type.137_181_all.table",quote=FALSE,sep="\t")

########################################################################

#140 vs 181

res.transcripts.type.140_181_all <- results(dds.deseq, contrast=c("type","s140","s181"))
summary(res.transcripts.type.140_181_all)
write.table(res.transcripts.type.140_181_all,file="results.transcripts.type.140_181_all.table",quote=FALSE,sep="\t")

##############################################################################




