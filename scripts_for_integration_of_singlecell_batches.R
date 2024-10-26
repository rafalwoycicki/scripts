library(dplyr)
library(Seurat)
library(patchwork)

#Loading the data

SoyBjap1_v4.data <- Read10X(data.dir = "./SoyBjap1_v4_filtered_feature_bc_matrix/")
SoyBjap2_v4.data <- Read10X(data.dir = "./SoyBjap2_v4_filtered_feature_bc_matrix/")
SoyBjap3_v4.data <- Read10X(data.dir = "./SoyBjap3_v4_filtered_feature_bc_matrix/")
SoyMock1_v4.data <- Read10X(data.dir = "./SoyMock1_v4_filtered_feature_bc_matrix/")
SoyMock2_v4.data <- Read10X(data.dir = "./SoyMock2_v4_filtered_feature_bc_matrix/")
SoyMock3_v4.data <- Read10X(data.dir = "./SoyMock3_v4_filtered_feature_bc_matrix/")

#####Ensembl test ONLY
#SoyBjap1_v4.data <- Read10X(data.dir = "./SoyBjap1_Ensembl_filtered_feature_bc_matrix/")
#SoyBjap2_v4.data <- Read10X(data.dir = "./SoyBjap2_Ensembl_filtered_feature_bc_matrix/")
#SoyBjap3_v4.data <- Read10X(data.dir = "./SoyBjap3_Ensembl_filtered_feature_bc_matrix/")
#SoyMock1_v4.data <- Read10X(data.dir = "./SoyMock1_Ensembl_filtered_feature_bc_matrix/")
#SoyMock2_v4.data <- Read10X(data.dir = "./SoyMock2_Ensembl_filtered_feature_bc_matrix/")
#SoyMock3_v4.data <- Read10X(data.dir = "./SoyMock3_Ensembl_filtered_feature_bc_matrix/")
############

# Creating Seurat Objects
SoyBjap1_v4 <- CreateSeuratObject(counts = SoyBjap1_v4.data, project = "SoyBjap1_v4", min.cells = 3, min.features = 200)
SoyBjap2_v4 <- CreateSeuratObject(counts = SoyBjap2_v4.data, project = "SoyBjap2_v4", min.cells = 3, min.features = 200)
SoyBjap3_v4 <- CreateSeuratObject(counts = SoyBjap3_v4.data, project = "SoyBjap3_v4", min.cells = 3, min.features = 200)
SoyMock1_v4 <- CreateSeuratObject(counts = SoyMock1_v4.data, project = "SoyMock1_v4", min.cells = 3, min.features = 200)
SoyMock2_v4 <- CreateSeuratObject(counts = SoyMock2_v4.data, project = "SoyMock2_v4", min.cells = 3, min.features = 200)
SoyMock3_v4 <- CreateSeuratObject(counts = SoyMock3_v4.data, project = "SoyMock3_v4", min.cells = 3, min.features = 200)

# Non nuclear features counting

SoyBjap1_v4[["percent.U"]] <- PercentageFeatureSet(SoyBjap1_v4, pattern = "^Glyma.U")
SoyBjap2_v4[["percent.U"]] <- PercentageFeatureSet(SoyBjap2_v4, pattern = "^Glyma.U")
SoyBjap3_v4[["percent.U"]] <- PercentageFeatureSet(SoyBjap3_v4, pattern = "^Glyma.U")
SoyMock1_v4[["percent.U"]] <- PercentageFeatureSet(SoyMock1_v4, pattern = "^Glyma.U")
SoyMock2_v4[["percent.U"]] <- PercentageFeatureSet(SoyMock2_v4, pattern = "^Glyma.U")
SoyMock3_v4[["percent.U"]] <- PercentageFeatureSet(SoyMock3_v4, pattern = "^Glyma.U")

head(SoyBjap1_v4@meta.data, 5)
head(SoyBjap2_v4@meta.data, 5)
head(SoyBjap3_v4@meta.data, 5)
head(SoyMock1_v4@meta.data, 5)
head(SoyMock2_v4@meta.data, 5)
head(SoyMock3_v4@meta.data, 5)

# exploratory plots

VlnPlot(SoyBjap1_v4, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)
plot1 <- FeatureScatter(SoyBjap1_v4, feature1 = "nCount_RNA", feature2 = "percent.U")
plot2 <- FeatureScatter(SoyBjap1_v4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("filtering_SoyBjap1_v4.pdf")
plot1 + plot2
dev.off()
plot1 + plot2

VlnPlot(SoyBjap2_v4, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)
plot1 <- FeatureScatter(SoyBjap2_v4, feature1 = "nCount_RNA", feature2 = "percent.U")
plot2 <- FeatureScatter(SoyBjap2_v4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("filtering_SoyBjap2_v4.pdf")
plot1 + plot2
dev.off()
plot1 + plot2

VlnPlot(SoyBjap3_v4, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)
plot1 <- FeatureScatter(SoyBjap3_v4, feature1 = "nCount_RNA", feature2 = "percent.U")
plot2 <- FeatureScatter(SoyBjap3_v4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("filtering_SoyBjap3_v4.pdf")
plot1 + plot2
dev.off()
plot1 + plot2

VlnPlot(SoyMock1_v4, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)
plot1 <- FeatureScatter(SoyMock1_v4, feature1 = "nCount_RNA", feature2 = "percent.U")
plot2 <- FeatureScatter(SoyMock1_v4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("filtering_SoyMock1_v4.pdf")
plot1 + plot2
dev.off()
plot1 + plot2

VlnPlot(SoyMock2_v4, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)
plot1 <- FeatureScatter(SoyMock2_v4, feature1 = "nCount_RNA", feature2 = "percent.U")
plot2 <- FeatureScatter(SoyMock2_v4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("filtering_SoyMock2_v4.pdf")
plot1 + plot2
dev.off()
plot1 + plot2

VlnPlot(SoyMock3_v4, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)
plot1 <- FeatureScatter(SoyMock3_v4, feature1 = "nCount_RNA", feature2 = "percent.U")
plot2 <- FeatureScatter(SoyMock3_v4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("filtering_SoyMock3_v4.pdf")
plot1 + plot2
dev.off()
plot1 + plot2

# Filtering datasets per min number of genes and with low number of non-nuclear genes

SoyBjap1_v4 <- subset(SoyBjap1_v4, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.U < 0.5)
SoyBjap2_v4 <- subset(SoyBjap2_v4, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.U < 0.5)
SoyBjap3_v4 <- subset(SoyBjap3_v4, subset = nFeature_RNA > 200 & nFeature_RNA < 8750 & percent.U < 0.5)
SoyMock1_v4 <- subset(SoyMock1_v4, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.U < 0.5)
SoyMock2_v4 <- subset(SoyMock2_v4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.U < 0.5)
SoyMock3_v4 <- subset(SoyMock3_v4, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.U < 0.5)

# Normalization

SoyBjap1_v4 <- NormalizeData(SoyBjap1_v4, verbose=FALSE)
SoyBjap2_v4 <- NormalizeData(SoyBjap2_v4, verbose=FALSE)
SoyBjap3_v4 <- NormalizeData(SoyBjap3_v4, verbose=FALSE)
SoyMock1_v4 <- NormalizeData(SoyMock1_v4, verbose=FALSE)
SoyMock2_v4 <- NormalizeData(SoyMock2_v4, verbose=FALSE)
SoyMock3_v4 <- NormalizeData(SoyMock3_v4, verbose=FALSE)

# Variable features selection

SoyBjap1_v4 <- FindVariableFeatures(SoyBjap1_v4, selection.method="vst", nfeatures=2000, verbose=FALSE)
SoyBjap2_v4 <- FindVariableFeatures(SoyBjap2_v4, selection.method="vst", nfeatures=2000, verbose=FALSE)
SoyBjap3_v4 <- FindVariableFeatures(SoyBjap3_v4, selection.method="vst", nfeatures=2000, verbose=FALSE)
SoyMock1_v4 <- FindVariableFeatures(SoyMock1_v4, selection.method="vst", nfeatures=2000, verbose=FALSE)
SoyMock2_v4 <- FindVariableFeatures(SoyMock2_v4, selection.method="vst", nfeatures=2000, verbose=FALSE)
SoyMock3_v4 <- FindVariableFeatures(SoyMock3_v4, selection.method="vst", nfeatures=2000, verbose=FALSE)

# Merging of batches

reference.list <- c(SoyBjap1_v4, SoyBjap2_v4, SoyBjap3_v4, SoyMock1_v4, SoyMock2_v4, SoyMock3_v4)

anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)

integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

library(ggplot2)
library(cowplot)
library(patchwork)

# Work on integrated data: scaling, PCA, Clustering, UMAP

DefaultAssay(integrated) <- "integrated"

all.genes <- rownames(integrated)

integrated <- ScaleData(integrated, features = all.genes, verbose = FALSE)

integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)

ElbowPlot(integrated, ndims=50) # it could have been 100

integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:50)

integrated <- FindClusters(integrated, resolution = 0.5) #23 clusters

#integrated <- RunTSNE(integrated, dims = 1:50, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)

integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:50)

# cluster visualisations

integrated$sample <- plyr::mapvalues(x = integrated$orig.ident, from = c("SoyBjap1_v4", "SoyBjap2_v4", "SoyBjap3_v4", "SoyMock1_v4", "SoyMock2_v4", "SoyMock3_v4"), to = c("SoyBjapv4", "SoyBjapv4", "SoyBjapv4", "SoyMockv4", "SoyMockv4", "SoyMockv4"))

saveRDS(integrated, file = "./integrated_v4_r05_c23.rds")

pdf("integratedUMAP_v4_spbSAMPLE.pdf")
DimPlot(integrated, reduction = "umap", split.by = "sample")
dev.off()
pdf("integratedUMAP_v4_grbSAMPLE.pdf")
DimPlot(integrated, reduction = "umap", group.by = "sample")
dev.off()
pdf("integratedUMAP_v4_grbSAMPLE_spbIDENTS.pdf")
DimPlot(integrated, reduction = "umap", group.by = "sample", split.by = "orig.ident")
dev.off()
pdf("integratedUMAP_v4_grbIDENT.pdf")
DimPlot(integrated, reduction = "umap", group.by = "orig.ident")
dev.off()
pdf("integratedUMAP_v4_spbIDENT.pdf")
DimPlot(integrated, reduction = "umap", split.by = "orig.ident")
dev.off()
pdf("integratedUMAP_v4_grbIDENTS_spbSAMPLE.pdf")
DimPlot(integrated, reduction = "umap", group.by = "orig.ident", split.by = "sample")
dev.off()

#Identify conserved markers:






