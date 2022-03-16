
library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)

dir.create("results-BM123barnyard_5_500_mod")

#cat gem_classification.csv | tail -n +2 | awk '{ print NR, $1 }' | grep 'soy4' | sed 's/ .*//' > soy4may.keepers

#comm -13 SoyBjap_Corn_ddold_prepathms.txt soy4may.keepers > soy4may.keepers.uniq

#zcat features.tsv.gz | sed 's/soy4_//g' | gzip - > features1.tsv.gz
############################################################################################

# Loading of data

BjapMay.data <- Read10X(data.dir = "../Soy4.May5_filtered_feature_bc_matrix/oryg")
BjapMay <- CreateSeuratObject(counts = BjapMay.data, project = "Bjap_May")

BjapMay

# removing dublets found earlier via Python script

keepers_BjapMay <- try(read.table("./soy4may.keepers.uniq"))
if (class(keepers_BjapMay) == 'try-error') {
keepers_BjapMay=NULL
} else {
BjapMay <- subset(BjapMay, cells=WhichCells(BjapMay, keepers_BjapMay[,1]), invert=FALSE)
}
####
BjapMay
BjapMay.data.dd <- BjapMay[["RNA"]]@data

BjapMay.Mays.index <- grep(pattern = "^may5-", x = rownames(BjapMay.data.dd), value = FALSE)
percent.BjapMay.Mays <- Matrix::colSums(BjapMay.data.dd[BjapMay.Mays.index, ])/Matrix::colSums(BjapMay.data.dd)
BjapMay[["percent.May"]] <- percent.BjapMay.Mays
head(BjapMay@meta.data, 5)

VlnPlot(BjapMay, features = c("nFeature_RNA", "nCount_RNA", "percent.May"), ncol = 3)

BjapMay_30 <- subset(BjapMay, subset = percent.May < 0.30)

BjapMay_30.data.dd <- BjapMay_30[["RNA"]]@data
BjapMay_30.Mays.index <- grep(pattern = "^may5-", x = rownames(BjapMay_30.data.dd), value = FALSE)
percent.BjapMay_30.Mays <- Matrix::colSums(BjapMay_30.data.dd[BjapMay_30.Mays.index, ])/Matrix::colSums(BjapMay_30.data.dd)
BjapMay_30[["percent.May"]] <- percent.BjapMay_30.Mays
head(BjapMay_30@meta.data, 5)

BjapNoMay.data.dd <- BjapMay_30.data.dd[c(-BjapMay_30.Mays.index), ]
BjapNoMay <- CreateSeuratObject(counts = BjapNoMay.data.dd, project = "BjapNoMay", min.cells=5)
BjapNoMay

v1 <- VlnPlot(BjapNoMay, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
f1 <- FeatureScatter(BjapNoMay, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
v1 + f1

BjapNoMayFiltered <- subset(BjapNoMay, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)
BjapNoMayFiltered

BjapNoMayFiltered$protocol <- "BjapNoMayFiltered"
BjapNoMayFiltered$treatment <- "Bjap"
BjapNoMayFiltered <- NormalizeData(BjapNoMayFiltered)
BjapNoMayFiltered <- FindVariableFeatures(BjapNoMayFiltered, selection.method = "vst", nfeatures = 2000)

##################################################################################################

MockMedi.data <- Read10X(data.dir = "../Soy4.Medicago_filtered_feature_bc_matrix/oryg")
MockMedi <- CreateSeuratObject(counts = MockMedi.data, project = "Mock_Medi")

MockMedi

####
keepers_MockMedi <- try(read.table("./soy4med.keepers.uniq"))
if (class(keepers_MockMedi) == 'try-error') {
keepers_MockMedi=NULL
} else {
MockMedi <- subset(MockMedi, cells=WhichCells(MockMedi, keepers_MockMedi[,1]), invert=FALSE)
}
####
MockMedi

MockMedi.data.dd <- MockMedi[["RNA"]]@data

MockMedi.Medi.index <- grep(pattern = "^med5-", x = rownames(MockMedi.data.dd), value = FALSE)
percent.MockMedi.Medi <- Matrix::colSums(MockMedi.data.dd[MockMedi.Medi.index, ])/Matrix::colSums(MockMedi.data.dd)
MockMedi[["percent.Medi"]] <- percent.MockMedi.Medi
head(MockMedi@meta.data, 5)

VlnPlot(MockMedi, features = c("nFeature_RNA", "nCount_RNA", "percent.Medi"), ncol = 3)

MockMedi_50 <- subset(MockMedi, subset = percent.Medi < 0.50)

MockMedi_50.data.dd <- MockMedi_50[["RNA"]]@data
MockMedi_50.Medi.index <- grep(pattern = "^med5-", x = rownames(MockMedi_50.data.dd), value = FALSE)
percent.MockMedi.Medi <- Matrix::colSums(MockMedi_50.data.dd[MockMedi_50.Medi.index, ])/Matrix::colSums(MockMedi_50.data.dd)
MockMedi_50[["percent.Medi"]] <- percent.MockMedi_50.Medi
head(MockMedi_50@meta.data, 5)

MockNoMedi.data.dd <- MockMedi_50.data.dd[c(-MockMedi_50.Medi.index), ]
MockNoMedi <- CreateSeuratObject(counts = MockNoMedi.data.dd, project = "MockNoMedi", min.cells=5)
MockNoMedi

v2 <- VlnPlot(MockNoMedi, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
f2 <- FeatureScatter(MockNoMedi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
v2 + f2

MockNoMediFiltered <- subset(MockNoMedi, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 4000)
MockNoMediFiltered

MockNoMediFiltered$protocol <- "MockNoMediFiltered"
MockNoMediFiltered$treatment <- "Mock"
MockNoMediFiltered <- NormalizeData(MockNoMediFiltered)
MockNoMediFiltered <- FindVariableFeatures(MockNoMediFiltered, selection.method = "vst", nfeatures = 2000)

############################################################################################

B1.data <-Read10X(data.dir = "../SoyBjap1.v4_filtered_feature_bc_matrix")

B1 <- CreateSeuratObject(counts = B1.data, project = "B1")
B1

doublets_B1_old <- try(read.table("./B1_ddold_prepathms.txt"))
if (class(doublets_B1_old) == 'try-error') {
doublets_B1_old=NULL
} else {
B1 <- subset(B1, cells=WhichCells(B1, doublets_B1_old[,1]), invert=TRUE)
}

B1
B1.data.dd <- B1[["RNA"]]@data

B1.U.index <- grep(pattern = "^Glyma.U", x = rownames(B1.data.dd), value = FALSE)

percent.B1.U <- Matrix::colSums(B1.data.dd[B1.U.index, ])/Matrix::colSums(B1.data.dd)

B1 <-CreateSeuratObject(counts = B1.data.dd, project = "B1", min.cells=5) #Andrews
B1

B1[["percent.U"]] <- percent.B1.U

head(B1@meta.data, 5)

postscript('results-BM123_5_500_mod/B1.0.vnlplot.ps')
VlnPlot(B1, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)  #mine
dev.off()

postscript('results-BM123_5_500_mod/B1.0.featurescatter.ps')
plot1 <- FeatureScatter(B1, feature1 = "nCount_RNA", feature2 = "percent.U")  #mine
plot2 <- FeatureScatter(B1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #mine
plot1 + plot2
dev.off()

###modify below!!!!!!!!!!!!!!!

B1 <- subset(B1, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)  #modified
B1

B1$protocol <- "B1"
B1$treatment <- "Bjap"
B1 <- NormalizeData(B1)
B1 <- FindVariableFeatures(B1, selection.method = "vst", nfeatures = 2000)

################################################################################################
################################################################################################

B2.data <-Read10X(data.dir = "../SoyBjap2.v4_filtered_feature_bc_matrix")

B2 <- CreateSeuratObject(counts = B2.data, project = "B2")
B2

doublets_B2_old <- try(read.table("./B2_ddold_prepathms.txt"))
if (class(doublets_B2_old) == 'try-error') {
doublets_B2_old=NULL
} else {
B2 <- subset(B2, cells=WhichCells(B2, doublets_B2_old[,1]), invert=TRUE)
}

B2
B2.data.dd <- B2[["RNA"]]@data

B2.U.index <- grep(pattern = "^Glyma.U", x = rownames(B2.data.dd), value = FALSE)

percent.B2.U <- Matrix::colSums(B2.data.dd[B2.U.index, ])/Matrix::colSums(B2.data.dd)

B2 <-CreateSeuratObject(counts = B2.data.dd, project = "B2", min.cells=5) #Andrews
B2

B2[["percent.U"]] <- percent.B2.U

head(B2@meta.data, 5)

postscript('results-BM123_5_500_mod/B2.0.vnlplot.ps')
VlnPlot(B2, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)  #mine
dev.off()

postscript('results-BM123_5_500_mod/B2.0.featurescatter.ps')
plot1 <- FeatureScatter(B2, feature1 = "nCount_RNA", feature2 = "percent.U")  #mine
plot2 <- FeatureScatter(B2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #mine
plot1 + plot2
dev.off()

###modify below!!!!!!!!!!!!!!!

B2 <- subset(B2, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)  #modified
B2

B2$protocol <- "B2"
B2$treatment <- "Bjap"
B2 <- NormalizeData(B2)
B2 <- FindVariableFeatures(B2, selection.method = "vst", nfeatures = 2000)

################################################################################################
################################################################################################

B3.data <-Read10X(data.dir = "../SoyBjap3.v4_filtered_feature_bc_matrix")

B3 <- CreateSeuratObject(counts = B3.data, project = "B3")
B3
###
doublets_B3_old <- try(read.table("./B3_ddold_prepathms.txt"))
if (class(doublets_B3_old) == 'try-error') {
doublets_B3_old=NULL
} else {
B3 <- subset(B3, cells=WhichCells(B3, doublets_B3_old[,1]), invert=TRUE)
}

###
B3
B3.data.dd <- B3[["RNA"]]@data

B3.U.index <- grep(pattern = "^Glyma.U", x = rownames(B3.data.dd), value = FALSE)

percent.B3.U <- Matrix::colSums(B3.data.dd[B3.U.index, ])/Matrix::colSums(B3.data.dd)

B3 <-CreateSeuratObject(counts = B3.data.dd, project = "B3", min.cells=5) #Andrews
B3

B3[["percent.U"]] <- percent.B3.U

head(B3@meta.data, 5)

postscript('results-BM123_5_500_mod/B3.0.vnlplot.ps')
VlnPlot(B3, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)  #mine
dev.off()

postscript('results-BM123_5_500_mod/B3.0.featurescatter.ps')
plot1 <- FeatureScatter(B3, feature1 = "nCount_RNA", feature2 = "percent.U")  #mine
plot2 <- FeatureScatter(B3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #mine
plot1 + plot2
dev.off()

B3 <- subset(B3, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)  #modified
B3

B3$protocol <- "B3"
B3$treatment <- "Bjap"
B3 <- NormalizeData(B3)
B3 <- FindVariableFeatures(B3, selection.method = "vst", nfeatures = 2000)

#################################################################################################
#################################################################################################

M1.data <- Read10X(data.dir = "../SoyMock1.v4_filtered_feature_bc_matrix")

M1 <- CreateSeuratObject(counts = M1.data, project = "M1")
M1
###
doublets_M1_old <- try(read.table("./M1_ddold_prepathms.txt"))
if (class(doublets_M1_old) == 'try-error') {
doublets_M1_old=NULL
} else {
M1 <- subset(M1, cells=WhichCells(M1, doublets_M1_old[,1]), invert=TRUE)
}
###
M1
M1.data.dd <- M1[["RNA"]]@data

M1.U.index <- grep(pattern = "^Glyma.U", x = rownames(M1.data.dd), value = FALSE)

percent.M1.U <- Matrix::colSums(M1.data.dd[M1.U.index, ])/Matrix::colSums(M1.data.dd)

M1 <-CreateSeuratObject(counts = M1.data.dd, project = "M1", min.cells=5) #Andrews
M1

M1[["percent.U"]] <- percent.M1.U

head(M1@meta.data, 5)

postscript('results-BM123_5_500_mod/M1.0.vnlplot.ps')
VlnPlot(M1, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)  #mine
dev.off()

postscript('results-BM123_5_500_mod/M1.0.featurescatter.ps')
plot1 <- FeatureScatter(M1, feature1 = "nCount_RNA", feature2 = "percent.U")  #mine
plot2 <- FeatureScatter(M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #mine
plot1 + plot2
dev.off()

M1 <- subset(M1, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)  #modified
M1

M1$protocol <- "M1"
M1$treatment <- "Mock"
M1 <- NormalizeData(M1)
M1 <- FindVariableFeatures(M1, selection.method = "vst", nfeatures = 2000)


#################################################################################################
#################################################################################################

M2.data <- Read10X(data.dir = "../SoyMock2.v4_filtered_feature_bc_matrix")

M2 <- CreateSeuratObject(counts = M2.data, project = "M2")
M2
###
doublets_M2_old <- try(read.table("./M2_ddold_prepathms.txt"))
if (class(doublets_M2_old) == 'try-error') {
doublets_M2_old=NULL
} else {
M2 <- subset(M2, cells=WhichCells(M2, doublets_M2_old[,1]), invert=TRUE)
}
###
M2
M2.data.dd <- M2[["RNA"]]@data

M2.U.index <- grep(pattern = "^Glyma.U", x = rownames(M2.data.dd), value = FALSE)

percent.M2.U <- Matrix::colSums(M2.data.dd[M2.U.index, ])/Matrix::colSums(M2.data.dd)

M2 <-CreateSeuratObject(counts = M2.data.dd, project = "M2", min.cells=5) #Andrews
M2

M2[["percent.U"]] <- percent.M2.U

head(M2@meta.data, 5)

postscript('results-BM123_5_500_mod/M2.0.vnlplot.ps')
VlnPlot(M2, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)  #mine
dev.off()

postscript('results-BM123_5_500_mod/M2.0.featurescatter.ps')
plot1 <- FeatureScatter(M2, feature1 = "nCount_RNA", feature2 = "percent.U")  #mine
plot2 <- FeatureScatter(M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #mine
plot1 + plot2
dev.off()

M2 <- subset(M2, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)  #modified
M2

M2$protocol <- "M2"
M2$treatment <- "Mock"
M2 <- NormalizeData(M2)
M2 <- FindVariableFeatures(M2, selection.method = "vst", nfeatures = 2000)

#################################################################################################
#################################################################################################

M3.data <- Read10X(data.dir = "../SoyMock3.v4_filtered_feature_bc_matrix")

M3 <- CreateSeuratObject(counts = M3.data, project = "M3")
M3
###
doublets_M3_old <- try(read.table("./M3_ddold_prepathms.txt"))
if (class(doublets_M3_old) == 'try-error') {
doublets_M3_old=NULL
} else {
M3 <- subset(M3, cells=WhichCells(M3, doublets_M3_old[,1]), invert=TRUE)
}
###
M3
M3.data.dd <- M3[["RNA"]]@data

M3.U.index <- grep(pattern = "^Glyma.U", x = rownames(M3.data.dd), value = FALSE)

percent.M3.U <- Matrix::colSums(M3.data.dd[M3.U.index, ])/Matrix::colSums(M3.data.dd)

M3 <-CreateSeuratObject(counts = M3.data.dd, project = "M3", min.cells=5) #Andrews
M3

M3[["percent.U"]] <- percent.M3.U

head(M3@meta.data, 5)

postscript('results-BM123_5_500_mod/M3.0.vnlplot.ps')
VlnPlot(M3, features = c("nFeature_RNA", "nCount_RNA", "percent.U"), ncol = 3)  #mine
dev.off()

postscript('results-BM123_5_500_mod/M3.0.featurescatter.ps')
plot1 <- FeatureScatter(M3, feature1 = "nCount_RNA", feature2 = "percent.U")  #mine
plot2 <- FeatureScatter(M3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #mine
plot1 + plot2
dev.off()

M3 <- subset(M3, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)  #mine
M3

M3$protocol <- "M3"
M3$treatment <- "Mock"
M3 <- NormalizeData(M3)
M3 <- FindVariableFeatures(M3, selection.method = "vst", nfeatures = 2000)

#################################################################################################
################################################################################################
#####################################################################################################
protocol.anchors <- FindIntegrationAnchors(object.list = list(BjapNoMayFiltered,MockNoMediFiltered,B1,B2,B3,M1,M2,M3), dims = 1:30)

protocol.combined <- IntegrateData(anchorset = protocol.anchors, dims = 1:30)

DefaultAssay(protocol.combined) <- "integrated"

protocol.combined <- ScaleData(protocol.combined, verbose = FALSE)
protocol.combined <- RunPCA(protocol.combined, npcs = 50, verbose = FALSE)

#DimHeatmap(protocol.combined, dims = 1:50, cells = 500, balanced = TRUE)

ElbowPlot(protocol.combined, ndims=50)

#########check no of PCs

protocol.combined <- RunUMAP(protocol.combined, reduction = "pca", dims = 1:50)

protocol.combined <- FindNeighbors(protocol.combined, reduction = "pca", dims = 1:50) #Andrews
##################################################################################################
DefaultAssay(protocol.combined) <- "integrated"

protocol.combined <- FindClusters(protocol.combined, resolution = 0.5) #Andrews

pdf("results-BM123barnyard_5_500_mod/groupdimplot-umap_05.pdf")
DimPlot(protocol.combined, reduction = "umap", group.by = "treatment")
dev.off()

pdf("results-BM123barnyard_5_500_mod/splitdimplotgroup-umap_05.pdf")
DimPlot(protocol.combined, reduction = "umap", split.by = "treatment", label = TRUE, label.size=5)
dev.off()

pdf("results-BM123barnyard_5_500_mod/dimplotclusters-umap_05.pdf")
DimPlot(protocol.combined, reduction = "umap", label = TRUE)
dev.off()

pdf("results-BM123barnyard_5_500_mod/dimplotclustersNolabel-umap_05.pdf")
DimPlot(protocol.combined, reduction = "umap")
dev.off()

pdf("results-BM123barnyard_5_500_mod/gsoy4-Glyma.02G149100-umap_05.pdf")
FeaturePlot(protocol.combined, features = "Glyma.02G149100")
dev.off()

write.table(table(protocol.combined@meta.data$seurat_clusters, protocol.combined@meta.data$orig.ident), "results-BM123barnyard_5_500_mod/cluster_counts_05.txt", sep="\t", quote=F)

write.table(prop.table(table(protocol.combined@meta.data$seurat_clusters, protocol.combined@meta.data$orig.ident),2), "results-BM123barnyard_5_500_mod/cluster_proportions.by_sample_05.txt", sep="\t", quote=F)

saveRDS(protocol.combined, "results-BM123barnyard_5_500_mod/protocol.combined_05.rds")

write.table(FetchData(protocol.combined, "ident"), "results-BM123barnyard_5_500_mod/seurat_bc_clustermap_05")

write.csv(protocol.combined@reductions$umap@cell.embeddings, file = "results-BM123barnyard_5_500_mod/umap_05.csv")

################################################
DefaultAssay(protocol.combined) <- "RNA"

pdf("results-BM123barnyard_5_500_mod/gsoy4-Glyma.02G149100dimplot-umap_05.pdf")
f <- FeaturePlot(protocol.combined, features = "Glyma.02G149100")
d <- DimPlot(protocol.combined, reduction = "umap", split.by = "treatment")
d + f
dev.off()

pdf("results-BM123barnyard_5_500_mod/soybeanroothairresponse-umap_05.pdf")
FeaturePlot(protocol.combined, features = c("Glyma.09G187000","Glyma.11G244800","Glyma.20G153400","Glyma.01G174400","Glyma.17G066800","Glyma.10G177100","Glyma.05G051400","Glyma.08G120100"))
dev.off()

pdf("results-BM123barnyard_5_500_mod/arabidopsisroothair-umap_05.pdf")
FeaturePlot(protocol.combined, features = c("Glyma.01G214400","Glyma.05G051200","Glyma.11G027600","Glyma.17G133100","Glyma.09G047700","Glyma.15G155100"))
dev.off()

pdf("results-BM123barnyard_5_500_mod/soybeanroothair-umap_05.pdf")
FeaturePlot(protocol.combined, features = c("Glyma.18G025200","Glyma.05G157400","Glyma.15G020700","Glyma.08G115000","Glyma.03G188300"))
dev.off()

pdf("results-BM123barnyard_5_500_mod/soybeanroothairresponseMAIN-umap_05.pdf")
FeaturePlot(protocol.combined, features = c("Glyma.11G244800","Glyma.20G153400","Glyma.01G174400","Glyma.05G051400","Glyma.08G120100"))
dev.off()

FeaturePlot(protocol.combined, features = c("Glyma.11G244800","Glyma.08G120100"))

FeaturePlot(protocol.combined, features = c("Glyma.02G149100","Glyma.19G065600","Glyma.10G024600"))

pdf("results-BM123barnyard_5_500_mod/UMIperNuclei-umap_05.pdf")
FeaturePlot(allSoy, features = "nCount_RNA", cols = c("yellow", "blue"), label=TRUE, label.size=2)
dev.off()

pdf("results-BM123barnyard_5_500_mod/GenesperNuclei-umap_05.pdf")
FeaturePlot(allSoy, features = "nFeature_RNA", cols = c("yellow", "blue"), label=TRUE, label.size=2)
dev.off()

pdf("results-BM123barnyard_5_500_mod/UMIandGeneperNuclei-umap_05.pdf")
u <- FeaturePlot(allSoy, features = "nCount_RNA", cols = c("yellow", "blue"), label=TRUE, label.size=2)
g <- FeaturePlot(allSoy, features = "nFeature_RNA", cols = c("yellow", "blue"), label=TRUE, label.size=2)
u + g
dev.off()

########################################################################################################
########################################################################################################
#########Finding Conserved Markers
###still for 
DefaultAssay(protocol.combined) <- "RNA"

c0.markers <- FindConservedMarkers(protocol.combined, ident.1 = 0, grouping.var = "treatment", verbose = FALSE)
write.table(c0.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c0_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c1.markers <- FindConservedMarkers(protocol.combined, ident.1 = 1, grouping.var = "treatment", verbose = FALSE)
write.table(c1.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c1_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c2.markers <- FindConservedMarkers(protocol.combined, ident.1 = 2, grouping.var = "treatment", verbose = FALSE)
write.table(c2.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c2_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c3.markers <- FindConservedMarkers(protocol.combined, ident.1 = 3, grouping.var = "treatment", verbose = FALSE)
write.table(c3.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c3_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c4.markers <- FindConservedMarkers(protocol.combined, ident.1 = 4, grouping.var = "treatment", verbose = FALSE)
write.table(c4.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c4_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c5.markers <- FindConservedMarkers(protocol.combined, ident.1 = 5, grouping.var = "treatment", verbose = FALSE)
write.table(c5.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c5_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c6.markers <- FindConservedMarkers(protocol.combined, ident.1 = 6, grouping.var = "treatment", verbose = FALSE)
write.table(c6.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c6_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c7.markers <- FindConservedMarkers(protocol.combined, ident.1 = 7, grouping.var = "treatment", verbose = FALSE)
write.table(c7.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c7_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c8.markers <- FindConservedMarkers(protocol.combined, ident.1 = 8, grouping.var = "treatment", verbose = FALSE)
write.table(c8.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c8_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c9.markers <- FindConservedMarkers(protocol.combined, ident.1 = 9, grouping.var = "treatment", verbose = FALSE)
write.table(c9.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c9_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c10.markers <- FindConservedMarkers(protocol.combined, ident.1 = 10, grouping.var = "treatment", verbose = FALSE)
write.table(c10.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c10_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c11.markers <- FindConservedMarkers(protocol.combined, ident.1 = 11, grouping.var = "treatment", verbose = FALSE)
write.table(c11.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c11_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c12.markers <- FindConservedMarkers(protocol.combined, ident.1 = 12, grouping.var = "treatment", verbose = FALSE)
write.table(c12.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c12_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c13.markers <- FindConservedMarkers(protocol.combined, ident.1 = 13, grouping.var = "treatment", verbose = FALSE)
write.table(c13.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c13_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c14.markers <- FindConservedMarkers(protocol.combined, ident.1 = 14, grouping.var = "treatment", verbose = FALSE)
write.table(c14.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c14_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c15.markers <- FindConservedMarkers(protocol.combined, ident.1 = 15, grouping.var = "treatment", verbose = FALSE)
write.table(c15.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c15_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c16.markers <- FindConservedMarkers(protocol.combined, ident.1 = 16, grouping.var = "treatment", verbose = FALSE)
write.table(c16.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c16_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c17.markers <- FindConservedMarkers(protocol.combined, ident.1 = 17, grouping.var = "treatment", verbose = FALSE)
write.table(c17.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c17_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c18.markers <- FindConservedMarkers(protocol.combined, ident.1 = 18, grouping.var = "treatment", verbose = FALSE)
write.table(c18.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c18_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c19.markers <- FindConservedMarkers(protocol.combined, ident.1 = 19, grouping.var = "treatment", verbose = FALSE)
write.table(c19.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c19_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c20.markers <- FindConservedMarkers(protocol.combined, ident.1 = 20, grouping.var = "treatment", verbose = FALSE)
write.table(c20.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c20_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c21.markers <- FindConservedMarkers(protocol.combined, ident.1 = 21, grouping.var = "treatment", verbose = FALSE)
write.table(c21.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c21_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c22.markers <- FindConservedMarkers(protocol.combined, ident.1 = 22, grouping.var = "treatment", verbose = FALSE)
write.table(c22.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c22_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)

c23.markers <- FindConservedMarkers(protocol.combined, ident.1 = 23, grouping.var = "treatment", verbose = FALSE)
write.table(c23.markers, "results-BM123barnyard_5_500_mod/findconservedmarkers_c23_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)


###################################################################################
############################################
c12 - Meristem
FeaturePlot(protocol.combined, features = c("Glyma.12G018900","Glyma.08G021700","Glyma.15G276000","Glyma.05G215500","Glyma.03G123600","Glyma.02G135700","Glyma.06G290300","Glyma.01G229600","Glyma.12G116700","Glyma.02G214300","Glyma.02G277600","Glyma.19G127200"))







#######################################

#ls -1 findconservedmarkers_c*_05.txt | while read file; do echo $file | sed 's/.*c\(.*\)_.*/\1/' > tmp; cluster=$(cat tmp); cat $file | tail -n +2 | head -n 1 | cut -f1 > tmp2; marker=$(cat tmp2); echo "'$cluster'" = "$marker"; done > conservedmarkers2clusters.txt

library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)

protocol.combined <- readRDS("results-BM123barnyard_5_500_mod/protocol.combined_05.rds")

DefaultAssay(protocol.combined) <- "RNA"

pdf("results-BM123barnyard_5_500_mod/conservedmarkers-feature_05.pdf")
FeaturePlot(protocol.combined, features = c("Glyma.09G047200","Glyma.18G054000","Glyma.14G145900","Glyma.13G091600","Glyma.07G222500","Glyma.18G043200","Glyma.10G185800","Glyma.17G180400","Glyma.10G028300","Glyma.15G274600","Glyma.09G131500","Glyma.20G151700","Glyma.18G078900","Glyma.15G008500","Glyma.06G259400","Glyma.08G125800","Glyma.07G129900","Glyma.15G129100","Glyma.02G080900","Glyma.16G106800","Glyma.08G021700"))
dev.off()

#protocol.combined <- RenameIdents(protocol.combined, '0' = "09G047200",'10' = "18G054000",'1' = "14G145900",'11' = "13G091600",'12' = "07G222500",'13' = "18G043200",'14' = "10G185800",'15' = "17G180400",'16' = "10G028300",'17' = "15G274600",'18' = "09G131500",'19' = "20G151700",'20' = "18G078900",'2' = "15G008500",'3' = "06G259400",'4' = "08G125800",'5' = "07G129900",'6' = "15G129100",'7' = "02G080900",'8' = "16G106800",'9' = "08G021700")

pdf("results-BM123barnyard_5_500_mod/conservedmarkers-umap_05.pdf")
DimPlot(protocol.combined, label = TRUE) + NoLegend()
dev.off()

pdf("results-BM123barnyard_5_500_mod/roothairmarkers-feature_05.pdf")
FeaturePlot(protocol.combined, features = c("Glyma.01G214400","Glyma.05G051200","Glyma.09G047700","Glyma.11G027600","Glyma.15G155100","Glyma.17G133100"))
dev.off()

protocol.combined$celltype.treatment <- paste(Idents(protocol.combined), protocol.combined$treatment, sep = "_")
protocol.combined$celltype <- Idents(protocol.combined)
Idents(protocol.combined) <- "celltype.treatment"

treatment.response <- FindMarkers(protocol.combined, ident.1 = "7_Bjap", ident.2 = "7_Mock", verbose = FALSE)
head(treatment.response, n = 15)
write.table(treatment.response, "results-BM123barnyard_5_500_mod/treatment.response_roothair_05.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
###repeat
pdf("results-BM123barnyard_5_500_mod/roothairdiffgenes-feature_05.pdf")
FeaturePlot(protocol.combined, features = c("Glyma.07G007102","Glyma.11G031400","Glyma.01G186600"), split.by = "treatment", max.cutoff = 3, cols = c("grey", "red"))
dev.off()

FeaturePlot(protocol.combined, features = c("Glyma.03G054100","Glyma.19G169400","Glyma.13G266100"), split.by = "treatment", max.cutoff = 3, cols = c("grey", "red"))


pdf("results-BM123barnyard_5_500_mod/roothairdiffgenes-vlnplot_05.pdf")
plots <- VlnPlot(protocol.combined, features = c("Glyma.07G007102","Glyma.06G324300","Glyma.11G031400","Glyma.08G238100","Glyma.01G186600"), split.by = "treatment", group.by = "celltype", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

#commonv1v2
pdf("results-BM123barnyard_5_500_mod/roothairdiffgenes-vlnplot_05.pdf")
plots <- VlnPlot(protocol.combined, features = c("Glyma.01G186600","Glyma.03G054100","Glyma.04G076000","Glyma.13G170600","Glyma.13G266100","Glyma.18G226500","Glyma.19G169400"), split.by = "treatment", group.by = "celltype", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

###################################################################################



