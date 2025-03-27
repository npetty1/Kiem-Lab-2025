#load relevant packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
#load in datasets
WBC1.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_wbc_1_Oct2024/outs/filtered_feature_bc_matrix")
WBC2.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_wbc_2_Oct2024/outs/filtered_feature_bc_matrix")
HSC1.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_hsc_1_Oct2024/outs/filtered_feature_bc_matrix")
HSC2.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_hsc_2_Oct2024/outs/filtered_feature_bc_matrix")
Cer1.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_cerebellum_1_Oct2024/outs/filtered_feature_bc_matrix")
Cer2.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_cerebellum_2_Oct2024/outs/filtered_feature_bc_matrix")
CTX1.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_frontal_cortex_1_Oct2024/outs/filtered_feature_bc_matrix")
CTX2.data <- Read10X(data.dir = "~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/Z21226_frontal_cortex_2_Oct2024/outs/filtered_feature_bc_matrix")
#Create Seurat objects
WBC1 <- CreateSeuratObject(counts = WBC1.data, project = "Z21226.WBC1", min.cells = 3, min.features = 200)
WBC2 <- CreateSeuratObject(counts = WBC2.data, project = "Z21226.WBC2", min.cells = 3, min.features = 200)
HSC1 <- CreateSeuratObject(counts = HSC1.data, project = "Z21226.HSC1", min.cells = 3, min.features = 200)
HSC2 <- CreateSeuratObject(counts = HSC2.data, project = "Z21226.HSC2", min.cells = 3, min.features = 200)
Cer1 <- CreateSeuratObject(counts = Cer1.data, project = "Z21226.Cer1", min.cells = 3, min.features = 200)
Cer2 <- CreateSeuratObject(counts = Cer2.data, project = "Z21226.Cer2", min.cells = 3, min.features = 200)
CTX1 <- CreateSeuratObject(counts = CTX1.data, project = "Z21226.CTX1", min.cells = 3, min.features = 200)
CTX2 <- CreateSeuratObject(counts = CTX2.data, project = "Z21226.CTX2", min.cells = 3, min.features = 200)

#Analyze WBC1
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%WBC1@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data quality
WBC1[["percent.mt"]] <- PercentageFeatureSet(WBC1, features = mito.gene.list)
VlnPlot(WBC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(WBC1)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(WBC1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WBC1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#normalize the data
WBC1 <- subset(WBC1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
#After filtering, normalize data and find top 10 variable features
WBC1 <- NormalizeData(WBC1, normalization.method = "LogNormalize", scale.factor = 10000)
WBC1 <- FindVariableFeatures(WBC1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(WBC1), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(WBC1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes <- rownames(WBC1)
WBC1 <- ScaleData(WBC1, features = all.genes)
#dimensional reduction
WBC1 <- RunPCA(WBC1, features = VariableFeatures(object = WBC1))
#Visualization of principal components
print(WBC1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WBC1, dims = 1:2, reduction = "pca")
DimPlot(WBC1, reduction = "pca") + NoLegend()
DimHeatmap(WBC1, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(WBC1)
#Perform clustering
WBC1 <- FindNeighbors(WBC1, dims = 1:15)
WBC1 <- FindClusters(WBC1, resolution = 0.3)
#Calculate UMAP clustering
WBC1 <- RunUMAP(WBC1, dims = 1:15)
DimPlot(WBC1, reduction = "umap")
#Fina all markers for each cluster
cluster0.markers <- FindMarkers(WBC1, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(WBC1, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(WBC1, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(WBC1, ident.1 = 3)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(WBC1, ident.1 = 4)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(WBC1, ident.1 = 5)
head(cluster5.markers, n = 5)
cluster6.markers <- FindMarkers(WBC1, ident.1 = 6)
head(cluster6.markers, n = 5)
cluster7.markers <- FindMarkers(WBC1, ident.1 = 7)
head(cluster7.markers, n = 5)
#Violin plots of differentially expressed genes
VlnPlot(WBC1, features = c("CD79A", "CD14","CD8A","CD74","NKG7"))
#Featureplot of differentially expressed genes
FeaturePlot(WBC1, features = c("CD79A", "CD14","CD8A","CD63","CD4"))
#Determine most differential features
WBC.markers <- FindAllMarkers(WBC1, only.pos = TRUE)
WBC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
WBC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(WBC1, features = top10$gene) + NoLegend()

#Analyze WBC2
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%WBC2@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data quality
WBC2[["percent.mt"]] <- PercentageFeatureSet(WBC2, features = mito.gene.list)
VlnPlot(WBC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(WBC2)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(WBC2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WBC2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#normalize the data
WBC2 <- subset(WBC2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
#After filtering, normalize data and find top 10 variable features
WBC2 <- NormalizeData(WBC2, normalization.method = "LogNormalize", scale.factor = 10000)
WBC2 <- FindVariableFeatures(WBC2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(WBC2), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(WBC2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes <- rownames(WBC2)
WBC2 <- ScaleData(WBC2, features = all.genes)
#dimensional reduction
WBC2 <- RunPCA(WBC2, features = VariableFeatures(object = WBC2))
#Visualization of principal components
print(WBC2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WBC2, dims = 1:2, reduction = "pca")
DimPlot(WBC2, reduction = "pca") + NoLegend()
DimHeatmap(WBC2, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(WBC2)
#Perform clustering
WBC2 <- FindNeighbors(WBC2, dims = 1:15)
WBC2 <- FindClusters(WBC2, resolution = 0.3)
#Calculate UMAP clustering
WBC2 <- RunUMAP(WBC2, dims = 1:15)
DimPlot(WBC2, reduction = "umap")
FeaturePlot(WBC2,reduction="umap",features="CD14")
#Fina all markers for each cluster
cluster0.markers <- FindMarkers(WBC2, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(WBC2, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(WBC2, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(WBC2, ident.1 = 3)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(WBC2, ident.1 = 4)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(WBC2, ident.1 = 5)
head(cluster5.markers, n = 5)
cluster6.markers <- FindMarkers(WBC2, ident.1 = 6)
head(cluster6.markers, n = 5)
cluster7.markers <- FindMarkers(WBC2, ident.1 = 7)
head(cluster7.markers, n = 5)
#Violin plots of differentially expressed genes
VlnPlot(WBC1, features = c("EGFP"))
#Featureplot of differentially expressed genes
FeaturePlot(WBC1, features = c("EGFP"))
#Determine most differential features
WBC.markers <- FindAllMarkers(WBC2, only.pos = TRUE)
WBC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
WBC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(WBC2, features = top10$gene) + NoLegend()


#Analyzing Cerebellum 1
#Define mitochondrial genes
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%Cer1@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data quality
Cer1[["percent.mt"]] <- PercentageFeatureSet(Cer1, features = mito.gene.list)
VlnPlot(Cer1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(Cer1)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(Cer1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cer1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter the data
Cer1 <- subset(Cer1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
#After filtering, normalize data and find top 10 variable features
Cer1 <- NormalizeData(Cer1, normalization.method = "LogNormalize", scale.factor = 10000)
Cer1 <- FindVariableFeatures(Cer1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Cer1), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(Cer1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes.Cer1 <- rownames(Cer1)
Cer1 <- ScaleData(Cer1, features = all.genes.Cer1)
#dimensional reduction
Cer1 <- RunPCA(Cer1, features = VariableFeatures(object = Cer1))
#Visualization of principal components
print(Cer1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Cer1, dims = 1:2, reduction = "pca")
DimPlot(Cer1, reduction = "pca") + NoLegend()
DimHeatmap(Cer1, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(Cer1)
#Perform clustering
Cer1 <- FindNeighbors(Cer1, dims = 1:15)
Cer1 <- FindClusters(Cer1, resolution = 0.3)
#Calculate UMAP clustering
Cer1 <- RunUMAP(Cer1, dims = 1:15)
DimPlot(Cer1, reduction = "umap")
#Fina all markers for each cluster
cluster0.markers <- FindMarkers(Cer1, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(Cer1, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(Cer1, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(Cer1, ident.1 = 3)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(Cer1, ident.1 = 4)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(Cer1, ident.1 = 5)
head(cluster5.markers, n = 5)

#Analyzing Cerebellum 2
#Define mitochondrial genes
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%Cer2@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data quality
Cer2[["percent.mt"]] <- PercentageFeatureSet(Cer2, features = mito.gene.list)
VlnPlot(Cer2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(Cer2)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(Cer2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cer2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter the data
Cer2 <- subset(Cer2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#After filtering, normalize data and find top 10 variable features
Cer2 <- NormalizeData(Cer2, normalization.method = "LogNormalize", scale.factor = 10000)
Cer2 <- FindVariableFeatures(Cer2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Cer2), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(Cer2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes <- rownames(Cer2)
Cer2 <- ScaleData(Cer2, features = all.genes)
#dimensional reduction
Cer2 <- RunPCA(Cer2, features = VariableFeatures(object = Cer2))
#Visualization of principal components
print(Cer2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Cer2, dims = 1:2, reduction = "pca")
DimPlot(Cer2, reduction = "pca") + NoLegend()
DimHeatmap(Cer2, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(Cer2)
#Perform clustering
Cer2 <- FindNeighbors(Cer2, dims = 1:15)
Cer2 <- FindClusters(Cer2, resolution = 0.3)
#Calculate UMAP clustering
Cer2 <- RunUMAP(Cer2, dims = 1:15)
DimPlot(Cer2, reduction = "umap")
#Fina all markers for each cluster
cluster0.markers <- FindMarkers(Cer2, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(Cer2, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(Cer2, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(Cer2, ident.1 = 3)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(Cer2, ident.1 = 4)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(Cer2, ident.1 = 5)
head(cluster5.markers, n = 5)
#Violin plots of differentially expressed genes
VlnPlot(Cer2, features = c("EGFP"))
#Featureplot of differentially expressed genes
FeaturePlot(Cer2, features = c("P2RY12","CX3CR1","CSF1R"))


#Analyzing Frontal Cortex 1
#Define mitochondrial genes
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%CTX1@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data qualitCTX1#Add mitochondria as a variable and visualize data quality
CTX1[["percent.mt"]] <- PercentageFeatureSet(CTX1, features = mito.gene.list)
VlnPlot(CTX1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(CTX1)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(CTX1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CTX1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter the data
CTX1 <- subset(CTX1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 2.5)
#After filtering, normalize data and find top 10 variable features
CTX1 <- NormalizeData(CTX1, normalization.method = "LogNormalize", scale.factor = 10000)
CTX1 <- FindVariableFeatures(CTX1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(CTX1), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(CTX1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes <- rownames(CTX1)
CTX1 <- ScaleData(CTX1, features = all.genes)
#dimensional reduction
CTX1 <- RunPCA(CTX1, features = VariableFeatures(object = CTX1))
#Visualization of principal components
print(CTX1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CTX1, dims = 1:2, reduction = "pca")
DimPlot(CTX1, reduction = "pca") + NoLegend()
DimHeatmap(CTX1, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(CTX1)
#Perform clustering
CTX1 <- FindNeighbors(CTX1, dims = 1:15)
CTX1 <- FindClusters(CTX1, resolution = 0.3)
#Calculate UMAP clustering
CTX1 <- RunUMAP(CTX1, dims = 1:15)
DimPlot(CTX1, reduction = "umap")
#Fina all markers for each cluster
cluster0.markers <- FindMarkers(CTX1, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(CTX1, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(CTX1, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(CTX1, ident.1 = 3)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(CTX1, ident.1 = 4)
head(cluster4.markers, n = 5)

#Violin plots of differentially expressed genes
VlnPlot(CTX1, features = c("EGFP"))
#Featureplot of differentially expressed genes
FeaturePlot(CTX1, features = c("P2RY12","CX3CR1","CSF1R"))

#Analyzing Frontal Cortex 2
#Define mitochondrial genes
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%CTX2@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data qualitCTX1#Add mitochondria as a variable and visualize data quality
CTX2[["percent.mt"]] <- PercentageFeatureSet(CTX2, features = mito.gene.list)
VlnPlot(CTX2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(CTX2)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(CTX2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CTX2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter the data
CTX2 <- subset(CTX2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
#After filtering, normalize data and find top 10 variable features
CTX2 <- NormalizeData(CTX2, normalization.method = "LogNormalize", scale.factor = 10000)
CTX2 <- FindVariableFeatures(CTX2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(CTX2), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(CTX2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes <- rownames(CTX2)
CTX2 <- ScaleData(CTX2, features = all.genes)
#dimensional reduction
CTX2 <- RunPCA(CTX2, features = VariableFeatures(object = CTX2))
#Visualization of principal components
print(CTX2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CTX2, dims = 1:2, reduction = "pca")
DimPlot(CTX2, reduction = "pca") + NoLegend()
DimHeatmap(CTX2, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(CTX2)
#Perform clustering
CTX2 <- FindNeighbors(CTX2, dims = 1:15)
CTX2 <- FindClusters(CTX2, resolution = 0.3)
#Calculate UMAP clustering
CTX2 <- RunUMAP(CTX2, dims = 1:15)
DimPlot(CTX2, reduction = "umap")
#Finding all markers for each cluster
cluster0.markers <- FindMarkers(CTX2, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(CTX2, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(CTX2, ident.1 = 2)
head(cluster2.markers, n=10)

#Violin plots of differentially expressed genes
VlnPlot(CTX2, features = c("EGFP"))
#Featureplot of differentially expressed genes
FeaturePlot(CTX2, features = c("P2RY12","CX3CR1","CSF1R"))


#Analyzing HSC 1
#Define mitochondrial genes
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%HSC1@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data qualitCTX1#Add mitochondria as a variable and visualize data quality
HSC1[["percent.mt"]] <- PercentageFeatureSet(HSC1, features = mito.gene.list)
VlnPlot(HSC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(HSC1)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(HSC1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HSC1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter the data
HSC1 <- subset(HSC1, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
#After filtering, normalize data and find top 10 variable features
HSC1 <- NormalizeData(HSC1, normalization.method = "LogNormalize", scale.factor = 10000)
HSC1 <- FindVariableFeatures(HSC1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(HSC1), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(HSC1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes <- rownames(HSC1)
HSC1 <- ScaleData(HSC1, features = all.genes)
#dimensional reduction
HSC1 <- RunPCA(HSC1, features = VariableFeatures(object = HSC1))
#Visualization of principal components
print(HSC1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HSC1, dims = 1:2, reduction = "pca")
DimPlot(HSC1, reduction = "pca") + NoLegend()
DimHeatmap(HSC1, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(HSC1)
#Perform clustering
HSC1 <- FindNeighbors(HSC1, dims = 1:20)
HSC1 <- FindClusters(HSC1, resolution = 0.3)
#Calculate UMAP clustering
HSC1 <- RunUMAP(HSC1, dims = 1:20)
DimPlot(HSC1, reduction = "umap")
#Fina all markers for each cluster
cluster0.markers <- FindMarkers(HSC1, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(HSC1, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(HSC1, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(HSC1, ident.1 = 3)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(HSC1, ident.1 = 4)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(HSC1, ident.1 = 5)
head(cluster5.markers, n = 5)
cluster6.markers <- FindMarkers(HSC1, ident.1 = 6)
head(cluster6.markers, n = 5)
cluster7.markers <- FindMarkers(HSC1, ident.1 = 7)
head(cluster7.markers, n = 10)
cluster8.markers <- FindMarkers(HSC1, ident.1 = 8)
head(cluster8.markers, n = 10)

#Violin plots of differentially expressed genes
VlnPlot(HSC1, features = c("EGFP"))
#Featureplot of differentially expressed genes
FeaturePlot(HSC1, features = c("P2RY12","CX3CR1","CSF1R"))

#Analyzing HSC 2
#Define mitochondrial genes
mito.gene.list <- read.csv(paste0("~/Desktop/Kiem Lab/Cell Stem Cell Single Cell/NHP_mitochondrial_genes.csv"), header=F)
mito.gene.list <- mito.gene.list$V1
mito.gene.list <- mito.gene.list[which(mito.gene.list%in%HSC2@assays$RNA$counts@Dimnames[[1]])]
#Add mitochondria as a variable and visualize data qualitCTX1#Add mitochondria as a variable and visualize data quality
HSC2[["percent.mt"]] <- PercentageFeatureSet(HSC2, features = mito.gene.list)
VlnPlot(HSC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(HSC2)
#plot feature to determine appropirate thresholds 
plot1 <- FeatureScatter(HSC2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HSC2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter the data
HSC2 <- subset(HSC2, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
#After filtering, normalize data and find top 10 variable features
HSC2 <- NormalizeData(HSC2, normalization.method = "LogNormalize", scale.factor = 10000)
HSC2 <- FindVariableFeatures(HSC2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(HSC2), 10)
#Plot features and annotate top 10 variable features
plot1 <- VariableFeaturePlot(HSC2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scale data
all.genes <- rownames(HSC2)
HSC2 <- ScaleData(HSC2, features = all.genes)
#dimensional reduction
HSC2 <- RunPCA(HSC2, features = VariableFeatures(object = HSC2))
#Visualization of principal components
print(HSC2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HSC2, dims = 1:2, reduction = "pca")
DimPlot(HSC2, reduction = "pca") + NoLegend()
DimHeatmap(HSC2, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(HSC2)
#Perform clustering
HSC2 <- FindNeighbors(HSC2, dims = 1:20)
HSC2 <- FindClusters(HSC2, resolution = 0.3)
#Calculate UMAP clustering
HSC2 <- RunUMAP(HSC2, dims = 1:20)
DimPlot(HSC2, reduction = "umap")
#Fina all markers for each cluster
cluster0.markers <- FindMarkers(HSC2, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(HSC2, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(HSC2, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(HSC2, ident.1 = 3)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(HSC2, ident.1 = 4)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(HSC2, ident.1 = 5)
head(cluster5.markers, n = 5)
cluster6.markers <- FindMarkers(HSC2, ident.1 = 6)
head(cluster6.markers, n = 5)

#Violin plots of differentially expressed genes
VlnPlot(HSC2, features = c("EGFP"))
#Featureplot of differentially expressed genes
FeaturePlot(HSC2, features = c("EGFP"))

#Merge brain and WBC datasets
brain.wbc.merged <-merge(CTX1, y=c(CTX2, Cer1, Cer2, WBC1, WBC2), add.cell.ids=c("CTX1", "CTX2", "Cer1", "Cer2", "WBC1","WBC2"), project="Z21226", merge.data=TRUE)
brain.wbc.merged <-NormalizeData(brain.wbc.merged)
brain.wbc.merged <-FindVariableFeatures(brain.wbc.merged)
brain.wbc.merged <-ScaleData(brain.wbc.merged)
brain.wbc.merged <-RunPCA(brain.wbc.merged)
brain.wbc.merged <-FindNeighbors(brain.wbc.merged, dims=1:30, reduction="pca")
brain.wbc.merged <-FindClusters(brain.wbc.merged, resolution = 2, cluster.name = "unintegrated clusters")
brain.wbc.merged <-RunUMAP(brain.wbc.merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(brain.wbc.merged, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
DoHeatmap(brain.wbc.merged, features = top10$gene) + NoLegend()
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="EGFP")
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="CD4")
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="CD8A")
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="CD79A")
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="NKG7")
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="AIF1")
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="P2RY12")
FeaturePlot(brain.wbc.merged, reduction = "umap.unintegrated", features="CSF1R")

#Merge all datasets
all.merged <- merge(CTX1, y = c(CTX2, Cer1, Cer2), add.cell.ids = c("CTX1", "CTX2", "Cer1", "Cer2"), project = "Z21226", merge.data=TRUE)
#normalize merged datasets, scale, and calculate PCA
all.merged <- NormalizeData(all.merged)
all.merged <- FindVariableFeatures(all.merged)
all.merged <- ScaleData(all.merged)
all.merged <- RunPCA(all.merged)
#Visalized merged, non-integrated data
all.merged <- FindNeighbors(all.merged, dims = 1:30, reduction = "pca")
all.merged <- FindClusters(all.merged, resolution = 2, cluster.name = "unintegrated_clusters")
all.merged <- RunUMAP(all.merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by sample and cluster
DimPlot(all.merged, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
#integrate samples
all.integrated <- IntegrateLayers(
  object=all.merged, method = CCAIntegration,
  new.reduction = "integrated.cca",
  verbose = FALSE
)
all.integrated <- JoinLayers(all.integrated)
#Process & visualize integrated data
all.integrated <- FindNeighbors(all.integrated, reduction = "integrated.cca", dims = 1:30)
all.integrated <- FindClusters(all.integrated, resolution = 0.1, cluster.name = "cca_clusters")
all.integrated <- RunUMAP(all.integrated, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(
  all.integrated,
  reduction = "umap.cca",
  group.by = c("orig.ident", "cca_clusters"),
  combine = FALSE, label.size = 2
)
#Show newly merged microglia data
DimPlot(all.integrated, reduction = "umap.cca")
#Assess GFP by cluster
FeaturePlot(all.integrated, reduction = "umap.cca", features="EGFP")
#Create violin plots investigating common microglia markers
VlnPlot(all.integrated, features=c("ITGAM","P2RY12","CX3CR1","AIF1","CSF1R", "SALL1"), ncol = 3)
#Create violin plots for common activated microglia markers
VlnPlot(all.integrated, features=c("CD74","MAMU-DRA","F13A1","IRF8","CCR1", "ADGRE1"), ncol = 3)
#Create violin plots for common disease associated genes
VlnPlot(all.integrated, features=c("APOE","PZP","CD36","IL18","TGFB1", "CAPG"), ncol = 3)
#Heatmap
integrated.microglia.markers <- FindAllMarkers(all.integrated, only.pos = TRUE)
integrated.microglia.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
integrated.microglia.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(all.integrated, features = top10$gene) + NoLegend()
#Look at genes defining cluster ID
cluster2.markers <- FindMarkers(all.integrated, ident.1 = 2)
head(cluster2.markers, n = 40)
cluster0.markers <- FindMarkers(all.integrated, ident.1 = 0)
head(cluster0.markers, n = 40)
#find differentially expressed genes
DiffEX <- FindMarkers(all.integrated, ident.1=2, ident.2=0)
head(DiffEX,n=40)
#Create volcano plot with differentially expressed genes
DiffEX$diffexpressed<-"NO"
DiffEX$diffexpressed[DiffEX$avg_log2FC>1 & DiffEX$p_val<0.0000000000000000000000005] <- "UP"       
DiffEX$diffexpressed[DiffEX$avg_log2FC< -1 & DiffEX$p_val<0.0000000000000000000000005] <- "DOWN"
DiffEX$label <- "NA"
#Create column of gene IDs
DiffEX$names <- rownames(DiffEX)
DiffEX$label[DiffEX$diffexpressed !="NO"] <- DiffEX$names[DiffEX$diffexpressed != "NO"] 
#Attempt at labeling points 
ggplot(data=DiffEX, aes(x=avg_log2FC,y=-log10(p_val),col=diffexpressed, label=label))+ geom_point()+theme_minimal() + geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(0.0000000000000000000000005), col="red") + scale_color_manual(values=c("blue", "black", "red")) + geom_text()
#Volcano plot sans labels
ggplot(data=DiffEX, aes(x=avg_log2FC,y=-log10(p_val),col=diffexpressed, label=label))+ geom_point()+theme_minimal() + geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(0.0000000000000000000000005), col="red") + scale_color_manual(values=c("blue", "black", "red"))
