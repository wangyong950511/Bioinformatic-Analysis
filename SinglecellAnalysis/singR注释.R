library(Seurat)
# 进行singleR注释
pbmc.data <- Read10X(data.dir = "/Users/wangyong/Documents/CASH/CASHscRNA/CASH_N3")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3,
                           min.features = 200)
pbmc
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 2700
## Number of edges: 97892
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8719
## Number of communities: 9
## Elapsed time: 0 seconds
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
pbmc
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 2000 variable features)
##  3 dimensional reductions calculated: pca, umap, tsne


plot1 <- DimPlot(pbmc, reduction = "pca", label = TRUE)
plot2 <- DimPlot(pbmc, reduction = "tsne", label = TRUE)
plot3 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
plot1 + plot2 + plot3



library(SingleR)
library(celldex)
### 数据库加载 
hpca.se <- celldex::HumanPrimaryCellAtlasData()
bpe.se<-celldex::BlueprintEncodeData()
save(hpca.se,file='HumanPrimaryCellAtlas_hpca.se_human.RData')
save(bpe.se,file='BlueprintEncode_bpe.se_human.RData')
load("HumanPrimaryCellAtlas_hpca.se_human.RData")
load("BlueprintEncode_bpe.se_human.RData")

##获取标准化矩阵
pbmc4SingleR <- GetAssayData(pbmc, slot = "data")  
pbmc.hesc <- SingleR(test = pbmc4SingleR, ref = hpca.se, labels = hpca.se$label.main)  #
pbmc.hesc

# seurat 和 singleR的table表
table(pbmc.hesc$labels, pbmc$seurat_clusters)
pbmc@meta.data$labels <- pbmc.hesc$labels
DimPlot(pbmc, group.by = c("seurat_clusters", "labels"), reduction = "pca")

DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "tsne")

DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "umap")



# 多个数据库注释
pbmc.hesc <- SingleR(test = pbmc4SingleR, ref = list(BP = bpe.se, HPCA = hpca.se),
                     labels = list(bpe.se$label.main, hpca.se$label.main))
table(pbmc.hesc$labels, pbmc$seurat_clusters)
pbmc@meta.data$labels <- pbmc.hesc$labels
DimPlot(pbmc, group.by = c("seurat_clusters", "labels"), reduction = "pca")

DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "tsne")

DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "umap")

### Annotation diagnostics
plotScoreHeatmap(pbmc.hesc)
plotDeltaDistribution(pbmc.hesc, ncol = 3)

tab <- table(label = pbmc.hesc$labels, cluster = pbmc$seurat_clusters)
tab
library(pheatmap)
pheatmap(log10(tab + 10))