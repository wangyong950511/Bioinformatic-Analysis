library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(SingleR)
library(SeuratData)



####多样本同种数据读取方式#####
dir_name1 = list.files("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/NGB/")
scRNAlist <- list()
for (i in 1:length(dir_name1)) {
  counts <- Read10X(data.dir = paste("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/NGB/", dir_name1[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(
    counts,
    project = dir_name1[i],
    min.cells = 3,
    min.features = 300
  )
}

dir_name2 = list.files("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IMP/")
for (i in 1:length(dir_name2)) {
  counts <- Read10X(data.dir = paste("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IMP/", dir_name2[i], sep = ""))
  scRNAlist[[i+3]] <- CreateSeuratObject(
    counts,
    project = dir_name2[i],
    min.cells = 3,
    min.features = 300
  )
}

dir_name3 = list.files("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IWP/")
for (i in 1:length(dir_name3)) {
  counts <- Read10X(data.dir = paste("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IWP/", dir_name3[i], sep = ""))
  scRNAlist[[i+7]] <- CreateSeuratObject(
    counts,
    project = dir_name3[i],
    min.cells = 3,
    min.features = 300
  )
}

dir_name4 = list.files("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IMR/")
for (i in 1:length(dir_name4)) {
  counts <- Read10X(data.dir = paste("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IMR/", dir_name4[i], sep = ""))
  scRNAlist[[i+11]] <- CreateSeuratObject(
    counts,
    project = dir_name4[i],
    min.cells = 3,
    min.features = 300
  )
}

dir_name5 = list.files("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IWR/")
for (i in 1:length(dir_name5)) {
  counts <- Read10X(data.dir = paste("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/IWR/", dir_name5[i], sep = ""))
  scRNAlist[[i+17]] <- CreateSeuratObject(
    counts,
    project = dir_name5[i],
    min.cells = 3,
    min.features = 300
  )
}


dir_name <- c(dir_name1, dir_name2, dir_name3, dir_name4, dir_name5)
scRNAlist[]
dir_name[]




####质控#####
for (i in 1:length(scRNAlist)) {
  sc <- scRNAlist[[i]]
  # 鼠为mt
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  # 计算红细胞比例
  HB_genes <- c("HBA1",
                "HBA2",
                "HBB",
                "HBD",
                "HBE1",
                "HBG1",
                "HBG2",
                "HBM",
                "HBQ1",
                "HBZ")
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))
  HB_genes <- rownames(sc@assays$RNA)[HB_m]
  HB_genes <- HB_genes[!is.na(HB_genes)]
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features = HB_genes)
  # # 小鼠红细胞比列
  # sc[["HB_percent"]] <- PercentageFeatureSet(sc, pattern = "^Hb[^(p)]")
  # 赋值
  scRNAlist[[i]] <- sc
  rm(sc)
}


####小提琴图#####
violin_befor <- list()
for (i in 1:length(scRNAlist)) {
  violin_befor[[i]] <- VlnPlot(
    scRNAlist[[i]],
    layer = "counts",
    features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"),
    pt.size = 0.01,
    ncol = 3
  )
}
violin_befor
violin_befor[[2]]



####质控#####
# nFeature_RNA：每个细胞检测表达的基因数大于200，小于7000；
# nCount_RNA：每个细胞测序的UMI count含量大于1000，且剔除最大的前3的细胞；
# mt_RNA每个：细胞的线粒体基因表达量占总体基因的比例小于10%；
# HB_percent：每个细胞的红细胞基因表达量占总体基因的比例小于3%；
scRNAlist <- lapply(
  X = scRNAlist,
  FUN = function(x) {
    x <- subset(
      x,
      subset = nFeature_RNA > 200 &
        nFeature_RNA < 3000 &
        mt_percent < 15 &
        HB_percent < 3 &
        nCount_RNA < quantile(nCount_RNA, 0.97) &
        nCount_RNA > 1000
    )
  }
)


####合并#####
# 合并
scRNAlist_merge <- merge(x = scRNAlist[[1]],
                         y = scRNAlist[-1],
                         add.cell.ids = dir_name)
VlnPlot(
  scRNAlist_merge,
  features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"),
  split.by = "orig.ident",
  pt.size = 0.01,
  ncol = 4,
  layer = "counts",
)
table(scRNAlist_merge[[]]$orig.ident)

# 添加新的分组信息到 meta.data 中
scRNAlist_merge$group <- ifelse(grepl("NGB", scRNAlist_merge$orig.ident), "normal", "tumor")
# 检查新分组情况
table(scRNAlist_merge$group)

####标准化&高变基因&归一化&PCA#####
scRNAlist_merge <- NormalizeData(scRNAlist_merge)
scRNAlist_merge <- FindVariableFeatures(scRNAlist_merge)
scRNAlist_merge <- ScaleData(scRNAlist_merge, vars.to.regress = c("mt_percent"))
scRNAlist_merge <- RunPCA(scRNAlist_merge, verbose = T)



####整合去批次#####
# harmony
scRNA_harmony <- IntegrateLayers(
  object = scRNAlist_merge, method = HarmonyIntegration, orig.reduction = "pca",
  new.reduction = 'harmony', verbose = TRUE
)


####合并样本#####
scRNA_harmony[["RNA"]] <- JoinLayers(scRNA_harmony[["RNA"]])



####降维聚类#####
ElbowPlot(scRNA_harmony, ndims = 50)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims =
                                 1:20)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.1)
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:20, reduction = "harmony")
scRNA_harmony <- RunTSNE(scRNA_harmony, dims = 1:20, reduction = "harmony")
# 保存降维后的数据
save(scRNA_harmony, file = "scRNA_harmony.Rdata")
# 聚类分辨率选择
clustree(scRNA_harmony)
# 整合情况分析
DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Harmony")
Idents(scRNA_harmony) <- "RNA_snn_res.0.1"
DimPlot(scRNA_harmony, reduction = "umap")


####绘图#####
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated1
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group")
umap_integrated2
tsne_integrated3 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "seurat_clusters",label = FALSE)+
  ggtitle("")
tsne_integrated3
umap_integrated4 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
umap_integrated4
umap_tsne_integrated <- umap_integrated1 + umap_integrated2 + tsne_integrated3 + umap_integrated4 +
  plot_layout(ncol = 2)
umap_tsne_integrated




plots <- VlnPlot(
  scRNA_harmony,
  features = c("CD19", "CD3E","CD3G"),
  group.by = "celltype",
  pt.size = 0,
  combine = FALSE
)
# 修改x轴文本方向为水平
plots <- lapply(plots, function(p) {
  p + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
})
wrap_plots(plots=plots,ncol = 1)




celltype = data.frame(ClusterID = 0:7, celltype = 'unkown')
celltype[celltype$ClusterID %in% c(0,6), 2] = 'Microglial Cells'
celltype[celltype$ClusterID %in% c(2), 2] = 'Conventional dendritic cells'
celltype[celltype$ClusterID %in% c(3), 2] = 'CD8 T cell'
celltype[celltype$ClusterID %in% c(1), 2] = 'T cell'
celltype[celltype$ClusterID %in% c(5), 2] = 'Cancer cell'
celltype[celltype$ClusterID %in% c(4), 2] = 'B cell'
celltype[celltype$ClusterID %in% c(7), 2] = 'Neural progenitor cells'
celltype
table(celltype$celltype)




scRNA_harmony@meta.data$celltype="NA"
Idents(scRNA_harmony) <- "RNA_snn_res.0.1"
for (i in 1:nrow(celltype)) {
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$RNA_snn_res.0.1==celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
table(scRNA_harmony@meta.data$celltype)



p.cell <- DimPlot(scRNA_harmony, reduction = "tsne", group.by="celltype",label = TRUE)+
  ggtitle("")
p.cell
p.cell <- DimPlot(scRNA_harmony, reduction = "umap", group.by="celltype",label = TRUE)+
  ggtitle("")
p.cell


VlnPlot(scRNA_harmony, features = c("TRIP6"), group.by = "celltype",split.by="group")
# 分组小提琴图
plots <- VlnPlot(
  scRNA_harmony,
  features = c("GZMB"),
  group.by = "celltype",
  split.by="group",
  pt.size = 0,
  combine = FALSE
)
wrap_plots(plots=plots,ncol = 1)


# 保存聚类后的数据
save(scRNA_harmony, file = "scRNA_harmony_celltype.Rdata")

load("~/R/scRNA/Data/Glioma/Glioma-GSE222520-FACS/scRNA_harmony.Rdata")


