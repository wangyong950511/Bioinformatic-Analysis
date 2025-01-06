library(Seurat)
library(multtest)
library(dplyr)
library(tidyverse)
library(harmony)
library(clustree)
library(patchwork)
library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)




# 设置工作目录
setwd("/home/drwang/R/scRNA/workdata")



####多样本读取#####
# 读取CASH
parent_dir <- "~/R/scRNA/Data/新辅助化疗后结直肠癌及癌旁/CASH/"
dir_name <- list.files(parent_dir)
dir_name
scRNAlist <- list()
for (i in 1:length(dir_name)) {
  counts <- Read10X(data.dir = paste("~/R/scRNA/Data/新辅助化疗后结直肠癌及癌旁/CASH/", dir_name[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(
    counts,
    project = dir_name[i],
    min.cells = 3,
    min.features = 300,
  )
}
# 读取NORMAL
seurat_data <- read.csv(gzfile("~/R/scRNA/Data/Normal2/GSE115469_Data.csv.gz"),
                        row.names = 1)
scRNAlist[[5]] <- CreateSeuratObject(
  counts = seurat_data,
  min.features = 200,
  min.cells = 3,
  project = "GSE115469"
)
new_name <- "GSE115469"
dir_name <- c(dir_name, new_name)
####多样本同种数据读取方式#####
dir_name2 = list.files("~/R/scRNA/Data/Normal1/")
for (i in 1:length(dir_name2)) {
  file_path <- paste0("~/R/scRNA/Data/Normal1/", dir_name2[i])  
  counts <- read.csv(gzfile(file_path), row.names = 1)
  scRNAlist[[i+5]] <- CreateSeuratObject(
    counts,
    project = sub("_.*$", "", dir_name2[i]),
    min.cells = 3,
    min.features = 300
  )
}
dir_name <- c(dir_name, dir_name2)
scRNAlist[]


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
# scRNAlist_merge@assays$RNA@layers$counts.CASH_N3
# 分组
scRNAlist_merge$group <- ifelse(
  scRNAlist_merge$orig.ident %in% c("CASH_N3", "CASH_N4", "CASH_N5", "CASH_N6", "CASH_N7"),
  "Exp",
  "Ctrl"
)
table(scRNAlist_merge[[]]$group)

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
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = seq(from = 0.1, to =
                                                                1.0, by = 0.1))
# scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.1)
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:20, reduction = "harmony")
scRNA_harmony <- RunTSNE(scRNA_harmony, dims = 1:20, reduction = "harmony")
# 保存降维后的数据
save(scRNA_harmony, file = "scRNA_harmony.Rdata")
# 聚类分辨率选择
clustree(scRNA_harmony)
# 整合情况分析
DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Harmony")
# 更改分辨率
# scRNA_harmony$RNA_snn_res.0.1
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
# 存图
ggsave("umap_tsne_integrated.pdf",
       umap_tsne_integrated,
       wi = 25,
       he = 15)

DimPlot(scRNA_harmony, reduction = "umap", split.by = "group")


####手动注释#####
# 计算所有marker
all.markers <- FindAllMarkers(scRNA_harmony,, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
# 查看标记基因
top10_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% 
                                top_n(n = 10, wt = avg_log2FC))

# 内皮细胞 "PECAM1", "VWF","CDH5","KDR","NOS3"
# 肝细胞   "ALB", "CYP3A4","HNF4A","APOE"
# 巨噬细胞 "CD68", "CD163","CSF1R","MRC1","FCGR3A","MARCO"
# T细胞 "CD3E", "CD4","CD8A","FOXP3","CD25","CD28","CD45RO","CCR7"
# NK细胞 "NCAM1", "KLRD1","KLRC1","KLRK1","FCGR3A","NCR1","PRF1"
# B细胞 "CD19","CD38","CD20","CD79A","CD79B","CD19","CD24"
# 胆管   "KRT19", "KRT7","EPCAM","SOX9"
# HSCs    "ACTA2","GFAP","Col1A1 ","MMP-2" 
#    "KRT19", "KRT7","EPCAM","SOX9"





plots <- VlnPlot(
  scRNA_harmony,
  features = c("CD3E", "CD4","CD8A","FOXP3","CD25","CD28","CD45RO","CCR7"),
  group.by = "seurat_clusters",
  pt.size = 0,
  combine = FALSE
)
wrap_plots(plots=plots,ncol = 1)




celltype = data.frame(ClusterID = 0:12, celltype = 'unkown')
celltype[celltype$ClusterID %in% c(6,5), 2] = 'Macrophages'
celltype[celltype$ClusterID %in% c(11), 2] = 'Cholangiocyte'
celltype[celltype$ClusterID %in% c(7), 2] = 'B cell'
celltype[celltype$ClusterID %in% c(0,1), 2] = 'T cell'
celltype[celltype$ClusterID %in% c(3), 2] = 'Endothelial cell'
celltype[celltype$ClusterID %in% c(4), 2] = 'Live cell'
celltype[celltype$ClusterID %in% c(2), 2] = 'NK cell'
celltype[celltype$ClusterID %in% c(8), 2] = 'Monocyte'
celltype[celltype$ClusterID %in% c(9), 2] = 'Plasma cell'
celltype[celltype$ClusterID %in% c(10,12), 2] = 'HSCs'
celltype
table(celltype$celltype)


scRNA_harmony@meta.data$celltype="NA"
Idents(scRNA_harmony) <- "seurat_clusters"
for (i in 1:nrow(celltype)) {
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters==celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
table(scRNA_harmony@meta.data$celltype)



p.cell <- DimPlot(scRNA_harmony, reduction = "tsne", group.by="celltype",label = TRUE)+
  ggtitle("")
p.cell
p.cell <- DimPlot(scRNA_harmony, reduction = "umap", group.by="celltype",label = TRUE)+
  ggtitle("")
p.cell



save(scRNA_harmony, file = "scRNA_harmony2.Rdata")


####热图#####
top5 = all.markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
top5 <- as.data.frame(top5)
markerdata <- ScaleData(scRNA_harmony,features = as.character(unique(top5$gene)),assay = "RNA")
DoHeatmap(markerdata,features = as.character(unique(top5$gene)),
          group.by = "RNA_snn_res.0.1",
          assay = 'RNA')



####FeaturePlot图#####
FeaturePlot(scRNA_harmony, features = c("ALB", "CYP3A4","HNF4A","APOE"))
FeaturePlot(scRNA_harmony, features = c("CDKN1A", "MDM2","BAX","BBC3"))
FeaturePlot(scRNA_harmony, features = c("SREBF1", "FASN","PPARA","CPT1A"))


####VlnPlot图#####
VlnPlot(scRNA_harmony, features = c("ALB", "CYP3A4","HNF4A","APOE"),pt.size = 0)
VlnPlot(scRNA_harmony, features = c("CDKN1A", "MDM2","BAX","BBC3"), group.by = "group")
VlnPlot(scRNA_harmony, features = c("ALB"), group.by = "celltype")


####差异分析#####
Idents(scRNA_harmony)="celltype"
# 寻找标记基因
cell.markers <- FindAllMarkers(
  object = scRNA_harmony,
  test.use = "wilcox",
  only.pos = FALSE,
  logfc.threshold = 0.25,
  min.pct = 0.5
)
table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$seurat_clusters)

Tcells_degs=FindMarkers(scRNA_harmony,
                        test.use = "wilcox",
                        only.pos = FALSE,
                        logfc.threshold = 0.25,
                        min.pct = 0.5,ident.1="肝细胞",ident.2 = "NK细胞")
mutate(gene=rownames(.))




####细胞比例图#####
Idents(scRNA_harmony)="celltype"
table(scRNA_harmony$orig.ident)
prop.table(table(Idents(scRNA_harmony)))
table(Idents(scRNA_harmony),scRNA_harmony$group)
Cellratio <- prop.table(table(Idents(scRNA_harmony),scRNA_harmony$orig.ident),margin=2)
Cellratio <-as.data.frame(Cellratio)                        
colnames(Cellratio) <- c("celltype","sample","ratio")
# 反转
colourCount=length(unique(Cellratio$celltype))
ggplot(Cellratio)+
  geom_bar(aes(x=sample,y=ratio,fill = celltype),stat = "identity",width = 0.7,size=0.5,colour='#222222')+
  theme_classic()+
  labs(x='sample',y='Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill = NA,color = "black",size=0.5,linetype = "solid"))
# 正常
ggplot(Cellratio)+
  geom_bar(aes(x=sample,y=ratio,fill = celltype),stat = "identity",width = 0.7,size=0.5,colour='#222222')+
  theme_classic()+
  labs(x='sample',y='Ratio')+
  theme(panel.border = element_rect(fill = NA,color = "black",size=0.5,linetype = "solid"))



Idents(scRNA_harmony)="celltype"
Liver_cells <- subset(scRNA_harmony, idents = "Live cell")
Idents(Liver_cells)="group"
Liver_cells_cash <- subset(Liver_cells, subset = group == "Exp")
Liver_cells_ctrl <- subset(Liver_cells, subset = group == "Ctrl")
####两个cluster分析####
Tcells_degs = FindMarkers(
  Liver_cells,
  only.pos = FALSE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  ident.1 = "Exp", ident.2 = "Ctrl"
) %>% mutate(gene = rownames(.)) 

Tcells_degs_fil = Tcells_degs %>% filter(pct.1 > 0.1 & pct.2 > 0.1 &
                                           p_val_adj < 0.05) %>% filter(abs(avg_log2FC) > 0.5)
colnames(Tcells_degs_fil)
table(abs(Tcells_degs_fil$avg_log2FC)>2)
plotdt=Tcells_degs_fil %>% mutate(gene=ifelse(abs(avg_log2FC)>=2,gene,NA))
#火山图绘制
# plotdt$p_val_adj[plotdt$p_val_adj == 0] <- 1e-10
ggplot(plotdt,
       aes(
         x = avg_log2FC,
         y = -log10(p_val_adj),
         size = pct.1,
         color = avg_log2FC
       )) +
  geom_point() +
  ggtitle(label = "NK cell VS T cell", subtitle = "Antibody") +
  geom_text_repel(aes(label = gene), size = 3, color = "black") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent", colour = NA)
  ) +
  scale_color_gradient2(
    low = "olivedrab",
    high = "salmon2",
    mid = "grey",
    midpoint = 0
  ) +
  scale_size(range = c(1, 3))



####GO&KEGG富集分析####
# 数据整理
Tcells_degs_fil$gene <- rownames(Tcells_degs_fil)
ids=bitr(Tcells_degs_fil$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Tcells_degs_fil=merge(Tcells_degs_fil,ids,by.x='gene',by.y='SYMBOL')
Tcells_degs_fil <- Tcells_degs_fil[order(Tcells_degs_fil$avg_log2FC,decreasing = T),]
Tcells_degs_list <- as.numeric(Tcells_degs_fil$avg_log2FC)
names(Tcells_degs_list) <- Tcells_degs_fil$ENTREZID
cluster3_de <- names(Tcells_degs_list)[abs(Tcells_degs_list)>1]
# GO
cluster3_ego <- enrichGO(cluster3_de,OrgDb = "org.Hs.eg.db",ont="BP",readable = TRUE)
head(cluster3_ego)
dotplot(cluster3_ego,showCategory=10,title="CASH VS Normal GO")
# KEGG
cluster3_ekg <- enrichKEGG(cluster3_de,organism = "hsa",pvalueCutoff = 0.05)
head(cluster3_ekg)
dotplot(cluster3_ekg,showCategory=10,title="CASH VS Normal KEGG")





####GSEA富集分析####
Tcells_degs$gene <- rownames(Tcells_degs)
ids=bitr(Tcells_degs$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Tcells_degs=merge(Tcells_degs,ids,by.x='gene',by.y='SYMBOL')
head(Tcells_degs)
Tcells_degs <- Tcells_degs[order(Tcells_degs$avg_log2FC,decreasing = T),]
cluster3.markers_list <- as.numeric(Tcells_degs$avg_log2FC)
names(cluster3.markers_list) <- Tcells_degs$ENTREZID
head(cluster3.markers_list)

cluster3_gsekg <- gseKEGG(cluster3.markers_list,organism = "hsa",pvalueCutoff = 0.05)
head(cluster3_gsekg)
cluster3_gsekg_arrange <- arrange(cluster3_gsekg,desc(abs(NES)))
head(cluster3_gsekg_arrange)

color <- c("#f7ca64","#43a5bf","#86c697","#a670d6","#ef998a")
gsekp1 <- gseaplot2(cluster3_gsekg_arrange,1:5,color = color,pvalue_table = F,base_size = 14)
gsekp1
gsekp2 <- upsetplot(cluster3_gsekg_arrange,n=5)
gsekp2 
