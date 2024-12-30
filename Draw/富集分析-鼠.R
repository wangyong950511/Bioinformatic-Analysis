# 加载必要的包
library(clusterProfiler)
library(org.Mm.eg.db)  # 对于人类基因
library(enrichplot)
library(ggplot2)
library(pathview)



####GO富集分析#####
# 
ego <- enrichGO(
  gene = symbols,                 # 基因列表
  OrgDb = org.Mm.eg.db,           # 物种数据库（小鼠为org.Mm.eg.db）
  keyType = "SYMBOL",             # 基因标识符类型
  ont = "ALL",                    # 富集的GO分类（ALL、BP、MF、CC）
  pAdjustMethod = "BH",           # p值调整方法（Benjamini-Hochberg）
  pvalueCutoff = 0.05,            # p值阈值
  qvalueCutoff = 0.2              # q值阈值
)
# 查看富集结果
head(ego)
# 可视化结果：柱状图
barplot(ego, showCategory = 10, title = "GO Enrichment Analysis")
# 可视化结果：气泡图
dotplot(ego, showCategory = 10, title = "GO Enrichment Analysis")



####KEGG富集分析#####
# 转换基因符号为ENTREZ ID
gene_entrez <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
head(gene_entrez)
# 执行KEGG富集分析
ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,    # 使用ENTREZ ID
  organism = "mmu",               # 人类物种代码为"hsa"，小鼠的物种代码为 "mmu"
  pAdjustMethod = "BH",           # p值调整方法
  pvalueCutoff = 1,            # p值阈值
  qvalueCutoff = 0.2              # q值阈值
)
# 提取 KEGG 结果为数据框
kegg_results <- as.data.frame(ekegg)
# 查看富集结果
head(ekegg)
# 可视化结果：柱状图
barplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
# 可视化结果：气泡图
dotplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")



# 筛选信号通路相关结果
signal_pathways <- ekegg@result[grep("signaling", ekegg@result$Description, ignore.case = TRUE), ]
# 将筛选的信号通路重新创建为 enrichResult 对象
signal_ekegg <- ekegg
signal_ekegg@result <- signal_pathways
# 绘制气泡图
dotplot(
  signal_ekegg, 
  showCategory = 20,  # 显示最多20个通路
  title = "KEGG Signaling Pathways Enrichment",
  color = "qvalue",  # 根据调整后的 Q 值着色
  size = "Count"       # 根据基因计数调整气泡大小
)

# #自定义气泡图
# dotplot(
#   signal_ekegg,
#   showCategory = 20,
#   title = "KEGG Signaling Pathways Enrichment",
#   color = "p.adjust",
#   size = "Count"
# ) +
#   scale_color_gradientn(
#     colors = c("blue", "yellow")  # 设置蓝色-黄色-红色的三色渐变
#   )  + # 设置蓝色-黄色-红色的三色渐变
#   theme(
#     panel.background = element_rect(fill = "white"), # 设置面板背景色
#     panel.grid.major = element_line(color = "lightgray"),    # 设置主网格线颜色
#     panel.grid.minor = element_line(color = "lightgray"), # 设置次网格线颜色
#     plot.title = element_text(size = 16, face = "bold"),  # 设置标题字体
#     axis.title = element_text(size = 14),                 # 设置坐标轴标题字体大小
#     axis.text = element_text(size = 12),                  # 设置坐标轴刻度字体大小
#     legend.position = "bottom"                            # 设置图例位置
#   )


# 指定KEGG通路图
pathview(
  gene.data = gene_entrez$ENTREZID,
  pathway.id = "mmu04657",        # 替换为感兴趣的通路ID
  species = "mmu",
  limit = list(gene = 2, cpd = 1)
)



####GSEA富集分析####
Tcells_degs$gene <- rownames(Tcells_degs)
ids=bitr(Tcells_degs$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
Tcells_degs=merge(Tcells_degs,ids,by.x='gene',by.y='SYMBOL')
head(Tcells_degs)
Tcells_degs <- Tcells_degs[order(Tcells_degs$avg_log2FC,decreasing = T),]
cluster3.markers_list <- as.numeric(Tcells_degs$avg_log2FC)
names(cluster3.markers_list) <- Tcells_degs$ENTREZID
head(cluster3.markers_list)

cluster3_gsekg <- gseKEGG(cluster3.markers_list,organism = "mmu",pvalueCutoff = 0.05)
head(cluster3_gsekg)
cluster3_gsekg_arrange <- arrange(cluster3_gsekg,desc(abs(NES)))
head(cluster3_gsekg_arrange)

color <- c("#f7ca64","#43a5bf","#86c697","#a670d6","#ef998a")
gsekp1 <- gseaplot2(cluster3_gsekg_arrange,1:5,color = color,pvalue_table = F,base_size = 14)
gsekp1
gsekp2 <- upsetplot(cluster3_gsekg_arrange,n=5)
gsekp2 
