# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(ggplot2)
library(pathview)



####GO富集分析#####
# 
ego <- enrichGO(
  gene = symbols,                 # 基因列表
  OrgDb = org.Hs.eg.db,           # 物种数据库（人类为org.Hs.eg.db）
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
gene_entrez <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(gene_entrez)
# 执行KEGG富集分析
ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,    # 使用ENTREZ ID
  organism = "hsa",               # 小鼠的物种代码为 "mmu"
  pAdjustMethod = "BH",           # p值调整方法
  pvalueCutoff = 1,            # p值阈值
  qvalueCutoff = 0.9              # q值阈值
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
  showCategory = 10,  # 显示最多20个通路
  title = "KEGG Signaling Pathways Enrichment",
  color = "qvalue",  # 根据调整后的 Q 值着色
  size = "Count"      # 根据基因计数调整气泡大小
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


# # 指定KEGG通路图
# pathview(
#   gene.data = gene_entrez$ENTREZID,
#   pathway.id = "hsa04024",        # 替换为感兴趣的通路ID
#   species = "hsa",
#   limit = list(gene = 2, cpd = 1)
# )
# 
# 
gene_list <-df3
# 设置 log2FC 为基因表达值
gene_list_unique <- gene_list %>%
  dplyr::group_by(ENTREZID.x) %>%
  dplyr::summarise(logFC = max(logFC.x)) %>%
  dplyr::ungroup()
# 根据实际情况调整
gene_list_unique$logFC <- -gene_list_unique$logFC
geneList <- gene_list_unique$logFC
names(geneList) <- gene_list_unique$ENTREZID.x

# 运行 pathview，生成带颜色的 KEGG 通路图
pathview(gene.data = geneList,           # 基因表达值
         pathway.id = "hsa04210",            # cAMP 信号通路的 KEGG ID
         species = "hsa",                 # 物种（hsa 代表人类）
         limit = list(gene = c(-3, 3)),   # 表达值的上下限
         low = "blue", high = "red", mid = "white",) # 颜色设置：上调为红，下调为蓝







####GSEA富集分析#####
library(msigdbr)
gene_list <-df3
# 假设 gene_list 数据框包含 SYMBOL 和 log2FC 列
# 对基因名去重，保留 log2FC 值最大的行
gene_list_unique2 <- gene_list %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(logFC = max(logFC.x)) %>%
  dplyr::ungroup()
gene_list_unique2$logFC <- -gene_list_unique2$logFC
# 构建 geneList：一个命名的数值向量
geneList <- gene_list_unique2$logFC
names(geneList) <- gene_list_unique2$SYMBOL
geneList <- geneList[!is.na(names(geneList))]  # 去掉 NA 的基因
# 排序 geneList，确保从高到低排列
geneList <- sort(geneList, decreasing = TRUE)
# 检查结果
head(geneList)


# 加载 KEGG
kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_gene_sets <- kegg_gene_sets[, c("gs_name", "gene_symbol")]
# 加载 Reactome 基因集
reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_gene_sets <- reactome_gene_sets[, c("gs_name", "gene_symbol")]
# 加载 GO:BP 基因集
go_bp_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_bp_gene_sets <- go_bp_gene_sets[, c("gs_name", "gene_symbol")]
# 合并 KEGG、Reactome 和 GO:BP 基因集
combined_gene_sets <- rbind(kegg_gene_sets, reactome_gene_sets, go_bp_gene_sets)
# 提取合并基因集中的基因符号
combined_genes <- unique(combined_gene_sets$gene_symbol)
# 检查 geneList 与合并基因集的重叠
combined_overlap_genes <- intersect(names(geneList), combined_genes)
# 输出重叠基因数量
cat("重叠基因数（KEGG + Reactome + GO:BP）：", length(combined_overlap_genes), "\n")
gsea_result_combined <- GSEA(geneList = geneList, 
                             TERM2GENE = combined_gene_sets, 
                             minGSSize = 5,  # 最小基因集大小
                             maxGSSize = 1000,  # 最大基因集大小
                             pvalueCutoff = 0.5)  # p值阈值

# 查看 GSEA 结果
head(gsea_result_combined@result)
# 绘制点图
dotplot(gsea_result_combined, 
        showCategory = 10, 
        color = "pvalue",
        x = "NES",
        title = "Top Enriched Pathways")+  # 设置标题
  xlim(c(-3, 3))
# 绘制条形图
top_pathways <- gsea_result_combined[1:10, ]  # 选择前 10 个通路
ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Enriched Pathways",
       x = "Pathways",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal()
# 绘制显著通路的富集曲线
gseaplot2(
  gsea_result_combined, 
  geneSetID = gsea_result_combined$ID[2],  # 替换为实际通路 ID
  title = gsea_result_combined$Description[2]
)


# 提取某个显著通路的核心基因信息
pathway_id <- "GOBP_DEVELOPMENTAL_CELL_GROWTH"  # 替换为显著通路 ID
core_genes <- gsea_result_combined@result[pathway_id, "core_enrichment"]
# 分割核心基因字符串为单独基因
core_genes_list <- strsplit(core_genes, "/")[[1]]
# 查看核心基因列表
print(core_genes_list)
