library(DESeq2)
library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)



####导入count矩阵和样本信息####
countData <- read.table("~/R/RNA/data/CASH/gene_counts.txt", header = TRUE, row.names = 1)
# 提取样本计数数据
count_matrix <- as.matrix(countData[, 6:ncol(countData)])  # 去除前五列

# 修改列名：去掉 "_sorted.bam"
colnames(count_matrix) <- gsub("_sorted.bam$", "", colnames(count_matrix))

# 创建 colData，确保包含 condition 列
colData <- data.frame(
  row.names = colnames(count_matrix),  # 使用 count_matrix 的列名
  condition = c("CASH", "CASH", "CASH", 
                "CASH_6W", "CASH_6W", "CASH_6W", 
                "NC", "NC", "NC", 
                "NCnet", "NCnet", "NCnet", "NCnet", "NCnet"),  # 根据样本数填充
  batch = c("Batch1", "Batch1", "Batch1", 
            "Batch1", "Batch1", "Batch1", 
            "Batch1", "Batch1", "Batch1", 
            "Batch2", "Batch2", "Batch2", "Batch2", "Batch2")  # 根据实际情况填充
)
# 筛选非 NCnet 的样本
colData_filtered <- colData[!(colData$condition %in% c("NCnet", "CASH_6W")), ]
count_matrix_filtered <- count_matrix[, rownames(colData_filtered)]
# 创建 DESeq2 数据集对象
dds_filtered <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered, 
  colData = colData_filtered, 
  design = ~ condition
)
# 运行差异表达分析
dds_filtered <- DESeq(dds_filtered)
# 获取比较结果，例：CASH_4W 和 NC 的差异表达分析
results_filtered <- results(dds_filtered, contrast = c("condition", "CASH", "NC"))
# 查看结果
head(results_filtered[order(results_filtered$pvalue), ])




####进行 PCA 分析####
rld_filtered <- vst(dds_filtered, blind = FALSE)
pcaData_filtered <- plotPCA(rld_filtered, intgroup = "condition", returnData = TRUE)
percentVar_filtered <- round(100 * attr(pcaData_filtered, "percentVar"))
# 绘制 PCA 图
pca_plot <- ggplot(pcaData_filtered, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_filtered[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_filtered[2], "% variance")) +
  theme_bw()
print(pca_plot)


####绘制热图####

# 指定感兴趣的基因集
gene_set <- c("Atm", "Atr", "Brca1", "Brca2", "Chk1", "Chk2",  
              "Rad51", "Mlh1", "Msh2", "Parp1", "Fancd2", 
              "Wrn","Xpc", "Rpa1", "Marco")  # 替换为实际基因名
# 检查基因是否存在于表达矩阵中
gene_set_filtered <- gene_set[gene_set %in% rownames(count_matrix_filtered)]
if (length(gene_set_filtered) == 0) {
  stop("None of the genes in the gene set are present in the count matrix.")
}
# 提取基因表达数据
heatmap_data <- count_matrix_filtered[gene_set_filtered, ]
# 标准化基因表达（例如 Z-score）
heatmap_data_scaled <- t(scale(t(heatmap_data)))  # 按行（基因）标准化
# 确保标准化后没有 NA 值（可能由于单基因或样本全零导致）
heatmap_data_scaled[is.na(heatmap_data_scaled)] <- 0
# 按分组排序样本（根据 colData_filtered 的 condition）
colData_filtered$condition <- factor(colData_filtered$condition, levels = unique(colData_filtered$condition))  # 确保条件顺序
sorted_indices <- order(colData_filtered$condition)  # 根据分组排序索引
# 重新排序表达矩阵和注释数据
heatmap_data_scaled <- heatmap_data_scaled[, sorted_indices]
annotation_col <- data.frame(Condition = colData_filtered$condition[sorted_indices])
rownames(annotation_col) <- colnames(heatmap_data_scaled)
# 检查数据一致性
if (!all(rownames(annotation_col) == colnames(heatmap_data_scaled))) {
  stop("Column names in heatmap data do not match annotation data.")
}
# 绘制热图
library(pheatmap)
pheatmap(
  heatmap_data_scaled,  # 标准化的表达矩阵
  cluster_rows = TRUE,  # 对基因聚类
  cluster_cols = FALSE,  # 不对样本聚类
  annotation_col = annotation_col,  # 样本分组注释
  fontsize_row = 10,    # 调整基因字体大小
  fontsize_col = 10,    # 调整样本字体大小
  main = "Heatmap of Genes Related to （DDR）"  # 图标题
)


####特定基因分析####
gene_of_interest <- "Slc27a6"  # 替换为您感兴趣的基因名称
if (gene_of_interest %in% rownames(results_filtered)) {
  gene_results <- results_filtered[gene_of_interest, ]
  log2_fold_change <- gene_results$log2FoldChange
  p_value <- gene_results$pvalue
  padj <- gene_results$padj
  # 打印结果
  cat("Log2 Fold Change:", log2_fold_change, "\n")
  cat("P-value:", p_value, "\n")
  cat("Adjusted P-value:", padj, "\n")
  # 获取该基因在不同样本中的表达数据
  gene_expression <- count_matrix_filtered[gene_of_interest, ]
  expression_data <- data.frame(
    sample = colnames(count_matrix_filtered),
    expression = gene_expression,
    condition = colData_filtered$condition
  )
  # 绘制箱线图
  boxplot <- ggplot(expression_data, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot() +
    labs(title = paste("Expression of", gene_of_interest), x = "Condition", y = "Expression") +
    theme_minimal()
  print(boxplot)
} else {
  # 基因不在结果中，输出提示信息
  cat("The gene", gene_of_interest, "is not found in the filtered results.\n")
}

####差异基因####
# 提取显著差异表达基因
# 移除 NA 值
filtered_results <- results_filtered[
  !is.na(results_filtered$padj) & !is.na(results_filtered$log2FoldChange), 
]
# 筛选显著差异基因：padj < 0.05 且 |log2FoldChange| > 1
sig_results <- filtered_results[
  filtered_results$padj < 0.05 & abs(filtered_results$log2FoldChange) > 1, 
]
# 转换为数据框便于查看和导出
sig_results_df <- as.data.frame(sig_results)
# 查看结果
head(sig_results_df)
# 导出为 CSV 文件（逗号分隔）
write.csv(sig_results_df, file = "differential_genes.csv", row.names = TRUE, quote = FALSE)



####绘制火山图####
Firehill<-results_filtered
# 将 DESeqResults 转换为 data.frame
Firehill <- as.data.frame(Firehill)
# 去除 log2FoldChange 或 padj 中包含 NA 的行
Firehill <- Firehill %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))
Firehill$gene <- rownames(Firehill)
# 校正转换
Firehill$logP <- -log10(Firehill$padj)
# 限制 log2FoldChange 范围
Firehill <- Firehill %>%
  mutate(log2FoldChange = ifelse(log2FoldChange > 5, 5, 
                                 ifelse(log2FoldChange < -5, -5, log2FoldChange)),
         logP = ifelse(logP > 30, 30, logP)  # 限制 logP 最大值为 40
         )

Firehill <- Firehill %>%
  mutate(Group = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "up",
    padj < 0.05 & log2FoldChange < -1 ~ "down",
    TRUE ~ "none"
  ))
table(Firehill$Group)
# 添加一个新列 Label
Firehill$Label <- ""
Firehill <- Firehill %>%
  arrange(padj)
up_genes <- head(Firehill$gene[which(Firehill$Group == "up")], 10)
down_genes <- head(Firehill$gene[which(Firehill$Group == "down")], 10)
top10_genes <- c(as.character(up_genes), as.character(down_genes))
Firehill$Label[match(top10_genes, Firehill$gene)] <- top10_genes

ggscatter(
  Firehill,
  x = "log2FoldChange",
  y = "logP",
  color = "Group",
  palette = c("#2f5688", "#BBBBBB", "#cc0000"),
  size = 1,
  label = Firehill$Label,
  font.label = 8,
  xlab="log2FC",
  ylab = "-log10(Adjust P-value)"
) + theme_base()+
  geom_hline(yintercept = 1.30,linetype="dashed")+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+
  geom_text_repel(
    data = subset(Firehill, Label != ""),
    aes(label = Label),
    size = 3,
    box.padding = 0.5,   # 标签框与点之间的间距
    point.padding = 0.5, # 点与标签之间的间距
    max.overlaps = 15,   # 最多允许多少个标签
    force = 1.5          # 调整标签之间的排斥力
  )






####富集分析####
# 加载必要的包
library(clusterProfiler)
library(org.Mm.eg.db)  # 对于人类基因"org.Hs.eg.db".
library(enrichplot)
library(ggplot2)
library(pathview)

symbols <- rownames(sig_results_df)
# 执行GO富集分析
ego <- enrichGO(
  gene = symbols,                 # 基因列表
  OrgDb = org.Mm.eg.db,           # 物种数据库（人类为org.Hs.eg.db）
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


# 转换基因符号为ENTREZ ID
gene_entrez <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
head(gene_entrez)
# 执行KEGG富集分析
ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,    # 使用ENTREZ ID
  organism = "mmu",               # 人类物种代码为"hsa"，小鼠为"mmu"
  pAdjustMethod = "BH",           # p值调整方法
  pvalueCutoff = 0.05,            # p值阈值
  qvalueCutoff = 0.2              # q值阈值
)
# 查看富集结果
head(ekegg)
# 可视化结果：柱状图
barplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
# 可视化结果：柱状图
svg("pathway_barplot.svg", width = 8, height = 6, family = "Times New Roman")  # 设置文件路径和图形大小
barplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
dev.off()  # 关闭图形设备，保存图形

# 可视化结果：气泡图
dotplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")


# 指定通路图
pathview(
  gene.data = gene_entrez$ENTREZID,
  pathway.id = "mmu00071",        # 替换为感兴趣的通路ID
  species = "mmu",
  limit = list(gene = 2, cpd = 1)
)
