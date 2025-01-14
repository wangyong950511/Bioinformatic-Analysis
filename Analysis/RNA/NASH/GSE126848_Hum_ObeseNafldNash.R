# -------------------------------
# Script Name: NC_OBESE_NAFLD_NASH_GSE126848.R
# Data：GSE126848
# Purpose: 分析NC-Obese-NAFLD-NASH
# Author: WangYong
# Date: 2025-1-10



# 加载库和设置 ------------------------------------
library(DESeq2)
library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(pheatmap)
library(org.Hs.eg.db)


####导入count矩阵和样本信息####
countData <- read.table("~/R/RNA/data/NASH/GSE126848/GSE126848_raw_counts_GRCh38.p13_NCBI.tsv.gz", header = TRUE, row.names = 1)
# 提取样本计数数据
count_matrix <- as.matrix(countData)  # 修剪数据


# 样本信息表
colData <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(c(
    rep("NAFLD", 15),
    rep("NASH", 16),
    rep("NC", 14),
    rep("Obese", 12)
  ), levels = c("NC", "NAFLD", "NASH", "Obese"))
)


# 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) > 10, ]

# 运行差异表达分析
dds <- DESeq(dds)

res_all <- results(dds)
summary(res_all)

res_sig <- subset(res_all, padj < 0.05 & abs(log2FoldChange) > 1)

# 数据标准化转换
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# 计算样本之间的距离矩阵
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# 绘制热图
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = colData,
         main = "Sample Distance Heatmap")
# 两两组差异表达分析
res_NAFLD_NC <- results(dds, contrast = c("condition", "NAFLD", "NC"))
res_sig_NAFLD_NC <- subset(res_NAFLD_NC, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(res_sig_NAFLD_NC, "DEGs_NAFLD_vs_NC.csv")

res_NASH_NC <- results(dds, contrast = c("condition", "NASH", "NC"))
res_sig_NASH_NC <- subset(res_NASH_NC, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(res_sig_NASH_NC, "DEGs_NASH_vs_NC.csv")

res_Obese_NC <- results(dds, contrast = c("condition", "Obese", "NC"))
res_sig_Obese_NC <- subset(res_Obese_NC, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(res_sig_Obese_NC, "DEGs_Obese_vs_NC.csv")



####实例分析####
# 单基因分析
analyze_gene_expression("PPARA", vsd, colData)
# 热图绘制
gene_set <- c(
  # 双链断裂修复 (DSB Repair)
  "ATM", "ATR", "BRCA1", "BRCA2", "CHEK1", "CHEK2", 
  "MRE11", "RAD50", "NBN", "RAD51", "RAD52", "XRCC2", "XRCC3",
  # 碱基切除修复 (BER)
  "PARP1", "APEX1", "OGG1", "MPG", "NTHL1", "UNG", "SMUG1",
  # 核苷酸切除修复 (NER)
  "XPC", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "DDB1", "DDB2",
  # 不匹配修复 (MMR)
  "MLH1", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2",
  # 同源重组修复 (HRR)
  "FANCD2", "FANCA", "FANCC", "FANCE", "FANCF", "FANCG", "FANCI",
  # 非同源末端连接 (NHEJ)
  "PRKDC", "LIG4", "XRCC4", "XRCC5", "XRCC6", "DCLRE1C", "NHEJ1",
  # 其他 DDR 关键基因
  "WRN", "BLM", "MUS81", "RAD18", "POLQ", "RPA1", "RPA2", "RPA3",
  "MARCO"
)
generate_ddr_heatmap(count_matrix, colData, gene_set, "Heatmap of DDR Genes")

####特定基因分析（预运行）####
# 定义函数
analyze_gene_expression <- function(GeneX, vsd, colData) {
  # Step 1: 查询基因的 NCBI Gene ID (ENTREZID)
  gene_info <- select(org.Hs.eg.db,
                      keys = GeneX,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
  
  # 检查查询结果
  if (nrow(gene_info) == 0 || is.na(gene_info$ENTREZID[1])) {
    stop("目标基因未找到对应的 ENTREZID，请检查基因名是否正确。")
  }
  
  # 提取 ENTREZID
  gene_id <- gene_info$ENTREZID[1]
  
  # Step 2: 检查基因是否存在于表达矩阵
  if (!(gene_id %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 3: 获取目标基因的表达值
  gene_expression <- assay(vsd)[gene_id, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NAFLD", "NASH"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 5: 绘制箱线图
  custom_colors_fill <- c("NC" = "white", "Obese" = "white", "NAFLD" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#2ca02c", "Obese" = "#FFC07C", "NAFLD" = "#ff7f0e", "NASH" = "#d62728")
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("NC", "Obese"), 
                                          c("NC", "NAFLD"), 
                                          c("NC", "NASH"), 
                                          c("Obese", "NASH"), 
                                          c("Obese", "NAFLD"), 
                                          c("NASH", "NAFLD")),
                       method = "t.test",
                       label = "p.signif") + # 两两比较，显示显著性标记
    labs(
      title = paste("Expression of", GeneX, "across Four Groups"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # 返回绘图
  return(plot)
}






####热图绘制（预运行）####
# 定义函数
generate_ddr_heatmap <- function(count_matrix, colData, gene_set, title = "Heatmap of DDR Genes") {
  # Step 1: 转换基因名为 ENTREZ ID
  gene_info <- select(
    org.Hs.eg.db, 
    keys = gene_set, 
    columns = c("ENTREZID", "SYMBOL"), 
    keytype = "SYMBOL"
  )
  
  # 去除没有匹配到的基因
  gene_info <- gene_info[!is.na(gene_info$ENTREZID), ]
  
  # 提取 ENTREZ ID
  gene_ids <- gene_info$ENTREZID
  
  # Step 2: 筛选基因 ID 存在于表达矩阵中的部分
  gene_ids_filtered <- gene_ids[gene_ids %in% rownames(count_matrix)]
  
  if (length(gene_ids_filtered) == 0) {
    stop("None of the genes in the provided gene set were found in the count matrix.")
  }
  
  # Step 3: 提取对应基因的表达数据
  heatmap_data <- count_matrix[gene_ids_filtered, ]
  
  # Step 4: Z-score 标准化
  heatmap_data_scaled <- t(scale(t(heatmap_data)))  # 按行标准化
  heatmap_data_scaled[is.na(heatmap_data_scaled)] <- 0  # 处理 NA 值
  
  # Step 5: 样本排序
  colData$condition <- factor(colData$condition, levels = unique(colData$condition))
  sorted_indices <- order(colData$condition)
  
  # 重新排序表达矩阵和分组信息
  heatmap_data_scaled <- heatmap_data_scaled[, sorted_indices]
  annotation_col <- data.frame(Condition = colData$condition[sorted_indices])
  rownames(annotation_col) <- colnames(heatmap_data_scaled)
  
  # Step 6: 替换表达矩阵的行名为 SYMBOL
  id_to_symbol <- setNames(gene_info$SYMBOL, gene_info$ENTREZID)
  rownames(heatmap_data_scaled) <- id_to_symbol[rownames(heatmap_data_scaled)]
  
  # 检查替换后的行名
  if (any(is.na(rownames(heatmap_data_scaled)))) {
    stop("Some gene SYMBOLs could not be mapped. Please check your input data.")
  }
  
  # Step 7: 自定义颜色
  annotation_colors <- list(
    Condition = c(
      "NC" = "#2ca02c",
      "Obese" = "#FFB347",
      "NAFLD" = "#ff7f0e",
      "NASH" = "#d62728"
    )
  )
  
  # Step 8: 绘制热图
  pheatmap(
    heatmap_data_scaled,
    cluster_rows = TRUE,  # 对基因进行聚类
    cluster_cols = FALSE,  # 不对样本进行聚类
    annotation_col = annotation_col,  # 样本分组注释
    annotation_colors = annotation_colors,  # 分组颜色映射
    fontsize_row = 10,  # 基因字体大小
    fontsize_col = 10,  # 样本字体大小
    show_rownames = TRUE,  # 显示基因名
    show_colnames = TRUE,  # 显示样本名
    main = title,  # 标题
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50)  # 颜色梯度
  )
}




