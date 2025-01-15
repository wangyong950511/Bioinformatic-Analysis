# -------------------------------
# Script Name: NC_NAFLD_BorderlineNASH_DefiniteNASH_GSE83596.R
# Data：GSE83596
# Purpose: 分析STAM模型随时间的变化规律
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
library(org.Mm.eg.db)
library(limma)




####导入count矩阵和样本信息####
# 文件路径
file_path <- "~/R/RNA/data/NASH/GSE83596/GSE83596_series_matrix.txt.gz"
# 读取所有行
lines <- readLines(file_path)
# 筛选非数据行（以 "!" 开头）
non_data_lines <- grepl("^!", lines)  # 返回布尔向量，TRUE 表示是非数据行
# 打印非数据行的总数
cat("非数据行的总数: ", sum(non_data_lines), "\n")
# 保存非数据行到文件
writeLines(lines[non_data_lines], "non_data_lines.txt")
# 提取数据行（非以 "!" 开头的行）
data_lines <- lines[!non_data_lines]
# 将数据行写入临时文件
writeLines(data_lines, "filtered_data.txt")
# 从临时文件读取数据
series_data <- read.table("filtered_data.txt", 
                          header = TRUE, 
                          sep = "\t", 
                          stringsAsFactors = FALSE)
# 转换行名为基因名
rownames(series_data) <- series_data$ID_REF
series_data <- series_data[, -1]
# 标准化数据
series_data_normalized <- normalizeBetweenArrays(series_data, method = "quantile")
# 样本分类信息
colData <- data.frame(
  row.names = colnames(series_data),  # 使用 count_matrix 或 series_data 的列名
  condition = factor(c(
    rep("Control", 14),   # Control
    rep("STAM_6w", 3),      # STAM 6 weeks
    rep("STAM_8w", 3),      # STAM 8 weeks
    rep("STAM_12w", 4),     # STAM 12 weeks
    rep("STAM_20w_Tumor", 4),    # STAM 20 weeks; tumor tissue
    rep("STAM_20w_NonTumor", 4)  # STAM 20 weeks; non-tumor tissue
  ), levels = c("Control",
                "STAM_6w", "STAM_8w", "STAM_12w", "STAM_20w_Tumor", "STAM_20w_NonTumor"))
)
# 假设 colData 是样本的分组信息
design <- model.matrix(~ 0 + colData$condition)
colnames(design) <- levels(colData$condition)




####实例分析####
# 单基因分析
plot_gene_expression("Cdkn1a", series_data, colData)
# 定义基因集合
gene_set_mouse <- c("Cdkn1a", "Trp53", "Myc", "Srebf1", "Srebf2", "Fasn", "Ppara")
# 调用函数
generate_mouse_heatmap(gene_set_mouse, series_data, colData, "Heatmap of Mouse DDR Genes")

####特定基因分析（预运行）####
# 定义函数
plot_gene_expression <- function(GeneX, series_data, colData) {
  # Step 1: 查询基因的 NCBI Gene ID (ENTREZID)
  gene_info <- select(org.Mm.eg.db,
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
  if (!(gene_id %in% rownames(series_data))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 3: 提取基因表达数据
  gene_expression <- series_data[gene_id, ]
  gene_expression <- as.numeric(gene_expression[1, ])
  
  # Step 4: 检查 colData 的行名是否与样本 ID 一致
  if (!all(rownames(colData) == colnames(series_data))) {
    colData <- colData[colnames(series_data), ]
  }
  
  # Step 5: 构建数据框
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition  # 样本分组信息
  )
  
  # Step 6: 确保分组是因子类型并按顺序排列
  gene_df$Condition <- factor(gene_df$Condition, levels = c("Control",
                                                            "STAM_6w", "STAM_8w", "STAM_12w",
                                                            "STAM_20w_Tumor", "STAM_20w_NonTumor"))
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 7: 自定义颜色
  custom_colors_fill <- c("Control" = "#F1C40F", "STAM_6w" = "#9B59B6", 
                          "STAM_8w" = "#2980B9", "STAM_12w" = "#1ABC9C", 
                          "STAM_20w_Tumor" = "#C0392B", "STAM_20w_NonTumor" = "#27AE60")
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    labs(
      title = paste("Expression of", GeneX, "across Groups"),
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
generate_mouse_heatmap <- function(gene_set_mouse, series_data, colData, title = "Heatmap of Mouse Genes") {
  # Step 1: 转换基因名为 ENTREZ ID
  gene_info <- select(
    org.Mm.eg.db, 
    keys = gene_set_mouse, 
    columns = c("ENTREZID", "SYMBOL"), 
    keytype = "SYMBOL"
  )
  
  # 去除没有匹配到的基因
  gene_info <- gene_info[!is.na(gene_info$ENTREZID), ]
  
  # 提取 ENTREZID 和 SYMBOL
  gene_ids <- gene_info$ENTREZID
  gene_symbols <- gene_info$SYMBOL
  
  # Step 2: 检查基因是否存在于表达矩阵
  missing_genes <- gene_ids[!(gene_ids %in% rownames(series_data))]
  if (length(missing_genes) > 0) {
    warning("以下基因不存在于表达矩阵中：", paste(missing_genes, collapse = ", "))
  }
  gene_ids <- gene_ids[gene_ids %in% rownames(series_data)]
  
  if (length(gene_ids) == 0) {
    stop("None of the genes in the provided gene set were found in the expression matrix.")
  }
  
  # Step 3: 提取表达矩阵
  heatmap_data <- series_data[gene_ids, ]
  
  # 使用基因符号替代行名
  rownames(heatmap_data) <- gene_symbols[gene_ids %in% gene_info$ENTREZID]
  
  # Step 4: 按行（基因）标准化
  heatmap_data <- t(scale(t(heatmap_data)))
  heatmap_data[is.na(heatmap_data)] <- 0  # 处理 NA
  
  # Step 5: 样本排序
  colData$condition <- factor(colData$condition, levels = unique(colData$condition))
  sorted_indices <- order(colData$condition)
  
  # 重新排序表达矩阵和分组信息
  heatmap_data <- heatmap_data[, sorted_indices]
  annotation_col <- data.frame(Condition = colData$condition[sorted_indices])
  rownames(annotation_col) <- colnames(heatmap_data)
  
  # Step 6: 自定义颜色
  annotation_colors <- list(
    Condition = c(
      "Control" = "#F1C40F", "STAM_6w" = "#9B59B6", 
      "STAM_8w" = "#2980B9", "STAM_12w" = "#1ABC9C", 
      "STAM_20w_Tumor" = "#C0392B", "STAM_20w_NonTumor" = "#27AE60"
    )
  )
  
  # Step 7: 绘制热图
  pheatmap(
    heatmap_data,
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


