# -------------------------------
# Script Name: NASH_Fibrosislevel_GSE162694.R
# Data：GSE162694
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
library(data.table)


####导入count矩阵和样本信息####
countData <- read.csv("~/R/RNA/data/NASH/GSE162694/GSE162694_raw_counts.csv.gz", header = TRUE, row.names = 1)
# 提取样本计数数据
count_matrix <- as.matrix(countData)  # 修剪数据


# 样本信息表
# 读取文件并跳过前 24 行
file_path <- "~/R/RNA/data/NASH/GSE162694/GSE162694_series_matrix.txt"
geo_data <- readLines(file_path)
geo_data <- geo_data[-(1:24)]
# 提取以 "!Sample_" 开头的行
sample_lines <- grep("^!Sample_", geo_data)
sample_metadata <- geo_data[sample_lines]
# 拆分为列
sample_split <- strsplit(sample_metadata, "\t")
field_names <- sapply(sample_split, function(x) x[1])
field_values <- lapply(sample_split, function(x) x[-1])
# 转置为数据框
sample_df <- as.data.frame(t(do.call(rbind, field_values)), stringsAsFactors = FALSE)
colnames(sample_df) <- gsub("^!Sample_", "", field_names)
# 去除所有列中的 `"` 符号
sample_df <- as.data.frame(lapply(sample_df, function(x) gsub('"', '', x)), stringsAsFactors = FALSE)
# 修改 title 列，去除 "nash" 之前的所有字符
sample_df$title <- gsub("^.*(?=nash)", "X548", sample_df$title, perl = TRUE)
# 提取感兴趣的列
selected_columns <- sample_df[, c("title", "characteristics_ch1.3")]
colnames(selected_columns) <- c("Title", "condition")
# 提取分组信息
selected_columns$Group <- sub(" .*", "", selected_columns$Title)
# 修改 condition 列，去除 "fibrosis stage: "
selected_columns$condition <- gsub("fibrosis stage: ", "", selected_columns$condition)
selected_columns$condition <- gsub("normal liver histology", "NC", selected_columns$condition)
# 确保 condition 列为因子类型
selected_columns$condition <- factor(selected_columns$condition)
# 创建 colData 数据框
colData <- data.frame(
  condition = selected_columns$condition,  # 提取 condition 列
  row.names = selected_columns$Title       # 将 Title 列作为行名
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

# 输出差异基因
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
analyze_gene_expression("P4HA1", vsd, colData)
gene_set <- c(
  # 细胞外基质 (ECM) 重塑
  "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", "COL6A1",
  "FN1", "MMP1", "MMP2", "MMP3", "MMP9", "MMP13",
  "TIMP1", "TIMP2", "LAMA2", "SPARC",
  
  # 成纤维细胞激活 (Fibroblast Activation)
  "ACTA2", "TAGLN", "POSTN", "PDGFRA", "PDGFRB",
  "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "CTGF",
  
  # 炎症 (Inflammation)
  "IL6", "IL1B", "TNF", "CXCL8", "CCL2", "CCL3", "CCL5", "CCR2",
  "ICAM1", "VCAM1", "PTGS2", "CD68", "CD163",
  
  # 肝星状细胞 (Hepatic Stellate Cells) 标志物
  "PDGFRB", "LRAT", "GFAP", "VIM", "RBP1", "ALDH1A1",
  
  # 纤维化因子 (Fibrogenic Factors)
  "TGFB1", "TGFB2", "TGFB3", "PDGFB", "VEGFA", "ANGPT2", 
  "LOX", "LOXL2", "ELN", "SERPINE1",
  
  # 其他 ECM 调控基因
  "ITGA5", "ITGB1", "THBS1", "THBS2", "FBN1", "HAS2"
)
generate_ddr_heatmap(count_matrix, colData, gene_set, "Heatmap of DDR Genes")

####特定基因分析（预运行）####
# 定义函数
analyze_gene_expression <- function(GeneX, vsd, colData) {
  # Step 1: 查询基因的 NCBI Gene ID (ENSEMBL)
  gene_info <- select(org.Hs.eg.db,
                      keys = GeneX,
                      columns = c("ENSEMBL", "SYMBOL"),
                      keytype = "SYMBOL")
  
  # 检查查询结果
  if (nrow(gene_info) == 0 || is.na(gene_info$ENSEMBL[1])) {
    stop("目标基因未找到对应的 ENSEMBL，请检查基因名是否正确。")
  }
  
  # 提取 ENSEMBL
  gene_id <- gene_info$ENSEMBL[1]
  
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
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "0", "1", "2", "3", "4"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 5: 绘制箱线图
  custom_colors_fill <- c("NC" = "white", "0" = "white", "1" = "white", "2" = "white", "3" = "white", "4" = "white")
  custom_colors_color <- c("NC" = "#98DF8A", "0" = "#FFB347", "1" = "#FF8C42", "2" = "#FF6F31", "3"  = "#D62728", "4" = "#A51220")
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  
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
    columns = c("ENSEMBL", "SYMBOL"), 
    keytype = "SYMBOL"
  )
  
  # 去除没有匹配到的基因
  gene_info <- gene_info[!is.na(gene_info$ENSEMBL), ]
  
  # 提取 ENTREZ ID
  gene_ids <- gene_info$ENSEMBL
  
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
  id_to_symbol <- setNames(gene_info$SYMBOL, gene_info$ENSEMBL)
  rownames(heatmap_data_scaled) <- id_to_symbol[rownames(heatmap_data_scaled)]
  
  # 检查替换后的行名
  if (any(is.na(rownames(heatmap_data_scaled)))) {
    stop("Some gene SYMBOLs could not be mapped. Please check your input data.")
  }
  
  # Step 7: 自定义颜色
  annotation_colors <- list(
    Condition = c(
      "NC" = "#98DF8A",  # 浅绿色
      "0"  = "#FFB347",  # 浅橙色
      "1"  = "#FF8C42",  # 橙色稍深
      "2"  = "#FF6F31",  # 深橙色
      "3"  = "#D62728",  # 红色
      "4"  = "#A51220"   # 深红色
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




