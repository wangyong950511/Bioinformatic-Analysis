library(DESeq2)
library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(pheatmap)
library(org.Hs.eg.db)
library(data.table)
library(myplotColors)




# 模块化代码
## 一、导入矩阵代码
# 阅读结构
readLines("~/R/RNA/data/NASH/GSE287943/GSE287943_countData.txt.gz", n = 50)
### 常规
countData <- read.table("~/R/RNA/data/NASH/GSE287943/GSE287943_countData.txt.gz", skip=0,header = TRUE, row.names = 1)
### 变种1-读取xlsx文件
countData <- read_excel("~/R/RNA/data/NASH/GSE253217/GSE253217_Liver_tmm_normalized_counts.xlsx")
### 变种1-读取tsv文件
countData <- read_tsv("~/R/RNA/data/NASH/GSE236832/GSE236832_de_summary.tsv.gz")
### 变种1-读取csv文件
countData <- read.csv(gzfile("~/R/RNA/data/NASH/GSE230745/GSE230745_kidney_chp_raw_counts.csv.gz"), header = TRUE, row.names = 1)
### 万能解决
countData <- read.delim(gzfile("~/R/RNA/data/NASH/GSE275069/GSM8465145_Matrix_HF_vs_ND.txt.gz"),header = TRUE,sep = "\t",fill = TRUE)
### 自己跑上游
rawData <- read.table("~/R/RNA/data/NASH/GSE160292/gene_counts.txt", skip=0,header = TRUE, row.names = 1)
countData <- rawData[, -(1:5)]
colnames(countData) <- sub("_.*", "", colnames(countData))



## 二、处理矩阵
# 获取基因名映射
gene_mapping <- AnnotationDbi::select(
  org.Mm.eg.db,                         # 人物种为Hs，鼠为Mm
  keys = rownames(countData),           # 使用 countData 的行名（ENSG ID）
  columns = c("SYMBOL"),                # 想要的列：基因名
  keytype = "ENSEMBL"                   # 关键列类型：ENTREZID，ACCNUM
)

### 常规
# 输出重复值
duplicated_symbols <- unique(gene_mapping$SYMBOL[duplicated(gene_mapping$SYMBOL)])
# 去除 SYMBOL 中为空值的行
gene_mapping <- gene_mapping[!is.na(gene_mapping$SYMBOL), ]
countData$SYMBOL <- gene_mapping$SYMBOL[match(rownames(countData), gene_mapping$ENSEMBL)]
# 去除 SYMBOL 为 NA 的行
countData <- countData[!is.na(countData$SYMBOL), ]
# 合并重复基因名（按 SYMBOL 分组叠加值）
countData <- aggregate(. ~ SYMBOL, data = countData, sum)
# 将 SYMBOL 列作为行名
rownames(countData) <- countData$SYMBOL
countData <- countData[, -1]  # 删除 SYMBOL 列


### 更改列名及行名（如果需要）
#常规
colnames(countData) <- gsub("Sample_", "S", colnames(countData))
##保留第一个_之前的文字
colnames(countData) <- sub("_.*", "", colnames(countData))
# 保留第二个_之前的文字
colnames(countData) <- sub("^([^_]*_[^_]*)_.*", "\\1", colnames(countData))
# 去除行名
rownames(geo_selected) <- sub("\\..*", "", rownames(geo_selected))




## 二、读取分组信息
# 1、阅读series_matrix文件
# 文件结构
readLines("~/R/RNA/data/NASH/GSE287943/GSE287943_series_matrix.txt.gz", n = 50)
# 读取文件
geo_data <- read.delim("~/R/RNA/data/NASH/GSE287943/GSE287943_series_matrix.txt.gz", header = TRUE, sep = "\t", skip = 28)
# 转置数据（确保是字符矩阵）
geo_transposed <- as.data.frame(t(geo_data), stringsAsFactors = FALSE)
geo_transposed <- geo_transposed[-1, ]  # 删除第一行
# 选取需要列（第一个为列名）
geo_selected <- geo_transposed[, c(9), drop = FALSE]
rownames(geo_selected) <- geo_selected[, 1]  # 将第一列设为行名
geo_selected <- geo_selected[, -1, drop = FALSE]
colnames(geo_selected) <- c("condition")
# 更改分组名
geo_selected$condition<- gsub(" \\+5% mannose", "_5_mannose", geo_selected$condition)
geo_selected$condition<- gsub(" \\+20% mannose", "_20_mannose", geo_selected$condition)
geo_selected$condition<- gsub("treatment: Normal Diet", "NC", geo_selected$condition)
geo_selected$condition<- gsub("treatment: FAT-MASH Diet", "MASH", geo_selected$condition)
geo_selected$condition<- gsub(" \\((therapy|Therapy)\\)", "_Therapy", geo_selected$condition)
geo_selected$condition <- factor(geo_selected$condition)


# a.获取 countData 的列名并对geo_selected进行排序（如果需要）
sample_order <- colnames(countData)
geo_selected <- geo_selected[match(sample_order, rownames(geo_selected)), , drop = FALSE]
# b.修改行名（如果需要）
rownames(geo_selected) <- paste0("x", rownames(geo_selected))




# 2、新建样本信息表
geo_selected <- data.frame(
  row.names = colnames(countData),
  condition = factor(c(
    rep("aaa", 4),
    rep("bbb", 11),
    rep("ccc", 5),
    rep("ddd", 6),
    rep("eee", 6)
  ), levels = c("aaa", "bbb", "ccc", "ddd", "eee"))
)






## 二、提取差异数组（均一化处理）
# 提取样本计数数据
count_matrix <- as.matrix(countData)  # 修剪数据
# 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = geo_selected,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) > 10, ]
# 运行差异表达分析
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

colData_GSEnnn <- geo_selected
dds_GSEnnn <- dds
save(dds_GSEnnn,colData_GSEnnn, file = "~/R/RNA/data/NASH/GSEnnn/Data_GSEnnn.RData")


## 带logFC的数据图
# 获取数据
load("~/R/RNA/data/NASH/GSEnnn/Data_GSEnnn.RData")
analyze_gene_expression_GSEnnn <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }

  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSEnnn"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}
analyze_gene_expression_GSEnnn("COL1A1",dds_GSEnnn,colData_GSEnnn)


