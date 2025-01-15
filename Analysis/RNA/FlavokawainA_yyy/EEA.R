# -------------------------------
# Script Name: EEA.R
# Data：Self
# Purpose: 分析Flavokawain A效果
# Author: WangYong
# Date: 2025-1-10



# 加载库和设置 ------------------------------------
library(readxl)
library(dplyr)
library(gdata)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(DESeq2)
library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)

setwd("/home/drwang/R/Workdata")
# 设置参数 ------------------------------------
logFC_limit=1
####导入样本信息####
# 指定文件路径
file1 <- "/home/drwang/Desktop/PC_NC_vs_PC_EFA_all.csv"
file2 <- "/home/drwang/Desktop/HCC_NC_vs_HCC_EFA_all.csv"
# 读取 CSV 文件
df1 <- read.csv(file1, header = TRUE, sep = ",")
df2 <- read.csv(file2, header = TRUE, sep = ",")


# 筛选满足条件的行：|logFC| ≥ 0.5 且 P.adj < 0.05
filtered_df1 <- df1 %>% filter(abs(logFC) >= logFC_limit & P.adj < 0.05)
filtered_df2 <- df2 %>% filter(abs(logFC) >= logFC_limit & P.adj < 0.05)


# 取交集（按 Gene 列）
common_genes <- inner_join(filtered_df1, filtered_df2, by = "SYMBOL")
symbols <- common_genes$SYMBOL




####绘制火山图准备####
## 手动参数
# 横纵坐标范围
ylim=30
xlim=6
# 显示显著基因个数
num_label=4
# 筛选FC范围
FoldChange=1
## 具体过程
# 取交集（按 Gene 列）
common_genes <- inner_join(df1, df2, by = "SYMBOL")
# 去除 SYMBOL 列中为空或 NA 的行
common_genes_df <- common_genes %>%
  filter(!is.na(SYMBOL) & SYMBOL != "-") %>%
  distinct(SYMBOL, .keep_all = TRUE)
# 创建 newdata 数据框
newdata <- data.frame(
  log2FoldChange = common_genes_df$logFC.x,  # 提取 logFC.x 列
  padj = common_genes_df$P.adj.x             # 提取 P.adj.x 列
)
# 将行名设置为 SYMBOL 列
rownames(newdata) <- common_genes_df$SYMBOL
# 数据路径
predata <- data.frame(newdata)
# 创建 Firehill 数据框
Firehill <- data.frame(
  padj = predata$padj,              # padj 列
  log2FoldChange = predata$log2FoldChange,  # log2FC 列
  gene = rownames(predata)          # 行名作为 gene 列
)
# 去除 log2FoldChange 或 padj 中包含 NA 的行
Firehill <- Firehill %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))
# 校正转换
Firehill$logP <- -log10(Firehill$padj)
# 使用 tanh 函数进行对称性变换
Firehill <- Firehill %>%
  mutate(log2FoldChange = log2FoldChange / (1 + abs(log2FoldChange) / xlim),
         logP = logP / (1 + abs(logP) / ylim),
  )
Firehill <- Firehill %>%
  mutate(Group = case_when(
    padj < 0.05 & log2FoldChange > FoldChange ~ "up",
    padj < 0.05 & log2FoldChange < -FoldChange ~ "down",
    TRUE ~ "none"
  ))
table(Firehill$Group)
# 添加一个新列 Label
Firehill$Label <- ""
Firehill <- Firehill %>%
  arrange(padj)
up_genes <- head(Firehill$gene[which(Firehill$Group == "up")], num_label)
down_genes <- head(Firehill$gene[which(Firehill$Group == "down")], num_label)
top10_genes <- c(as.character(up_genes), as.character(down_genes))
Firehill$Label[match(top10_genes, Firehill$gene)] <- top10_genes
# 绘制火山图
ggscatter(
  Firehill,
  x = "log2FoldChange",
  y = "logP",
  color = "Group",
  palette = c("#2f5688", "#BBBBBB", "#cc0000"),
  size = 1,
  xlab = "log2FC.A",
  ylab = "-log10(Adjust P-value).A"
) +
  theme_base() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-FoldChange, FoldChange), linetype = "dashed") +
  geom_text_repel(
    data = subset(Firehill, Label != ""),
    aes(label = Label, color = Group),  # 标签颜色与分组颜色一致
    size = 3,
    show.legend = FALSE, # 禁用图例标签
    box.padding = 0.5,   # 标签框与点之间的间距
    point.padding = 0.6, # 点与标签之间的间距
    max.overlaps = 50,   # 最多允许多少个标签
    force = 1.5          # 调整标签之间的排斥力
  )



####KEGG富集分析####
# 普通KEGG 
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


# 信号通路相关KEGG
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


####KEGG通路图####
gene_list <- common_genes
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
         pathway.id = "hsa04010",            # cAMP 信号通路的 KEGG ID
         species = "hsa",                 # 物种（hsa 代表人类）
         limit = list(gene = c(-3, 3)),   # 表达值的上下限
         low = "blue", high = "red", mid = "white",) # 颜色设置：上调为红，下调为蓝





