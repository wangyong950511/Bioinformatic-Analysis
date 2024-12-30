library(DESeq2)
library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)



####手动调节数据####
# 数据路径
predata <- data.frame(newdata)
# 横纵坐标范围
ylim=200
xlim=6
# 显示显著基因个数
num_label=10
# 筛选FC范围
FoldChange=0.8


####具体过程####
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
  xlab = "log2FC",
  ylab = "-log10(Adjust P-value)"
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
    point.padding = 0.5, # 点与标签之间的间距
    max.overlaps = 50,   # 最多允许多少个标签
    force = 1.5          # 调整标签之间的排斥力
  )
