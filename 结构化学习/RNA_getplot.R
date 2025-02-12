# 获取数据
load("~/R/RNA/data/NASH/GSE287943/Data_GSE287943.RData")


# vsd单组简易函数
analyze_gene_expression_GSEnnn <- function(GeneX, vsd, colData) {
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
  
  # Step 3: 仅保留 "NC" 和 "MASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MASH"))
  
  # Step 5: 按因子水平排序数据
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 6: 统计分析（仅 NC vs. MASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 7: 自定义颜色（动态生成）
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("Expression of", GeneX),
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
analyze_gene_expression_GSEnnn("Trip6", vsd_GSEnnn, colData_GSEnnn)
