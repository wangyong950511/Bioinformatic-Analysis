


analyze_gene_expression_GSE154892 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionCDAHFD_12W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "CDAHFD_12W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "CDAHFD_12W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_12W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("Expression of", GeneX),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  return(plot)
}
