

# 获取差异表达结果
res <- results(dds, contrast = c("condition", "CDAHFD_7w", "NC"))
# 过滤显著差异表达基因 (FDR < 0.05)
res_filtered <- res[which(res$padj < 0.05), ]
# 按照 log2FoldChange 排序
res_filtered <- res_filtered[order(res_filtered$log2FoldChange, decreasing = TRUE), ]
# 导出结果到 CSV 文件
write.csv(as.data.frame(res_filtered), file = "DEGs_results.csv", row.names = TRUE)
