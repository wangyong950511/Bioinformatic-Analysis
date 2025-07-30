# 1. 差异表达分析
sc.tl.rank_genes_groups(
    adata_CD8T,
    groupby="leiden_CD8T_0.5",
    method="wilcoxon",         # 你也可以选择 't-test', 'logreg'
    use_raw=True,
    key_added="rank_genes_CD8T"
)

# 2. 导出 top 20 差异基因为 CSV
sc.get.rank_genes_groups_df(
    adata_CD8T,
    group=None,                # 所有群
    key="rank_genes_CD8T",
    pval_cutoff=0.05
).groupby("group").head(20).to_csv("CD8T_marker_top20.csv", index=False)

print("Top 20 差异基因已保存为 CD8T_marker_top20.csv")
