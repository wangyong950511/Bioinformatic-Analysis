## 
genes = [
    "ISG15", "IFIT2", "IFIT3", "KLRB1", "SLC4A10", "TMIGD2", "IL7R",
    "TOB1", "TCF7", "GZMK", "NR4A1", "TNFSF9", "GZMH", "GZMA",
    "CD14", "ZNF683", "CXCR6", "XCL1", "XCL2", "PDCD1", "CXCL13",
    "HAVCR2", "TIGIT", "LAG3", "LAYN", "ENTPD1", "STMN1", "MKI67"
]

# 聚类标签（如有多个，请确认使用的列名）
groupby_key = "leiden_CD8T_0.5"  

# 绘制气泡图
sc.pl.dotplot(
    adata_CD8T,
    var_names=genes,
    groupby=groupby_key,
    standard_scale="var",     # 每行基因标准化（颜色更易比较）
    dot_max=0.7,              # 调整最大圆（默认是0.5），更接近图中视觉
    dot_min=0.05,             # 最小圆点
    color_map="RdBu_r",          # 使用蓝白红渐变色
    figsize=(12, 3),
    dendrogram=False,
    swap_axes=False,
    save="my_dotplot.pdf"
)
