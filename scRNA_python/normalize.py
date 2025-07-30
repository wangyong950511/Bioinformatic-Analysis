# Harmony去批次

# 降维聚类
sc.pp.neighbors(adata_CD8T, use_rep="X_pca_harmony")
sc.tl.umap(adata_CD8T)
sc.tl.leiden(adata_CD8T, resolution=0.5, key_added="leiden_CD8T_0.5")
sc.pl.umap(adata_CD8T, color=["leiden_CD8T_0.5"])
