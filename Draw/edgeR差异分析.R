# 安装statmod包
install.packages("statmod")

# 安装edgeR包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

# 载入所需包
library(statmod)
library(edgeR)

# 读取数据
rawdata <- read.csv("/Users/wangyong/Downloads/read_counts/case01data_filtered.csv", header = TRUE, row.names = 1)

# 查看前六行
head(rawdata)

# 创建edgeR储存数据的DGElist对象
y <- DGEList(counts = rawdata[, 5:10], genes = rawdata[, 1:4])
y

# 按照基因count之和，对数据进行降序排列
o <- order(rowSums(y$counts), decreasing = TRUE)
y <- y[o, ]

# 查看排序后的结果
head(y$counts)
tail(y$counts)

# 对数据进行TMMnormalization，计算标准化因子
y <- calcNormFactors(y)
y$samples

# 检查样本中的异常值
plotMDS(y)

# 生成试验设计矩阵
Patient <- factor(c(8, 8, 33, 33, 51, 51))
Tissue <- factor(c("N", "T", "N", "T", "N", "T"))
data.frame(Sample = colnames(y), Patient, Tissue)

# 指定blocking factor和比较组，生成design矩阵
# 这种加性模型适用于配对实验设计或具有批次效应的实验
design <- model.matrix( ~ Patient + Tissue)
rownames(design) <- colnames(y)
design

# 可以用BCV plot查看离散系数
plotBCV(y)

# 估计离散系数
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion

# 先进行glm拟合
fit <- glmFit(y, design)
colnames(fit)

# 然后执行likelihood ratio test，对肿瘤与正常组织进行差异性分析
lrt <- glmLRT(fit)

# 展示top10显著表达的差异基因
topTags(lrt)

# 根据p值对分析结果升序排列
o <- order(lrt$table$PValue)

# 查看top10基因的CPM值
CPM <- as.data.frame(cpm(y)[o[1:10], ])
CPM

# 尝试计算log2FC值，结果有细微差异
N <- (CPM$N8 + CPM$N33 + CPM$N51)
Tum <- (CPM$T8 + CPM$T33 + CPM$T51)
log2(Tum) - log2(N)

# 以FDR=0.05为阈值，统计差异基因个数
summary(decideTests(lrt))

# 绘制差异基因散点图
plotMD(lrt)
abline(h = c(-1, 1), col = "blue")