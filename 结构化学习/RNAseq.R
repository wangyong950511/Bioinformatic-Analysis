# 模块化代码
## 一、导入矩阵代码
# 阅读结构
readLines("~/R/RNA/data/NASH/GSE287943/GSE287943_countData.txt.gz", n = 50)
### 常规
countData <- read.table("~/R/RNA/data/NASH/GSE287943/GSE287943_countData.txt.gz", header = TRUE, row.names = 1)
### 变种1-跳过
countData <- read.table("~/R/RNA/data/NASH/GSE287943/GSE287943_countData.txt.gz", skip=1,header = TRUE, row.names = 1)
## 二、处理矩阵
### 常规
\# 输出重复值
duplicated_symbols <- unique(gene_mapping$SYMBOL[duplicated(gene_mapping$SYMBOL)])
\# 去除 SYMBOL 中为空值的行
gene_mapping <- gene_mapping[!is.na(gene_mapping$SYMBOL), ]  
countData$SYMBOL <- gene_mapping$SYMBOL[match(rownames(countData), gene_mapping$ENSEMBL)]
\# 去除 SYMBOL 为 NA 的行
countData <- countData[!is.na(countData$SYMBOL), ]
\# 合并重复基因名（按 SYMBOL 分组叠加值）
countData <- aggregate(. ~ SYMBOL, data = countData, sum)
\# 将 SYMBOL 列作为行名
rownames(countData) <- countData$SYMBOL
countData <- countData[, -1]  # 删除 SYMBOL 列
### 更改列名
#### 常规
colnames(countData) <- gsub("Sample_", "S", colnames(countData))  # 将 "Sample_" 替换为 "S"
#### 保留第一个_之前的文字
colnames(countData) <- sub("_.\*", "", colnames(countData))
#### 保留第二个_之前的文字
colnames(countData) <- sub("^([^_]\*_[^_]\*)_.\*", "\\1", colnames(countData))



# 阅读结构
readLines("~/R/RNA/data/NASH/GSE287943/GSE287943_countData.txt.gz", n = 50)



# 定义文件路径
file_path <- "~/R/RNA/data/NASH/GSE287943/GSE287943_series_matrix.txt.gz"
# 读取文件
geo_data <- read.delim(file_path, header = TRUE, sep = "\t", skip = 28)
# 转置数据（确保是字符矩阵）
geo_transposed <- as.data.frame(t(geo_data), stringsAsFactors = FALSE)
geo_transposed <- geo_transposed[-1, ]  # 删除第一行
# 选取需要列（第一个为列名）
geo_selected <- geo_transposed[, c(22,13)]
rownames(geo_selected) <- geo_selected[, 1]  # 将第一列设为行名
geo_selected <- geo_selected[, -1, drop = FALSE]
colnames(geo_selected) <- c("condition")
# 更改分组名
geo_selected$condition<- gsub(" \\+5% mannose", "_5_mannose", geo_selected$condition)
geo_selected$condition<- gsub(" \\+20% mannose", "_20_mannose", geo_selected$condition)
geo_selected$condition<- gsub("treatment: Normal Diet", "NC", geo_selected$condition)
geo_selected$condition<- gsub("treatment: FAT-MASH Diet", "MASH", geo_selected$condition)
geo_selected$condition<- gsub(" \\((therapy|Therapy)\\)", "_Therapy", geo_selected$condition)
