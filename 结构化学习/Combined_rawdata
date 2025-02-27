

# 设置文件路径（确保路径正确）
file_path <- "~/R/RNA/data/NASH/GSE114261/Rawdata"
# 获取所有的 .gz 文件
files <- list.files(file_path, pattern = ".*\\.gz$", full.names = TRUE)
# 初始化一个空数据框用于存储数据
combined_data <- NULL
# 读取第一个文件并设置行名
first_file <- files[1]
first_data <- read.table(first_file, header = FALSE, row.names = 1, sep = "\t", fill = TRUE)
# 提取样本名（第一个 `_` 之前的部分）
sample_name <- sub("_.*", "", basename(first_file))  # 获取文件名并去掉第一个 `_` 之后的部分
colnames(first_data) <- sample_name
# 将第一个文件的数据初始化到 combined_data 中
combined_data <- first_data
# 遍历其余文件并按基因名顺序拼接
for (file in files[-1]) {
  # 读取每个文件
  data <- read.table(file, header = FALSE, row.names = 1, sep = "\t", fill = TRUE)
  # 确保 data 是数据框
  if (is.vector(data)) {
    data <- as.data.frame(data)
  }
  # 提取样本名（第一个 `_` 之前的部分）
  sample_name <- sub("_.*", "", basename(file))  # 获取文件名并去掉第一个 `_` 之后的部分
  colnames(data) <- sample_name
  # 拼接数据
  combined_data <- cbind(combined_data, data)
}

