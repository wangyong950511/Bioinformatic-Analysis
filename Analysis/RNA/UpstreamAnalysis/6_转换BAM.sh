#!/bin/bash

# 临时文件目录
temp_dir="/media/drwang/新加卷/temp"

# 如果临时目录不存在，则创建它
mkdir -p "$temp_dir"

# 遍历所有的 SAM 文件，转换为排序后的 BAM 文件
for sam_file in *.sam; do
    # 提取文件名，去掉扩展名
    base_name="${sam_file%.sam}"
    
    # 转换为 BAM 并排序，使用指定的临时目录
    samtools view -bS "$sam_file" | samtools sort -o "${base_name}_sorted.bam" -T "$temp_dir/${base_name}_temp"
    
    echo "Converted and sorted: ${base_name}_sorted.bam"
done
