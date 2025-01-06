#!/bin/bash


# 设置索引路径
index_path="~/Documents/Refdata/Rnaseq/GRCm39_index"

# 遍历每个以 "_1.fastq.gz" 结尾的文件
for file in *"_1.fq.gz"; do
    # 生成对应的 _2 文件名
    file2="${file/_1/_2}"
    
    # 检查对应的 _2 文件是否存在
    if [[ -f "$file2" ]]; then
        # 生成输出文件名
        output_file="${file/_1.fq.gz/.sam}"
        
        # 运行 HISAT2
        hisat2 -x "$index_path" -1 "$file" -2 "$file2" -S "$output_file" -p 20
    else
        echo "对应的文件 $file2 不存在，跳过 $file"
    fi
done