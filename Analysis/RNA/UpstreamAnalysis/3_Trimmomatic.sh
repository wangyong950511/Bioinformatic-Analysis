#!/bin/bash

# Trimmomatic 安装路径（确保 Trimmomatic 在 PATH 中，或指定完整路径）
TRIMMOMATIC_PATH="/opt/Trimmomatic-0.39/trimmomatic-0.39.jar"

# 适配器序列文件路径
ADAPTERS="/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# 循环遍历所有输入文件，匹配 *_1.fq.gz 和 *_2.fq.gz 的配对
for forward_file in *_1.fq.gz; do
    # 确定相应的 reverse 文件名
    reverse_file=${forward_file/_1/_2}
    
    # 检查反向文件是否存在
    if [[ -f "$reverse_file" ]]; then
        # 设置输出文件名
        output_forward="${forward_file%.fq.gz}_trimmed.fq.gz"
        output_reverse="${reverse_file%.fq.gz}_trimmed.fq.gz"
        unpaired_forward="${forward_file%.fq.gz}_unpaired.fq.gz"
        unpaired_reverse="${reverse_file%.fq.gz}_unpaired.fq.gz"

        # 运行 Trimmomatic
        java -jar $TRIMMOMATIC_PATH PE -phred33 \
            "$forward_file" "$reverse_file" \
            "$output_forward" "$unpaired_forward" \
            "$output_reverse" "$unpaired_reverse" \
            ILLUMINACLIP:$ADAPTERS:2:30:10 \
            SLIDINGWINDOW:4:20 \
            MINLEN:50

        echo "Trimmed: $forward_file and $reverse_file"
    else
        echo "Reverse file for $forward_file not found!"
    fi
done