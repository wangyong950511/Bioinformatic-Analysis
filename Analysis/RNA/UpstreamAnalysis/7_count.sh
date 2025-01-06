#!/bin/bash

# 输出计数结果文件名
output_file="gene_counts.txt"

# 使用 featureCounts 计算基因计数，启用多线程
featureCounts -T 20 \
              -a ~/Documents/Refdata/Rnaseq/GCF_000001635.27_GRCm39_genomic.gtf \
              -o "$output_file" \
              -p --countReadPairs -t exon -g gene_id *_sorted.bam

echo "Counted features for all BAM files and saved to $output_file"
