#!/bin/bash

input_vcf="$1"
output_vcf="$2"

if [[ ! -f "$input_vcf" ]]; then
    echo "Error: Input file not found!"
    exit 1
fi

# 处理 VCF 文件，提取 INFO 字段最后一个 '-' 后的数字，筛选 >=50 的行
awk -F'\t' '
BEGIN { OFS="\t" }
{
    if ($0 ~ /^#/) { print; next }  # 直接打印 VCF 头部
    split($8, arr, "-")  # 按 '-' 拆分 INFO 字段
    num = arr[length(arr)]  # 获取最后一个数字
    if (num >= 50) print
}' "$input_vcf" > "$output_vcf"

echo "Filtering complete. Results written to: $output_vcf"
