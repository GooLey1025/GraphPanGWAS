import pandas as pd
import argparse

print("运行表型处理 process_phenotype.py")

# 解析命令行参数
parser = argparse.ArgumentParser(description="根据样本顺序排序数据文件，并格式化输出")
parser.add_argument('map_file', type=str, help="包含样本顺序的文件")
parser.add_argument('csv_file', type=str, help="要排序的CSV文件（带ID列）")
parser.add_argument('output_file', type=str, help="输出的排序并格式化后的文件路径")
args = parser.parse_args()

# 读取样本顺序
with open(args.map_file, "r") as f:
    sample_order = [line.strip() for line in f.readlines()]

# 读取表型数据
data = pd.read_csv(args.csv_file, sep='\t')

# 确保 ID 为字符串类型
data['ID'] = data['ID'].astype(str)

# 按样本顺序重新排序
ordered_data = data.set_index('ID').reindex(sample_order).reset_index()

# 写出文件，缺失值替换为 -9，浮点数保留 2 位，确保 -9 是整数
with open(args.output_file, 'w') as f:
    # 写入表头
    f.write('\t'.join(ordered_data.columns) + '\n')
    
    # 写入每一行
    for _, row in ordered_data.iterrows():
        line_values = []
        for col, val in row.items():
            if col == 'ID':
                line_values.append(str(val))
            elif pd.isna(val):
                line_values.append('-9')  # 缺失值：整数 -9
            elif isinstance(val, int) or val == -9:
                line_values.append(str(int(val)))  # 保持整数
            else:
                line_values.append(f"{val:.2f}")  # 浮点数保留两位小数
        f.write('\t'.join(line_values) + '\n')

print(f"✅ 数据已排序并保存为 '{args.output_file}'，缺失值为 -9，浮点保留两位小数")
