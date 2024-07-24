import subprocess
import pandas as pd
import argparse
import os

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="处理 emapper 注释结果并创建 OrgDB")
    parser.add_argument("emapper_anno", help="emapper 注释结果文件路径")
    parser.add_argument("output_dir", help="存储输出的目录")
    parser.add_argument("proteins_file", help="蛋白质 FASTA 文件路径")
    return parser.parse_args()

def preprocess_data(emapper_file):
    """读取并预处理 emapper 注释结果"""
    emapper = pd.read_csv(emapper_file, sep="\t", comment='#', header=None, 
                          usecols=[0, 6, 7, 8, 9, 11, 12],  # 仅选择必要的列
                          names=['GID', 'COG', 'Gene_Name', 'Gene_Symbol', 'GO', 'KO', 'Pathway'])
    print("数据加载完毕，显示前几行和列索引：")
    print(emapper.head())
    print("列索引：", emapper.columns)
    return emapper

def call_r_script(cleaned_data_path, output_dir, proteins_file):
    """调用 R 脚本处理 OrgDB 创建"""
    print("使用以下参数调用 R 脚本：", cleaned_data_path, output_dir, proteins_file)
    subprocess.run(["Rscript", "create_orgdb.R", cleaned_data_path, output_dir, proteins_file], check=True)

def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)
    cleaned_data_path = os.path.join(args.output_dir, "cleaned_data.csv")
    emapper_data = preprocess_data(args.emapper_anno)
    emapper_data.to_csv(cleaned_data_path, index=False)
    print(f"数据已保存至：{cleaned_data_path}")
    
    if os.path.exists(cleaned_data_path) and os.path.exists(args.proteins_file):
        print("文件存在，准备进行 R 处理。")
        call_r_script(cleaned_data_path, args.output_dir, args.proteins_file)
    else:
        print("文件创建失败或蛋白质文件不存在。")

if __name__ == "__main__":
    main()
