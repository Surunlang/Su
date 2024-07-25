import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def read_and_merge_data(classification_file, coding_potential_file):
    # 读取第一个文件
    data1 = pd.read_csv(classification_file, sep='\t')

    # 筛选出新发现的异构体
    new_isoforms = data1[data1['structural_category'].str.contains('novel_')]

    # 读取第二个文件
    data2 = pd.read_csv(coding_potential_file, sep='\t')

    # 合并两个数据表，匹配sequence_name
    merged_data = pd.merge(new_isoforms, data2, left_on='isoform', right_on='sequence_name')

    return merged_data

def plot_coding_noncoding(merged_data, output_dir):
    # 筛选出coding和noncoding的数据
    coding_data = merged_data[merged_data['classification'] == 'coding']
    noncoding_data = merged_data[merged_data['classification'] == 'noncoding']

    # 绘制编码和非编码异构体的数量
    labels = ['Coding', 'Noncoding']
    sizes = [len(coding_data), len(noncoding_data)]
    colors = ['#ff9999', '#66b3ff']
    explode = (0.1, 0)  # "explode" the 1st slice

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.title('Coding vs Noncoding Isoforms')
    plt_path = os.path.join(output_dir, 'coding_vs_noncoding_isoforms.pdf')
    plt.savefig(plt_path)
    plt.show()

    return coding_data, noncoding_data

def save_data(coding_data, noncoding_data, output_dir):
    # 保存编码和非编码异构体的数据到文件
    coding_path = os.path.join(output_dir, 'coding_isoforms.csv')
    noncoding_path = os.path.join(output_dir, 'noncoding_isoforms.csv')

    coding_data[['isoform', 'associated_gene']].to_csv(coding_path, index=False)
    noncoding_data[['isoform', 'associated_gene']].to_csv(noncoding_path, index=False)

def main():
    parser = argparse.ArgumentParser(description='Process and plot coding vs noncoding isoforms.')
    parser.add_argument('-i', '--input', required=True, help='Input classification file')
    parser.add_argument('-t', '--tsv', required=True, help='Input coding potential file')
    parser.add_argument('-o', '--output', required=True, help='Output directory for results')

    args = parser.parse_args()

    classification_file = args.input
    coding_potential_file = args.tsv
    output_dir = args.output

    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    merged_data = read_and_merge_data(classification_file, coding_potential_file)
    coding_data, noncoding_data = plot_coding_noncoding(merged_data, output_dir)
    save_data(coding_data, noncoding_data, output_dir)

if __name__ == "__main__":
    main()

