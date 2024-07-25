import pandas as pd
import matplotlib.pyplot as plt
import argparse

def process_files(input_files):
    # 初始化一个空的数据框
    combined_data = pd.DataFrame()

    # 读取每个文件并添加分组标签
    for i, input_file in enumerate(input_files):
        data = pd.read_csv(input_file, sep='\t')
        data['group'] = f'Group_{i+1}'
        combined_data = pd.concat([combined_data, data], ignore_index=True)
    
    return combined_data

def calculate_isoform_categories(data):
    # 计算每个基因的转录本数量
    gene_isoform_count = data.groupby(['group', 'associated_gene']).size().reset_index(name='isoform_count')

    # 将isoform数量分为四个类别
    gene_isoform_count['isoform_category'] = pd.cut(gene_isoform_count['isoform_count'], bins=[0, 1, 2, 3, float('inf')], labels=['1', '2', '3', '4 or more'])

    # 计算每个分组中每个类别的基因数量
    group_summary = gene_isoform_count.groupby(['group', 'isoform_category'], observed=True).size().reset_index(name='gene_count')

    return group_summary

def plot_summary(group_summary, output_file=None):
    # 创建条形图
    plt.figure(figsize=(10, 6))
    categories = ['1', '2', '3', '4 or more']
    groups = group_summary['group'].unique()
    
    # 将数据转换为绘图所需的格式
    plot_data = group_summary.pivot(index='group', columns='isoform_category', values='gene_count').reindex(columns=categories)

    # 绘制条形图
    plot_data.plot(kind='bar', stacked=True, figsize=(10, 6), width=0.7, color=['skyblue', 'lightgreen', 'orange', 'lightcoral'])

    # 添加数据标签
    for i, row in group_summary.iterrows():
        group_index = list(groups).index(row['group'])
        bar_positions = [p + group_index * 0.7 for p in range(len(categories))]
        bar_offset = bar_positions[categories.index(row['isoform_category'])]
        plt.text(bar_offset, row['gene_count'], f'{row["gene_count"]}', ha='center', va='bottom', fontsize=9)

    # 自定义图表
    plt.xlabel('Group')
    plt.ylabel('Number of Genes')
    plt.title('Number of Genes with Different Number of Isoforms by Group')
    plt.legend(title='Isoform Category', loc='upper right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # 保存或显示图表
    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Process multiple final_classification.txt files and plot results.")
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input file paths')
    parser.add_argument('-o', '--output', required=False, help='Output file path for the plot')
    parser.add_argument('-s', '--summary', required=True, help='Output file path for the summary statistics')
    args = parser.parse_args()

    # 处理文件
    combined_data = process_files(args.input)
    
    # 计算转录本类别
    group_summary = calculate_isoform_categories(combined_data)
    
    # 打印并保存统计结果
    print(group_summary)
    group_summary.to_csv(args.summary, index=False, sep='\t')
    
    # 绘制图表
    plot_summary(group_summary, args.output)

if __name__ == '__main__':
    main()

