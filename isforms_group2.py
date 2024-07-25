import pandas as pd
import matplotlib.pyplot as plt
import argparse

def process_files(input_files):
    # 初始化一个空的数据框
    combined_data = pd.DataFrame()

    # 读取每个文件并添加分组标签
    for i, input_file in enumerate(input_files):
        data = pd.read_csv(input_file, sep='\t')
        data['group'] = data['structural_category'].apply(lambda x: 'Annotated' if x in ['full-splice_match', 'incomplete-splice_match'] else 'Novel')
        data['source'] = f'Group_{i+1}'
        combined_data = pd.concat([combined_data, data], ignore_index=True)
    
    return combined_data

def calculate_isoform_count(data):
    # 计算每个分组中的转录本数量
    isoform_count = data.groupby(['source', 'group']).size().reset_index(name='isoform_count')

    return isoform_count

def plot_summary(isoform_count, output_file=None):
    # 创建条形图
    plt.figure(figsize=(12, 6))
    
    # 绘制条形图
    isoform_count.pivot(index='source', columns='group', values='isoform_count').plot(kind='bar', figsize=(12, 6), width=0.8, color=['skyblue', 'orange'])
    
    # 添加数据标签
    for i, row in isoform_count.iterrows():
        plt.text(i, row['isoform_count'], f'{int(row["isoform_count"])}', ha='center', va='bottom', fontsize=9, rotation=90)

    # 自定义图表
    plt.xlabel('Source')
    plt.ylabel('Number of Isoforms')
    plt.title('Number of Isoforms by Group')
    plt.legend(title='Group', loc='upper right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # 保存或显示图表
    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Process final_classification.txt files and plot results.")
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input file paths')
    parser.add_argument('-o', '--output', required=False, help='Output file path for the plot')
    parser.add_argument('-s', '--summary', required=True, help='Output file path for the summary statistics')
    args = parser.parse_args()

    # 处理文件
    data = process_files(args.input)
    
    # 计算转录本数量
    isoform_count = calculate_isoform_count(data)
    
    # 打印并保存统计结果
    print(isoform_count)
    isoform_count.to_csv(args.summary, index=False, sep='\t')
    
    # 绘制图表
    plot_summary(isoform_count, args.output)

if __name__ == '__main__':
    main()

