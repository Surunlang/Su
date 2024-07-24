import pandas as pd
import argparse

def main(input_file):
    # 读取classification.txt文件
    df = pd.read_csv(input_file, sep='\t')

    # 显示文件的前几行和列名称
    print(df.head())
    print(df.columns)

    # 检查是否存在FSM类型的同工型
    if 'structural_category' in df.columns:
        total_fsm = df[df['structural_category'] == 'full-splice_match'].shape[0]
        print(f'Total FSM isoforms: {total_fsm}')

        if total_fsm > 0:
            # 统计未匹配注释TSS的FSM Isoforms数量
            fsm_tss_not_matching = df[(df['structural_category'] == 'full-splice_match') & (df['diff_to_gene_TSS'] != 0)].shape[0]

            # 统计未匹配注释TTS的FSM Isoforms数量
            fsm_tts_not_matching = df[(df['structural_category'] == 'full-splice_match') & (df['diff_to_gene_TTS'] != 0)].shape[0]

            # 计算比例
            tss_not_matching_percentage = (fsm_tss_not_matching / total_fsm) * 100
            tts_not_matching_percentage = (fsm_tts_not_matching / total_fsm) * 100

            # 输出结果
            print(f'未匹配注释TSS的FSM Isoforms比例：{tss_not_matching_percentage:.2f}%')
            print(f'未匹配注释TTS的FSM Isoforms比例：{tts_not_matching_percentage:.2f}%')
        else:
            print("No FSM isoforms found in the classification.txt file.")
    else:
        print("The file does not contain a 'structural_category' column.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process SQANTI3 classification file.')
    parser.add_argument('-i', '--input', required=True, help='Input classification file')
    args = parser.parse_args()
    main(args.input)

