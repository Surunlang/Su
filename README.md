# Su
tts_tss_count.py这是用来统计classification文件的TSS和TTS的FSM情况，一般是过滤掉了10000以下序列，想降低标准可以自行调整
#
Q20_Q30_caculate.py这是用来批量计算fastp软件留下的信息，统计测序文件的Q20和Q30
#
hisat2_map_ratio.py 这是用来批量统计hisat2和测序数据的对比率
#
stat_hisat2bam_files.sh 这是用来统计Sample  Total   Unmapped (%)    Uniquely mapped (%)     Multiple mapped（%）Total mapped(%) 这几个指标的文件。

#orgdb.py和create_orgdb.R是用来生成非模式物种的Go包，只要通过emmaper.py软件注释成annotation文件后，将annotation文件
和该物种的蛋白文件传递至orgdb.py里即可生成自己物种的orgdb包，记得蛋白文件的序列号要换成gene_id

#Snakefile、config.yml、sample.txt、sample_info.txt、contracts.txt这四个加起来，在需要测序文件、参考基因组的gtf和fasta文件即可全自动运行Snakemake的RNA-seq分析流程
