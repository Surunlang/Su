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

#singularity exec ~/workspace/software/rnasamba_latest.sif rnasamba classify -p final.faa rna_classification.tsv final.fasta 
~/workspace/software/full_length_weights.hdf5
这是我自己的预测非编码RNA和编码RNA的目录，用这套命令传输至coding_novel.py
#final_classification.txt这个文件是SQANTI3运行的最终结果，可以去试试，这套组合针对的是iso-seq测序
#coding_novel.py是用来承接rnasamba软件预测的结果，-t或者--tsv 输入rnasamb软件输出的tsv文件（rna_classification.tsv）即可，-i 再输入final_classification.txt文件即可，再输入-o my_species_dir(自定义文件夹)

