import pandas as pd
import os
configfile: "config.yaml"


#SAMPLES = ["Cs1_rep1","Cs1_rep2","Cs2_rep1","Cs2_rep2","Cs3_rep1","Cs3_rep2","Cs4_rep1","Cs4_rep2","Cs5_rep1","Cs5_rep2"]
#SAMPLES = config["samples"]
#读取文件并获得第二列样本名称
#with open('samples.txt','r') as file:
  #逐行读取内容
# lines = file.readlines()
  #使用列表解析来从每行中提取样本名称
# SAMPLES = [line.split('\t').strip() for line in lines]
SAMPLES = config["samples"]
libType = config["libType"]


rule all:
   input:
    #用列表推导式生成
    #["clean/{}_{}.fastp.fastq.gz".format(sample,i) for sample in SAMPLES for i in [1,2]]
    #expand("clean/{}_{}fastp.fastq.gz", sample = SAMPLES, i = [1,2]),
    expand("aligned/{sample}.bam" if config["aligner"] == "hisat2" else "aligned/{sample}Aligned.out.bam", sample=SAMPLES),
    expand("quanti/{sample}.count",sample =SAMPLES),
    counts = "quanti/genes.counts.matrix",
    tpm = "quanti/genes.TMM.TPM.matrix",
    de_result = "quanti/de_result.tsv"
rule trim_fastq:
  input:
     r1 = lambda wildcards: "{}/{}.R1.fastq.gz".format(config["data_directory"], wildcards.sample),
     r2 = lambda wildcards: "{}/{}.R2.fastq.gz".format(config["data_directory"], wildcards.sample)
  output:
    r1 = "clean/{sample}_1.fastp.fastq.gz",
    r2 = "clean/{sample}_2.fastp.fastq.gz",
    html = "clean/{sample}.fastp.html",
    json = "clean/{sample}.json"
  threads: 15
  log: "clean/{sample}.fastp.log"
  shell:
    """
    fastp -w {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json} 1>{log} 2>&1

    """
if config["aligner"] == "hisat2":
#构建参考基因组
 rule build_hisat2_index:
  input:
    genome = config["genome"]
  output:
    expand("genome.hisat2.{i}.ht2",i = range(1,9))
  shell:
    """
       hisat2-build {input.genome} genome.hisat2  #输出文件写死，流程简单

    """

# hisat2 比对
 rule align_hisat2:
   input:
     r1 = "clean/{sample}_1.fastp.fastq.gz",
     r2 = "clean/{sample}_2.fastp.fastq.gz",
     idx = expand("genome.hisat2.{i}.ht2",i = range(1,9))
   output:
     sam = "aligned/{sample}.sam",
     log = "aligned/{sample}.hisat2.log"
   params:
     libType="RF"
   threads: 35
   shell:
    """
      hisat2 -x genome.hisat2  -p {threads} -1 {input.r1} -2 {input.r2} --new-summary --rna-strandness {params.libType}  -S {output.sam} 2>{output.log}
    """
#构建STAR
# 如果配置文件选择的是STAR，就加载这段代码
# STAR index building
lambda wildcards: "{}/{}_1.fastq.gz".format(config["data_directory"], wildcards.samples),
lambda wildcards: "{}/{}_2.fastq.gz".format(config["data_directory"], wildcards.samples),
rule build_star_index:
    input:
        fasta=config["genome"]
    output:
        complete=directory("star_index/")
    threads: 35
    run:
        if config["aligner"] == "star":
            shell("""
                ulimit -n 10000
                STAR --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeDir star_index/ \
                --genomeFastaFiles {input.fasta};\
                touch {output.complete}
            """)

# STAR align
rule align_star:
    input:
        r1 = "clean/{sample}_1.fastp.fastq.gz",
        r2 = "clean/{sample}_2.fastp.fastq.gz",
        idx = directory("star_index/")
    output:
        bam = "aligned/{sample}Aligned.out.bam"
    params:
        prefix = "aligned/{sample}"
    threads: 20
    shell:
        """
        if [ "{config[aligner]}" == "star" ]; then
            STAR --genomeDir {input.idx} \
                --readFilesIn {input.r1} {input.r2} \
                --runThreadN {threads} \
                --outSAMtype BAM Unsorted \
                --readFilesCommand zcat \
                --outFileNamePrefix {params.prefix}
                 echo "Finished STAR alignment for {wildcards.sample}"
        else
            echo "Aligner not set to STAR in the config. Exiting..."
            exit 1
        fi
        """

        

# sam转bam
if config["aligner"] == "hisat2":
 rule sam_to_bam:
    input:
        sam = "aligned/{sample}.sam" 
    output:
        bam = "aligned/{sample}.bam",
        idx = "aligned/{sample}.bam.bai"
    threads: 40
    run:
        if config["aligner"] == "hisat2":
            shell("""
                samtools sort -@ {threads} -o {output.bam} {input.sam}
                samtools index  {output.bam}
                rm {input.sam}
            """)

# featureCounts表达定量
def hisat2_to_featureCounts(strandness):
    if strandness == "FR":
        return 1
    elif strandness == "RF":
        return 2
    elif strandness == "unstranded":
        return 0
    else:
        raise ValueError("Invalid strandness value. Only 'FR' or ’RF‘ are allowed for HISAT2.")
#if config["counts"] == "featruecounts"
#使用featruecounts进行定量
rule quanti_featruecounts:
    input:
        bam = "aligned/{sample}.bam" if config["aligner"]=="hisat2" else "aligned/{sample}Aligned.out.bam",
        gtf = config["gtf"]
    output:
        counts = "quanti/{sample}.count"
    params:
        out_perfix = "quanti/{sample}",
        libType2 = hisat2_to_featureCounts(libType),
        featurecounts = config["featurecounts"]
    shell:
        """
        Rscript {params.featurecounts} \
            -b {input.bam} \
            -g {input.gtf} \
            -s {params.libType2} \
            -o {params.out_perfix}
        """
#使用htseq进行定量
#if config["counts"] == "htseq"
# rule count_htseq:
# input:
#        bam="aligned/{sample}.bam" if config["aligner"] == "hisat2" else "aligned/{sample}Aligned.out.bam", # 或者你的STAR的输出文件
#        gtf="path/to/your/annotation.gtf"
#    output:
#        counts="counts/{sample}.htseq.txt"
#    params:
#        yes_or_no = config["yes_or_no"]
#    shell:
#        """
#        htseq-count -f bam -s {params.yes_or_no} -t exon -i gene_id {input.bam} {input.gtf} > {output.counts}
#        """
#注意：

#你需要按照你的需求更改文件路径和规则参数。
#在HTSeq-count命令中，我使用了以下参数：
#-f bam: 输入是BAM格式。
#-s no: 假设没有链特异性。如果你的数据是链特异的，请使用yes。
#-t exon: 计数exons。
#-i gene_id: 用于标识基因的属性。
#请根据你的实际情况修改这些参数。如果你使用了其他的比对工具或有不同的文件结构，也请确保更改路径和文件名以匹配你的工作流。

#最后，确保你有一个配置的Snakemake执行环境，并在包含这个Snakefile的目录中执行它。


#合并表达矩阵
rule abundance_estimates_to_matrix:  
  input:
    counts_files = expand("quanti/{sample}.count",sample = SAMPLES)
  output:
    counts = "quanti/genes.counts.matrix",
    tpm = "quanti/genes.TMM.TPM.matrix"
  shell:
    """
    ls {input.counts_files} >quanti_files.txt
    perl ../software/RunFeatureCounts/abundance_estimates_to_matrix.pl \
      --est_method featureCounts \
      --out_prefix quanti/genes \
      --quant_files quanti_files.txt
    
    sed -i '1s/^/gene_id/' quanti/genes.counts.matrix
    sed -i '1s/^/gene_id/' quanti/genes.TMM.TPM.matrix
    """

#差异表达分析
rule DE_analysis:
  input:
    counts = "quanti/genes.counts.matrix",
  output:
    de_result = "quanti/de_result.tsv"
  params:
    TRINITY_HOME=config["TRINITY_HOME"],
    de_method = config["de_method"],
    samples_file=config["samples_file"],
    contrasts_file=config["contrasts_file"]

  shell:
   """
   perl {params.TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl\
         --matrix {input.counts} \
         --method {params.de_method} \
         --samples_file {params.samples_file} \
         --output DE_analysis \
         --contrasts {params.contrasts_file}
   
   awk 'FNR==1 && NR!=1{{next}}{{print}}' DE_analysis/*DE_results > {output.de_result}
   """

