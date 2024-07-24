#!/bin/bash

# 使用说明函数
usage() {
    echo "Usage: $0 -i <input_directory> -t <threads>"
    exit 1
}

# 解析命令行参数
while getopts "i:t:" opt; do
    case ${opt} in
        i )
            INPUT_DIR=$OPTARG
            ;;
        t )
            THREADS=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# 检查是否提供了输入路径和线程数
if [ -z "$INPUT_DIR" ] || [ -z "$THREADS" ]; then
    usage
fi

# 创建一个临时文件来保存统计信息
STATS_FILE="bam_stats.txt"
echo -e "Sample\tTotal\tUnmapped (%)\tUniquely mapped (%)\tMultiple mapped (%)\tTotal mapped (%)" > $STATS_FILE

# 统计 BAM 文件总数
bam_files=("$INPUT_DIR"/*.bam)
total_files=${#bam_files[@]}
current_file=1

# 遍历目录中的每个 BAM 文件
for bam_file in "${bam_files[@]}"; do
    # 获取文件名（不包含路径和扩展名）
    base_name=$(basename "$bam_file" .bam)
    
    echo "Processing file $current_file of $total_files: $bam_file ..."
    
    # 生成对齐统计信息
    flagstat_output=$(samtools flagstat -@ "$THREADS" "$bam_file")
    
    # 解析统计信息
    total=$(echo "$flagstat_output" | awk 'NR==1 {print $1}')
    mapped=$(echo "$flagstat_output" | awk 'NR==5 {print $1}')
    unmapped=$(echo "$total - $mapped" | bc)
    properly_paired=$(echo "$flagstat_output" | awk 'NR==9 {print $1}')
    singletons=$(echo "$flagstat_output" | awk 'NR==11 {print $1}')
    
    unmapped_percent=$(echo "scale=2; ($unmapped / $total) * 100" | bc)
    uniquely_mapped_percent=$(echo "scale=2; ($properly_paired / $total) * 100" | bc)
    multiple_mapped_percent=$(echo "scale=2; ($singletons / $total) * 100" | bc)
    total_mapped_percent=$(echo "scale=2; ($mapped / $total) * 100" | bc)
    
    # 将统计信息添加到表格
    echo -e "$base_name\t$total\t$unmapped_percent\t$uniquely_mapped_percent\t$multiple_mapped_percent\t$total_mapped_percent" >> $STATS_FILE
    
    echo "Finished processing file $current_file of $total_files: $bam_file"
    current_file=$((current_file + 1))
done

echo "All BAM files processed. Statistics saved to $STATS_FILE"

