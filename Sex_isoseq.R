# 加载所需的库
library(optparse)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Path to the input expression matrix file 输入表达矩阵文件的路径", metavar = "character"),
  make_option(c("-c", "--contract"), type = "character", default = NULL, help = "Path to the contract file 输入分组对比文件的路径", metavar = "character"),
  make_option(c("-s", "--sample"), type = "character", default = NULL, help = "Path to the sample file 输入样本文件的路径", metavar = "character"),
  make_option(c("-n", "--num_replicates"), type = "integer", default = 3, help = "Number of replicates 每个组的重复次数 [default= %default]", metavar = "integer"),
  make_option(c("-G", "--group_prefix_G"), type = "character", default = "G", help = "Prefix for group G 组 G 的前缀 [default= %default]", metavar = "character"),
  make_option(c("-X", "--group_prefix_X"), type = "character", default = "X", help = "Prefix for group X 组 X 的前缀 [default= %default]", metavar = "character"),
  make_option(c("-t", "--threshold"), type = "integer", default = 10, help = "Threshold for expression filtering 表达过滤的阈值 [default= %default]", metavar = "integer"),
  make_option(c("-f", "--final_classification"),
              type = "character",
              default = NULL,
              help = "Path to final_classification file",
              metavar = "FILE")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查必需的参数是否提供
if (is.null(opt$input) | is.null(opt$contract) | is.null(opt$sample)) {
  print_help(opt_parser)
  stop("Input file, contract file and sample file must be provided. 必须提供输入文件、分组文件和样本文件。", call. = FALSE)
}

# 定义处理数据的函数
process_data <- function(input_file, contract_file, sample_file, num_replicates, threshold = 10, group_prefix_G = "G", group_prefix_X = "X") {
  # 读取表达矩阵
  x <- read.csv(input_file, sep = "\t", row.names = 1)
  print("Expression matrix column names:")
  print(colnames(x))
  
  # 将表达矩阵转换为数值型
  x <- data.frame(lapply(x, as.numeric), row.names = row.names(x))
  print(nrow(x))
  # 读取分组文件
  contract <- read.csv(contract_file, sep = "\t", header = FALSE)
  colnames(contract) <- c("Group_G", "Group_X")
  print("Contract file:")
  print(contract)
  
  # 读取样本文件并跳过空行
  samples <- read.csv(sample_file, sep = "\t", header = FALSE, skipNul = TRUE, na.strings = "")
  colnames(samples) <- c("Group", "Sample")
  samples <- na.omit(samples)
  print("Sample file:")
  print(samples)
  
  
  # 匹配样本名并重命名列名
  matched_samples <- samples$Sample[samples$Sample %in% colnames(x)]
  unmatched_samples <- samples$Sample[!samples$Sample %in% colnames(x)]
  if (length(unmatched_samples) > 0) {
    cat("Unmatched samples found: ", paste(unmatched_samples, collapse = ", "), "\n")
    stop("Please check the sample names in your sample file and expression matrix.")
  }
  print("Filtered expression matrix:")
  print(head(x))
  
  
  # 过滤掉表达量总和小于阈值的行
  row_sums <- rowSums(x)
  print(colnames(x))
  print(nrow(x))
  print("Hello")
  x <- x[row_sums >= 10, ]
  print(dim(x))
  print(nrow(x))
  # 自动匹配列名并进行分组计算
  calculate_group_mean <- function(data, sample_names, num_replicates) {
    matched_samples <- sample_names[sample_names %in% colnames(data)]
    print(paste("Matched samples for group:", paste(matched_samples, collapse = ", ")))
    if (length(matched_samples) == 0) {
      stop("No valid samples found in data for group: ", paste(sample_names, collapse = ", "))
    }
    split_samples <- split(matched_samples, ceiling(seq_along(matched_samples) / num_replicates))
    means <- sapply(split_samples, function(samples) rowMeans(data[, samples, drop = FALSE]))
    return(means)
  }
  
  # 计算分组的平均值
 #group_means <- lapply(unique(samples$Group), function(group) {
    #sample_names <- samples$Sample[samples$Group == group]
    #valid_samples <- sample_names[sample_names %in% colnames(x)]
    #if (length(valid_samples) == 0) {
      #stop("No valid samples found in data for group: ", group)
   # }
   # print(paste("Processing group:", group))
   # print(sample_names)
   # calculate_group_mean(x, sample_names, num_replicates)
  #})
   group_means <- lapply(unique(samples$Group), function(group) {
    sample_names <- samples$Sample[samples$Group == group]
    valid_samples <- sample_names[sample_names %in% colnames(x)]
    if (length(valid_samples) == 0) {
      stop("No valid samples found in data for group: ", group)
    }
    print(paste("分组归类是:", group))
    print(sample_names)
    calculate_group_mean(x, sample_names, num_replicates)
  })
  # 合并结果
  group_means_df <- do.call(cbind, group_means)
  colnames(group_means_df) <- unique(samples$Group)
  row.names(group_means_df) <- row.names(x)
  print(head(group_means_df))
  
  # 初始化结果
  G_single_expression_counts <- numeric(nrow(contract))
  X_single_expression_counts <- numeric(nrow(contract))
  common_expression_counts <- numeric(nrow(contract))
  group_pairs <- paste0(contract$Group_G, "_vs_", contract$Group_X)
  print(group_pairs)
  # 存储表达ID的列表
  G_single_expression_ids <- list()
  X_single_expression_ids <- list()
  common_expression_ids <- list()
  
  # 比较组别
  for (i in 1:nrow(contract)) {
    G_group <- group_means_df[, contract$Group_G[i]]
    X_group <- group_means_df[, contract$Group_X[i]]
   print(paste("Comparing[比较]", contract$Group_G[i], "with", contract$Group_X[i])) 
    # 应用规则
    G_group[G_group < threshold] <- 0
    X_group[X_group < threshold] <- 0
    
    G_single_expression <- (G_group != 0 & X_group == 0) | (G_group > 400 & X_group < 1)
    X_single_expression <- (X_group != 0 & G_group == 0) | (X_group > 400 & G_group < 1)
    common_expression <- (G_group != 0 & X_group != 0) & !((G_group > 400 & X_group < 1) | (X_group > 400 & G_group < 1))
    
    G_single_expression_counts[i] <- sum(G_single_expression)
    X_single_expression_counts[i] <- sum(X_single_expression)
    common_expression_counts[i] <- sum(common_expression)
    
    G_single_expression_ids[[i]] <- data.frame(ID = row.names(group_means_df)[G_single_expression], Group = group_pairs[i])
    X_single_expression_ids[[i]] <- data.frame(ID = row.names(group_means_df)[X_single_expression], Group = group_pairs[i])
    common_expression_ids[[i]] <- data.frame(ID = row.names(group_means_df)[common_expression], Group = group_pairs[i])
  }
  
  # 合并所有组的结果并写入一个文件
  all_G_single_expression_ids <- do.call(rbind, G_single_expression_ids)
  all_X_single_expression_ids <- do.call(rbind, X_single_expression_ids)
  all_common_expression_ids <- do.call(rbind, common_expression_ids)
  
  write.csv(all_G_single_expression_ids, file = "all_G_single_expression_ids.csv", row.names = FALSE)
  write.csv(all_X_single_expression_ids, file = "all_X_single_expression_ids.csv", row.names = FALSE)
  write.csv(all_common_expression_ids, file = "all_common_expression_ids.csv", row.names = FALSE)
  
  # 生成统计结果
  results <- data.frame(
    Group_Pair = group_pairs,
    G_Single_Expression = G_single_expression_counts,
    X_Single_Expression = X_single_expression_counts,
    Common_Expression = common_expression_counts
  )
  print(results)
  
  # 作图
  results_long <- melt(results, id.vars = "Group_Pair", variable.name = "Expression_Type", value.name = "Count")
  
  # 自定义配色
  colors <- c("G_Single_Expression" = "#FF9999", "X_Single_Expression" = "#99CCFF", "Common_Expression" = "#66CC66")
  
  # 保存图片为PDF
  p <- ggplot(results_long, aes(x = Group_Pair, y = Count, fill = Expression_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Isoform Expression Comparison", x = "Group Pair", y = "Count") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.line = element_line(color = "black"), 
          axis.ticks = element_line(color = "black")) +
    geom_text(aes(label = Count), vjust = -0.3, position = position_dodge(0.9), size = 3)
  ggsave("Isoform Expression Comparison.pdf", plot = p, device = "pdf", width = 10, height = 6)
  
  # 分成两组进行同样的分析，只分G和X
  G_group_means <- rowMeans(group_means_df[, grepl(group_prefix_G, colnames(group_means_df))])
  X_group_means <- rowMeans(group_means_df[, grepl(group_prefix_X, colnames(group_means_df))])
  
  # 计算G组和X组的单独表达
  G_single_expression_total <- sum((G_group_means != 0 & X_group_means == 0) | (G_group_means > 400 & X_group_means < 10))
  X_single_expression_total <- sum((X_group_means != 0 & G_group_means == 0) | (X_group_means > 400 & G_group_means < 10))
  common_expression_total <- sum((G_group_means != 0 & X_group_means != 0) & !((G_group_means > 400 & X_group_means < 10) | (X_group_means > 400 & G_group_means < 10)))
  
  # 生成总统计结果
  results_total <- data.frame(
    Expression_Type = c("G_Single_Expression", "X_Single_Expression", "Common_Expression"),
    Count = c(G_single_expression_total, X_single_expression_total, common_expression_total)
  )
  
  print(results_total)
  
  # 保存总的表达量图像为PDF
  p1 <- ggplot(results_total, aes(x = Expression_Type, y = Count, fill = Expression_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.3) +
    labs(title = "Total Isoform Expression Comparison", x = "Expression Type", y = "Count") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.line = element_line(color = "black"), 
          axis.ticks = element_line(color = "black")) +
    geom_text(aes(label = Count), vjust = -0.3, position = position_dodge(0.5), size = 3)
  ggsave("Isoform Expression Comparison by Group_single_Common.pdf", plot = p1, device = "pdf", width = 10, height = 6)
  # 额外的柱状图，仅显示G和X的单独表达
  results_total_by_group <- data.frame(
    Group = c("G_Group", "X_Group"),
    G_Single_Expression = c(G_single_expression_total, NA),
    X_Single_Expression = c(NA, X_single_expression_total)
  )
  
  results_total_by_group_long <- melt(results_total_by_group, id.vars = "Group", variable.name = "Expression_Type", value.name = "Count")
  results_total_by_group_long <- results_total_by_group_long[results_total_by_group_long$Count != 0, ]
  results_total_by_group_long <- na.omit(results_total_by_group_long)
  
  # 保存G和X单独表达的图像为PDF
  p2 <- ggplot(results_total_by_group_long, aes(x = Group, y = Count, fill = Expression_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.3) +
    labs(title = "Isoform Expression Comparison by Group", x = "Group", y = "Count") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.line = element_line(color = "black"), 
          axis.ticks = element_line(color = "black")) +
    geom_text(aes(label = Count), vjust = -0.3, position = position_dodge(0.5), size = 3)
    
   ggsave("Isoform Expression Comparison by Two_Group_Single.pdf", plot = p2, device = "pdf", width = 10, height = 6)
  #真麻烦，不想写了a
  # 保存总的单独表达的ID到文件
  G_single_expression_ids_flat <- do.call(rbind, G_single_expression_ids)
  X_single_expression_ids_flat <- do.call(rbind, X_single_expression_ids)
  common_expression_ids_flat <- do.call(rbind, common_expression_ids)
  write.csv(G_single_expression_ids_flat, "G_single_expression_ids_total.csv", row.names = FALSE)
  write.csv(X_single_expression_ids_flat, "X_single_expression_ids_total.csv", row.names = FALSE)
  write.csv(common_expression_ids_flat, "common_expression_ids_total.csv", row.names = FALSE)
}
# 运行处理数据的函数
process_data(opt$input, opt$contract, opt$sample, opt$num_replicates, opt$threshold, opt$group_prefix_G, opt$group_prefix_X)
if (!is.null(opt$final_classification)) {
  # 确保对象 G_single_expression_ids_flat 和 X_single_expression_ids_flat 已生成
  if (file.exists("G_single_expression_ids_total.csv") && file.exists("X_single_expression_ids_total.csv")) {
    # 读取生成的文件
    G_single_expression_ids_flat <- read.csv("G_single_expression_ids_total.csv")
    X_single_expression_ids_flat <- read.csv("X_single_expression_ids_total.csv")
    
    # 读取classification文件
    classification <- read.delim(opt$final_classification, sep = "\t")
    classification$associated_gene <- gsub("gene-", "", classification$associated_gene)
    
    # 合并工蜂的独有isoform的id，转换gene_id
    G <- left_join(G_single_expression_ids_flat, classification, by = c("ID" = "isoform"))
    filtered_G <- G[grepl("^novel", G$associated_gene), ]
    G_novel_gene <- table(filtered_G$Group)
    G_novel_gene <- as.data.frame(G_novel_gene)
    G_novel_gene$Type <- "Worker"
    colnames(G_novel_gene) <- c("Group","Count","Type")
    
    # 合并雄蜂的独有isoform的id，转换gene_id
    X <- left_join(X_single_expression_ids_flat, classification, by = c("ID" = "isoform"))
    filtered_X <- X[grepl("^novel", X$associated_gene), ]
    X_novel_gene <- table(filtered_X$Group)
    X_novel_gene <- as.data.frame(X_novel_gene)
    X_novel_gene$Type <- "Drone"
    colnames(X_novel_gene) <- c("Group","Count","Type")
    
    # 合并矩阵
    singel_gene_or_AS <- rbind(filtered_G, filtered_X)
    write.csv(singel_gene_or_AS, file = "singel_gene_or_AS.cvs")
    combined_data <- rbind(G_novel_gene, X_novel_gene)
    
    # 绘制分体条形图
    colors <- c("Worker" = "#FF9999", "Drone" = "#99CCFF")
    p3 <- ggplot(combined_data, aes(x = Group, y = Count, fill = Type)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = Count), vjust = -0.3, position = position_dodge(0.9), size = 3) +
      labs(title = "Novel Genes in Worker and Drone Bees", 
           x = "Group", 
           y = "Count of Novel Genes and AS") +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(panel.grid = element_blank(), 
            axis.line = element_line(color = "black"), 
            axis.ticks = element_line(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("Novel_Gene_and_AS.pdf", plot = p3, device = "pdf", width = 10, height = 6)
  } else {
    stop("G_single_expression_ids_total.csv or X_single_expression_ids_total.csv not found")
  }
}
print("Your Work Done!")

