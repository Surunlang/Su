#!/usr/bin/env Rscript

#加载必要的库
library(readr)
library(scales)
library(dplyr)
library(tidyr)
library(AnnotationForge)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(BiocGenerics)
library(seqinr)  # 读取 FASTA 文件所需

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script_name.R <cleaned_data_path> <output_dir> <proteins_file>")
}
cleaned_data_path <- args[1]  # CSV 文件的完整路径
output_dir <- args[2]         # 输出目录的完整路径
proteins_file <- args[3]      # 蛋白质 FASTA 文件的路径

# 确认文件存在
if (!file.exists(cleaned_data_path)) {
    stop("文件不存在：", cleaned_data_path)
}
if (!file.exists(proteins_file)) {
    stop("蛋白质文件不存在：", proteins_file)
}

script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

# 读取数据
data <- read.csv(cleaned_data_path, stringsAsFactors = FALSE)

# 处理 GO 列，分拆成多行，并添加 EVIDENCE 列
go_info <- data %>%
  select(GID, GO) %>%
  filter(GO != "-") %>%
  separate_rows(GO, sep = ",", convert = FALSE) %>%
  filter(GO != "") %>%
  mutate(EVIDENCE = "IEA")

# 准备其他需要的数据表
gene_info <- data %>%
  select(GID, Gene_Name) %>%
  filter(Gene_Name != "-")
gene2go <- go_info
eggnog_anno <- length(gene_info$GID)
emapper <- data
# 创建 OrgDB
makeOrgPackage(
  gene_info = gene_info,
  go = go_info,
  maintainer = 'zhangsan <zhangsan@genek.cn>',
  author = 'zhangsan',
  outputDir = output_dir,
  tax_id = 0000,
  genus = 'M',
  species = 'y',
  goTable = "go",
  version = "1.0"
)

# 创建输出目录并设置安装路径
dir.create(file.path(output_dir, "R_Library"), recursive = TRUE, showWarnings = FALSE)
pkgbuild::build(file.path(output_dir, "org.My.eg.db"), dest_path = output_dir)
install.packages(file.path(output_dir, "org.My.eg.db_1.0.tar.gz"), repos = NULL, type = "source", lib = file.path(output_dir, "R_Library"))
library(org.My.eg.db, lib = file.path(output_dir, "R_Library"))

###############################################
# GO statistics and plot
###############################################
all_gene <- getName.list(read.fasta(file = proteins_file, seqtype = 'AA'))
go_anno = length(unique(gene2go$GID))
go_bp <- clusterProfiler::groupGO(gene     = all_gene,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "BP",
                 level    = 2,
                 readable = FALSE)

go_bp <- as.data.frame(go_bp)
go_bp$GO_Class <- "Biological Process"

go_cc <- clusterProfiler::groupGO(gene     = all_gene,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "CC",
                 level    = 2,
                 readable = FALSE)

go_cc <- as.data.frame(go_cc)
go_cc$GO_Class <- "Cellular Component"

go_mf <- clusterProfiler::groupGO(gene     = all_gene,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "MF",
                 level    = 2,
                 readable = FALSE)
go_mf <- as.data.frame(go_mf)
go_mf$GO_Class <- "Molecular Function"

go_all <- rbind(go_bp, go_cc, go_mf)
p <- ggplot(go_all) +
  geom_bar(aes(x = Description,
               y = Count,
               fill = GO_Class),
           stat = "identity") +
  facet_wrap(~GO_Class, scales = "free_x") +
  labs(title = "GO function classification", y = "Number of genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.title.x = element_blank(),
        legend.position = "none")

write.table(go_all, file.path(output_dir, "go.txt"), sep = "\t", quote = F)
ggsave(file.path(output_dir, "go.pdf"), p, width = 20, height = 7)

###############################################
# Pathway statistics and plot 
###############################################
gene2pathway <- dplyr::select(data, GID, Pathway) %>%
  separate_rows(Pathway, sep = ',', convert = FALSE) %>%
  filter(str_detect(Pathway, 'ko'))

# Assuming you have already loaded kegg_info.RData
load(file = paste(script_dir, "kegg_info.RData", sep = "/"))

gene2pathway <- gene2pathway %>%
  left_join(pathway2name, by = "Pathway") %>%
  select(GID, Pathway, Pathway_Name, Pathway_Class, Pathway_Subclass) %>%
  distinct() %>%
  na.omit()
pathway_anno <- length(unique(gene2pathway$GID))
pathway_stat <- dplyr::select(gene2pathway, GID, Pathway_Class, Pathway_Subclass) %>% 
  distinct() %>% 
  group_by(Pathway_Class, Pathway_Subclass) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = (Count / sum(Count)) * 100)  # 直接计算百分比

p <- ggplot(pathway_stat, aes(x = Pathway_Subclass, y = Percentage)) +
  geom_bar(aes(fill = Pathway_Class), stat = 'identity') +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), nudge_y = 0.5) +  # 显示百分比
  labs(y = "Percent of genes (%)", x = "", fill = "Class") +
  coord_flip() +
  theme_classic()

ggsave(file.path(output_dir, "pathway.pdf"), p, width = 20, height = 7)
write.table(pathway_stat, file.path(output_dir, "pathway_stat.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

###############################################
# COG statistics and plot
###############################################
# 确保读取 COG 分类数据
cog_funclass <- read_delim(paste(script_dir, "cog_funclass.tab", sep = "/"),
                           "\t", escape_double = FALSE, trim_ws = TRUE)

# 函数用于插入逗号，准备数据
insert_comma <- function(x) {
  str_c(x, sep = '', collapse = ',')
}

# 选择和准备 COG 数据
gene2cog <- dplyr::select(data, GID, COG) %>%
  filter(!is.na(COG) & COG != '-') %>%
  mutate(COG = sapply(str_split(COG, ''), insert_comma)) %>%
  separate_rows(COG, sep = ',', convert = FALSE) %>%
  left_join(cog_funclass, by = 'COG')
cog_anno <- length(unique(gene2cog$GID))
# 统计每个 COG 类别的基因数量
cog_stat <- gene2cog %>%
  group_by(COG_Name) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

# 绘制 COG 分类的图表
p <- ggplot(cog_stat, aes(x = COG_Name, y = Count, fill = COG_Name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.5) +
  labs(title = "COG/KOG Function Classification", x = "", y = "Number of genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.key.size = unit(1, "line"),
        legend.text = element_text(size = 7.5)) +
  guides(fill = guide_legend(ncol = 1))

# 保存图表和数据
ggsave(file.path(output_dir, "cog.pdf"), p, width = 16, height = 7)
write.table(gene2cog, file.path(output_dir, "cog.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

###############################################
# Number and percentage
###############################################
total_gene = length(all_gene)
anno_stat <- tibble(
  database = c("EggNOG", "GO", "COG/KOG", "KEGG Pathway"),
  number = c(eggnog_anno, go_anno, cog_anno, pathway_anno),
  percentage = (c(eggnog_anno, go_anno, cog_anno, pathway_anno) / total_gene) * 100
)

write.table(anno_stat, file.path(output_dir, "anno_stat.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

