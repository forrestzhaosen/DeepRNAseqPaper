source("~/deep_RNAseq/code/util.R")
library(ggplot2)
library(corrplot)
RNA_data_path <- "/storage/liu/home/u247123/deep_RNAseq/freeze2023Oct/qc_data/"

sample_list <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/gradiant/Gradient_list.tsv")
sample_list_represent <- sample_list %>% filter(gradient==1,deduped=='yes',samples_per_run==1,sample_type!='Blood')
sample_list_represent$total_reads <- sapply(sample_list_represent$runID,get_total_reads)
expression_matrix <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/count_matrix.tsv") 
expression_matrix <- expression_matrix[,c("Name","hgnc_id","hgnc_symbol",sample_list_represent$runID)]
colnames(expression_matrix)[match(sample_list_represent$runID, colnames(expression_matrix))] <- sample_list_represent$paper_ID

cor_matrix <- cor(expression_matrix[4:12], method = "pearson")


pdf("corrplot.pdf", width = 8, height = 6)
corp <- corrplot(cor_matrix, method = "color")
dev.off()
