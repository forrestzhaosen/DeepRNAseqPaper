source("~/deep_RNAseq/code/util.R")
RNA_data_path <- "/storage/liu/RNAseq/freze2024Sep/qc_data/"
sample_list <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/gradiant/all_list.tsv")
sample_list_represent <- sample_list %>% filter(gradient==1,deduped=='yes',samples_per_run==1)
sample_list_represent$total_reads <- sapply(sample_list_represent$runID,get_total_reads)
color_palette <- c("#999999","#E69F00","#56B4E9","#009E73",
                   "#F0E442","#0072B2","#D55E00","#CC79A7")

count_example <- read_tsv(paste0(RNA_data_path,"/rnaseqc/09bdfe1a-4b1c-40fd-9a2f-60f7a00c03c3.gene_reads.gct"),skip = 2)
expression_matrix_all <- count_example[,1] %>%
  filter(!str_detect(Name, "PAR")) %>%
  mutate(Name = gsub("\\..*", "", Name)) %>%
  left_join(hgnc_esemble,by= c("Name"="ensembl_gene_id")) %>%
  filter(gene_biotype=="protein_coding",hgnc_id!="HGNC:30554") %>%
  filter(chromosome_name!="MT")
for (i in 1:nrow(sample_list_represent)) {
  runID <- sample_list_represent$runID[i]
  sampleID <- sample_list_represent$runID[i]
  expression <- read_tsv(paste0(RNA_data_path,"/rnaseqc/",runID,".gene_reads.gct"),skip=2) %>%
    filter(!str_detect(Name, "PAR")) %>%
    mutate(Name = gsub("\\..*", "", Name)) %>%
    left_join(hgnc_esemble,by= c("Name"="ensembl_gene_id")) %>%
    filter(gene_biotype=="protein_coding",hgnc_id!="HGNC:30554") %>%
    filter(chromosome_name!="MT")
  expression[[runID]][expression[[runID]]<50] <- 0
  expression_matrix_all <- cbind(expression_matrix_all,expression[,sampleID])
}

#calculate_MRSD <- function(median,Nfrag,targetcount){
#  res <- log(1-targetcount/(median*Nfrag),base=(Nfrag-1)/Nfrag)
#  return(res)
#}
calculate_MRSD <- function(median,targetcount){
  #res <- log10(targetcount/median)
  res <- targetcount/median/1000000
  return(res)
}

plot_list=list()
barplot_list=list()
#count_vec <- vector()
for (i in 1:length(unique(sample_list_represent$sample_type))){
  #i=1
  #geneID="ENSG00000187634"
  color=color_palette[i]
  sampletype=unique(sample_list_represent$sample_type)[i]
  reference_list <- sample_list_represent %>% filter(gradient==1,deduped=='yes',sample_type==sampletype)
  expression_matrix <- expression_matrix_all[,c("Name","hgnc_id","hgnc_symbol",reference_list$runID)]
  expression_matrix[,reference_list$runID] <- expression_matrix[,reference_list$runID]/reference_list$total_reads
  expression_matrix$median <-  apply(expression_matrix[, reference_list$runID], 1, median)
  combined_matrix <- tibble()
  for (count in c(1000,5000,10000,50000)){
    matrix <- expression_matrix
    matrix$targetcount <- count
    matrix$sample_type <- sampletype
    matrix$mrsd <- sapply(expression_matrix$median, calculate_MRSD, targetcount = count)
    #count_vec <- c(count_vec,length(matrix$mrsd[matrix$mrsd!="NaN" & mrsd!=Inf]))
    combined_matrix <- rbind(combined_matrix,matrix)
  }
  p1 <- ggplot(combined_matrix, aes(x=factor(targetcount), y=mrsd)) + # fill=name allow to automatically dedicate a color for each group
    geom_violin(fill=color)+
    #ylim(0,100000000)+
    xlab(paste0("Targeted count - ",sampletype))+
    scale_y_log10() +
    ylab(NULL)+
    theme_classic()+
    theme(legend.position = "none")
  plot_list[[i]]=p1
  number_matrix <- combined_matrix %>% group_by(targetcount) %>%
    dplyr::summarise(feasible_gene=sum(ifelse(mrsd != "NaN" & mrsd != Inf, 1, 0)))
  p2 <- ggplot(number_matrix, aes(x=factor(targetcount), y=feasible_gene)) +
    geom_bar(stat = "identity",fill=color)+
    xlab(paste0("Targeted count - ",sampletype))+
    ylab("No. of feasible genes")+
    theme_classic()+
    theme(legend.position = "none")
  barplot_list[[i]]=p2
  write_csv(combined_matrix,paste0("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/",sampletype,"-expression.csv"))
}
p1=do.call(grid.arrange, c(plot_list, ncol=1))
#p2=do.call(grid.arrange, c(barplot_list, ncol=2))
p1
#p2
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/MRSD_expression.pdf",plot =p1,
       width = 150,height = 150,units ="mm")


############## get expression plot for expression analysis
combined_matrix <- tibble()
#count_vec <- vector()
for (i in 1:length(unique(sample_list_represent$sample_type))){
  #i=1
  #geneID="ENSG00000187634"
  color=color_palette[i]
  sampletype=unique(sample_list_represent$sample_type)[i]
  reference_list <- sample_list_represent %>% filter(gradient==1,deduped=='yes',sample_type==sampletype)
  expression_matrix <- expression_matrix_all[,c("Name","hgnc_id","hgnc_symbol",reference_list$runID)]
  expression_matrix$median <-  log10(apply(expression_matrix[, reference_list$runID], 1, median))
  expression_matrix$sample_type <- sampletype
  expression_matrix <- expression_matrix[,c("Name","hgnc_id","hgnc_symbol","median","sample_type")] %>%
    filter(median != 0)
  combined_matrix <- rbind(combined_matrix,expression_matrix)
}

p1 <- ggplot(combined_matrix, aes(x=sample_type, y=median)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin(aes(fill=sample_type))+
  #ylim(0,10000)+
  xlab("sample type")+
  ylab("log(gene count)")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 20))
p1
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/expression_big_ass.pdf",plot =p1,
       width = 200,height = 150,units ="mm")
