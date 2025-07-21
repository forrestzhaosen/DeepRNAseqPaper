source("~/deep_RNAseq/code/util.R")
RNA_data_path <- "/storage/liu/RNAseq/freze2024Sep/qc_data/"
sample_list <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/gradiant/all_list.tsv")
sample_list_represent <- sample_list %>% filter(gradient==1,deduped=='yes',samples_per_run==1)
sample_list_represent$total_reads <- sapply(sample_list_represent$runID,get_total_reads)
color_palette <- c("#999999","#E69F00","#56B4E9","#009E73",
                   "#F0E442","#0072B2","#D55E00","#CC79A7")
## read splicing matrix
#junction_matrix_all <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/junction_matrix.tsv")
junction_df <- tibble()
for (i in 1:nrow(sample_list)) {
  runID <- sample_list$runID[i]
  sampleID <- sample_list$sequenceID[i]
  df <- read_tsv(paste0(RNA_data_path,"star/",runID,".SJ.out.tab"),
                 col_names = c("chromosome","start","end","strand","motif","annotated","unique_map","multi_map","junction_alignment"))
  df$width <- df$end - df$start
  df[[sampleID]] <- df$unique_map + df$multi_map
  df <- df %>% 
    filter(.[[sampleID]] >= 3, junction_alignment >= 10, nchar(chromosome) < 6) %>%
    dplyr::select(chromosome, start, end, width, strand, sampleID)
  if (nrow(junction_df) == 0) {
    junction_df <- df
  } else {
    junction_df <- full_join(junction_df, df, by = c("chromosome", "start", "end", "width", "strand"))
  }
}
junction_df[is.na(junction_df)] <- 0

junction_df <- ensembl_transcripts %>% 
  left_join(junction_df,by=c("region_start"="start","region_end"="end")) %>%
  filter(gene_biotype=="protein_coding",transcript_is_canonical==1,region_type=="intron") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

#write_tsv(as.data.frame(unique(junction_df$hgnc_symbol)),"~/deep_RNAseq/paper_analysis/MRSD/genelist")

calculate_MRSD <- function(median,targetcount){
  #res <- log10(targetcount/median)
  res <- targetcount/median/1000000
  return(res)
}

## change coverage
plot_list=list()
barplot_list=list()
#count_vec <- vector()
for (i in 1:length(unique(sample_list_represent$sample_type))){
  #geneID="ENSG00000187634"
  color=color_palette[i]
  sampletype=unique(sample_list_represent$sample_type)[i]
  reference_list <- sample_list_represent %>% filter(gradient==1,deduped=='yes',sample_type==sampletype)
  junction_matrix <- junction_df[,c("chr","region_start","region_end","width","ensembl_transcript_id","ensembl_gene_id","hgnc_symbol",reference_list$sequenceID)]
  junction_matrix[,reference_list$sequenceID] <- junction_matrix[,reference_list$sequenceID]/reference_list$total_reads
  junction_matrix$median <-  apply(junction_matrix[, reference_list$sequenceID], 1, median)
  combined_matrix <- tibble()
  for (count in c(10,20,50,100)){
    percentage=0.05 #loop manually
    matrix <- junction_matrix %>% group_by(ensembl_transcript_id,ensembl_gene_id,hgnc_symbol) %>%
      dplyr::summarise(base_count=quantile(median,percentage))
    matrix$percentage <- 1-percentage
    matrix$targetcount <- count
    matrix$sample_type <- sampletype
    matrix$mrsd <- sapply(matrix$base_count, calculate_MRSD, targetcount = count)
    #count_vec <- c(count_vec,length(matrix$mrsd[matrix$mrsd!="NaN" & matrix$mrsd<10000000000]))
    combined_matrix <- rbind(combined_matrix,matrix)
  }
  p1 <- ggplot(combined_matrix, aes(x=factor(targetcount), y=mrsd)) + # fill=name allow to automatically dedicate a color for each group
    geom_violin(fill=color)+
    scale_y_log10() +
    #ylim(0,100000000)+
    xlab(paste0("Junction coverage - ",sampletype))+
    ylab(NULL)+
    theme_classic()+
    theme(legend.position = "none")
  plot_list[[i]]=p1
  number_matrix <- combined_matrix %>% group_by(targetcount) %>%
    dplyr::summarise(feasible_gene=sum(ifelse(mrsd != "NaN" & mrsd!=Inf, 1, 0)))
  p2 <- ggplot(number_matrix, aes(x=factor(targetcount), y=feasible_gene)) +
    geom_bar(stat = "identity",fill=color)+
    ylim(0,15500)+
    xlab(paste0("Target junction coverage - ",sampletype))+
    ylab("No. of feasible genes")+
    theme_classic()+
    theme(legend.position = "none")
  barplot_list[[i]]=p2
  ##save tables as source data
  write_csv(combined_matrix,paste0("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/",sampletype,"-95.csv"))
}
p1=do.call(grid.arrange, c(plot_list, ncol=1))
p1
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/MRSD_junction.pdf",plot =p1,
       width = 150,height = 150,units ="mm")
p2=do.call(grid.arrange, c(barplot_list, ncol=2))
p2
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/MRSD_junction_bar.pdf",plot =p2,
       width = 150,height = 150,units ="mm")


###### change percentage
plot_list=list()
barplot_list=list()
#count_vec <- vector()
for (i in 1:length(unique(sample_list_represent$sample_type))){
  #geneID="ENSG00000187634"
  color=color_palette[i]
  sampletype=unique(sample_list_represent$sample_type)[i]
  reference_list <- sample_list_represent %>% filter(gradient==1,deduped=='yes',sample_type==sampletype)
  junction_matrix <- junction_df[,c("chr","region_start","region_end","width","ensembl_transcript_id","ensembl_gene_id","hgnc_symbol",reference_list$sequenceID)]
  junction_matrix[,reference_list$sequenceID] <- junction_matrix[,reference_list$sequenceID]/reference_list$total_reads
  junction_matrix$median <-  apply(junction_matrix[, reference_list$sequenceID], 1, median)
  combined_matrix <- tibble()
  for (percentage in c(0.5,0.25,0.15,0.05)){
    count=10 #loop manually
    matrix <- junction_matrix %>% group_by(ensembl_transcript_id,ensembl_gene_id) %>%
      dplyr::summarise(base_count=quantile(median,percentage))
    matrix$percentage <- 1-percentage
    matrix$targetcount <- count
    matrix$mrsd <- sapply(matrix$base_count, calculate_MRSD, targetcount = count)
    #count_vec <- c(count_vec,length(matrix$mrsd[matrix$mrsd!="NaN" & matrix$mrsd<10000000000]))
    combined_matrix <- rbind(combined_matrix,matrix)
  }
  p1 <- ggplot(combined_matrix, aes(x=factor(percentage), y=mrsd)) + # fill=name allow to automatically dedicate a color for each group
    geom_violin(fill=color)+
    #ylim(0,100000000)+
    xlab(sampletype)+
    ylab("logMRSD")+
    theme_classic()+
    theme(legend.position = "none")
  plot_list[[i]]=p1
  number_matrix <- combined_matrix %>% group_by(percentage) %>%
    dplyr::summarise(feasible_gene=sum(ifelse(mrsd != "NaN" & mrsd!=Inf, 1, 0)))
  p2 <- ggplot(number_matrix, aes(x=factor(percentage), y=feasible_gene)) +
    geom_bar(stat = "identity",fill=color)+
    ylim(0,16000)+
    xlab(paste0("Target percentage - ",sampletype))+
    ylab("No. of feasible genes")+
    theme_classic()+
    theme(legend.position = "none")
  barplot_list[[i]]=p2
  ##save tables as source data
}
p1=do.call(grid.arrange, c(plot_list, ncol=1))
p1
p2=do.call(grid.arrange, c(barplot_list, ncol=2))
p2
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/MRSD_junction_bar.pdf",plot =p2,
       width = 180,height = 150,units ="mm")


###correlation with MRSD
test_gene_list <- read_csv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/33c179d8-f8b6-459c-ab06-5c5e1587033c_medium_expressing.csv")
mrsd <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/msrd_results_75_10.tsv")
mrsd_deep <- combined_matrix %>% filter(percentage==0.25,targetcount==10) %>%
  left_join(mrsd,by=c("gene_version_id"="ensembl_id")) %>%
  filter(`MRSD(Fibroblasts)`!="-",transcript_type=="MANE") 
mrsd_deep$`MRSD(Fibroblasts)`=as.numeric(mrsd_deep$`MRSD(Fibroblasts)`)
mrsd_deep$mrsd_reg=mrsd_deep$`MRSD(Fibroblasts)`*1000000
mrsd_deep <- filter(mrsd_deep,mrsd < 100000000,mrsd_reg < 10000000)
cor_matrix <- cor(mrsd_deep[,c("mrsd","mrsd_reg")], method = "pearson")
p <- ggplot(mrsd_deep,aes(x=mrsd,y=mrsd_reg))+
  geom_point()
p
