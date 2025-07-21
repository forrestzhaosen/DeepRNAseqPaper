source("~/deep_RNAseq/code/util.R")
RNA_data_path <- "/storage/liu/home/u247123/deep_RNAseq/freeze2023Oct/qc_data/"
sample_list <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/gradiant/Gradient_list.tsv")
sample_list_represent <- sample_list %>% filter(gradient==1,deduped=='yes',sample_type!='Blood')
sample_list_represent$total_reads <- sapply(sample_list_represent$runID,get_total_reads)
color_palette <- c("#999999","#E69F00","#56B4E9","#009E73",
                   "#F0E442","#0072B2","#D55E00","#CC79A7")
## read splicing matrix
junction_matrix_all <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/junction_matrix.tsv")
junction_matrix_all <- ensembl_transcripts %>% 
  mutate(gene_version_id = gsub("\\..*", "", gene_version_id)) %>%
  left_join(junction_matrix_all,by=c("region_start"="start","region_end"="end")) %>%
  left_join(hgnc_esemble,by= c("gene_version_id"="ensembl_gene_id")) %>% 
  filter(gene_biotype=="protein_coding",canonical==1,region_type=="intron") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

calculate_MRSD <- function(median,Nfrag,targetcount){
  res <- log(1-targetcount/(median*Nfrag),base=(Nfrag-1)/Nfrag)
  return(res)
}

plot_list=list()
barplot_list=list()
#count_vec <- vector()
for (i in 1:length(unique(sample_list_represent$sample_type))){
  i=3
  #geneID="ENSG00000187634"
  color=color_palette[i]
  sampletype=unique(sample_list_represent$sample_type)[i]
  reference_list <- sample_list_represent %>% filter(gradient==1,deduped=='yes',sample_type==sampletype)
  junction_matrix <- junction_matrix_all[,c("chrom","region_start","region_end","width","tx_id","gene_version_id","gene_name",reference_list$runID)]
  junction_matrix[,reference_list$runID] <- junction_matrix[,reference_list$runID]/reference_list$total_reads
  junction_matrix$median <-  apply(junction_matrix[, reference_list$runID], 1, median)
  combined_matrix <- tibble()
  for (count in c(5,10,15,20)){
    percentage=0.25 #loop manually
    matrix <- junction_matrix %>% group_by(tx_id,gene_version_id) %>%
      dplyr::summarise(base_count=quantile(median,percentage))
    matrix$percentage <- percentage
    matrix$targetcount <- count
    matrix$mrsd <- sapply(matrix$base_count, calculate_MRSD, Nfrag = 5000000000, targetcount = count)
    #count_vec <- c(count_vec,length(matrix$mrsd[matrix$mrsd!="NaN" & matrix$mrsd<10000000000]))
    combined_matrix <- rbind(combined_matrix,matrix)
  }
  p1 <- ggplot(combined_matrix, aes(x=factor(targetcount), y=mrsd)) + # fill=name allow to automatically dedicate a color for each group
    geom_violin(fill=color)+
    ylim(0,1000000000)+
    xlab(sampletype)+
    ylab("MRSD")+
    theme_classic()+
    theme(legend.position = "none")
  plot_list[[i]]=p1
  number_matrix <- combined_matrix %>% group_by(targetcount) %>%
    dplyr::summarise(feasible_gene=sum(ifelse(mrsd != "NaN" & mrsd < 10000000000, 1, 0)))
  p2 <- ggplot(number_matrix, aes(x=factor(targetcount), y=feasible_gene)) +
    geom_bar(stat = "identity",fill=color)+
    xlab(paste0("Desired junction coverage - ",sampletype))+
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
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/MRSD_junction.pdf",plot =p1,
       width = 200,height = 200,units ="mm")


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
