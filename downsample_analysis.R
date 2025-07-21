library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
library(ggsci)
library(RColorBrewer) # pallete
library(gridExtra)
source("~/deep_RNAseq/code/util.R")

OMIM_list <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/disease_gene_list/genemap2_10032023.txt",skip = 3)  %>% filter(!is.na(Phenotypes),!is.na(`Approved Gene Symbol`))
OMIM_gene_list <- OMIM_list$`Ensembl Gene ID` %>% unique() %>% unlist()

RNA_data_path <- "/storage/liu/RNAseq/freze2024Sep/qc_data/"
sample_list <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/gradiant/all_list.tsv")
sample_list$total_reads <- sapply(sample_list$runID,get_total_reads)

#metadata <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/ultima_illumina_overlapping_fastq/formal_comparison/metadata.txt") ### two version of metadata
#metadata_FBR <- metadata[metadata$sample_type=="FBR",]
#metadata_NFR <- metadata[metadata$sample_type=="NFR",]


##########################################
#### correlation analysis ####
##########################################
sample_pair <- sample_list[!is.na(sample_list$downsample_runID),]
cor <- character()
for (i in 1:nrow(sample_pair)){
  runID_ultima <- sample_pair$runID[i]
  runID_illumina <- sample_pair$downsample_runID[i]
  expression_ultima <- read_tsv(paste0(RNA_data_path,"rnaseqc/",runID_ultima,".gene_reads.gct"),skip=2) %>%
    filter(!str_detect(Name, "PAR")) %>%
    mutate(Name = gsub("\\..*", "", Name)) %>%
    left_join(hgnc_esemble,by= c("Name"="ensembl_gene_id")) %>%
    filter(gene_biotype=="protein_coding",hgnc_id!="HGNC:30554")
  names(expression_ultima)[3] <- "ultima_count"
  expression_illumina <- read_tsv(paste0(RNA_data_path,"rnaseqc/",runID_illumina,".gene_reads.gct"),skip=2) %>%
    filter(!str_detect(Name, "PAR")) %>%
    mutate(Name = gsub("\\..*", "", Name)) %>%
    left_join(hgnc_esemble,by= c("Name"="ensembl_gene_id")) %>%
    filter(gene_biotype=="protein_coding",hgnc_id!="HGNC:30554")
  names(expression_illumina)[3] <- "illumina_count"
  expression_comb <- cbind(expression_ultima,expression_illumina[,"illumina_count"])
  expression_comb <- filter(expression_comb,illumina_count>50 & ultima_count>50)
  text <- round(cor(expression_comb$ultima_count,expression_comb$illumina_count, method = "spearman"),3)
  cor[[i]]=paste0(text)
}
sample_list$cor <- cor
Nexpressed <- t(sapply(sample_list$runID,get_gene_reads))
sample_list$Nexpressed1 <- Nexpressed[,1]

Nexpressed <- t(sapply(sample_list$downsample_runID,get_gene_reads))
sample_list$Nexpressed2 <- Nexpressed[,1]
write_tsv(sample_list,"/storage/liu/home/u247123/deep_RNAseq/paper_analysis/gradiant/3sample_downsample.tsv")


## do the correlation plot
plot_list <- list()
for (i in 1:nrow(sample_pair)){
  runID_ultima <- sample_pair$runID[i]
  runID_illumina <- sample_pair$downsample_runID[i]
  expression_ultima <- read_tsv(paste0(RNA_data_path,"rnaseqc/",runID_ultima,".gene_reads.gct"),skip=2) %>%
    filter(!str_detect(Name, "PAR")) %>%
    mutate(Name = gsub("\\..*", "", Name)) %>%
    left_join(hgnc_esemble,by= c("Name"="ensembl_gene_id")) %>%
    filter(gene_biotype=="protein_coding",hgnc_id!="HGNC:30554")
  names(expression_ultima)[3] <- "ultima_count"
  expression_illumina <- read_tsv(paste0(RNA_data_path,"rnaseqc/",runID_illumina,".gene_reads.gct"),skip=2) %>%
    filter(!str_detect(Name, "PAR")) %>%
    mutate(Name = gsub("\\..*", "", Name)) %>%
    left_join(hgnc_esemble,by= c("Name"="ensembl_gene_id")) %>%
    filter(gene_biotype=="protein_coding",hgnc_id!="HGNC:30554")
  names(expression_illumina)[3] <- "illumina_count"
  expression_comb <- cbind(expression_ultima,expression_illumina[,"illumina_count"])
  expression_comb <- filter(expression_comb,illumina_count>50 & ultima_count>50)
  #current_cor <- as.character(cor[i])
  p <- ggplot(expression_comb,aes(x=ultima_count,y=illumina_count))+
    geom_point(alpha = 0.5) +
    scale_x_log10()+
    scale_y_log10()+
    #xlim(0,5000000000)+
    xlab("Gene read count - true")+
    ylab("Gene read count - downsample")+
    #geom_text(aes_(label=current_cor, x = 1000, y = 1000000))+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
    theme_classic()
  #p
  #geom_smooth(method=lm,formula=y ~ poly(x,2), color="red", fill="#69b3a2", se=TRUE)
  plot_list[[i]]=p
}
p1=do.call(grid.arrange, c(plot_list, ncol=3))
p1
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/bg1205_downsample.pdf",plot =p1,
       width = 200,height = 150,units ="mm")


###### qc  anlaysis
get_exonic_rate <- function(runID){
  if (runID!="" & !is.na(runID)){
    rnaseqc <- read_tsv(paste0(RNA_data_path,runID,".multiqc_data/","multiqc_rna_seqc.txt"))
    return(rnaseqc$`Exonic Rate`)
  } else{
    return("")
  }
}
sample_list$exonic_rate <- sapply(sample_list$runID,get_exonic_rate)

plot_list <- list()

for (i in 1:3){
  sample <- c("BG1205-LBR1","BG1476-ICR1","UDN761469-FBR2")[i]
  data <- sample_list[sample_list$deduped=="yes"&sample_list$sample_ID==sample,]
  p <- ggplot(data,aes(x=total_reads,y=exonic_rate)) +
    geom_point() +
    geom_line() +
    ylim(0,1)+
    ylab("Exonic rate")+
    xlab("Sequencing depth")+
    scale_x_log10() +
    #geom_line() +
    guides(color = "none") +
    theme_classic()+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 25, hjust = 1))
  plot_list[[i]]=p
}
#plot_list[[1]]
p1=do.call(grid.arrange, c(plot_list, ncol=3))
p1
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/exonic_rate.pdf",plot =p1,
       width = 400,height = 130,units ="mm")
