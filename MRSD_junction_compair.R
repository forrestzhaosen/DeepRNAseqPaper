library(scales)

source('/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/code/RNA-seq/util.R')
RNA_data_path <- "/Volumes/zhaosen/freze2024Sep/qc_data/"

sample_list <-  read_tsv("/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/code/deep_paper/gradiant/all_list.tsv") %>%
  filter(samples_per_run==1)
sample_list_represent <- sample_list %>% filter(gradient==1,deduped=='yes',samples_per_run==1)
sample_list_represent$total_reads <- sapply(sample_list_represent$runID,get_total_reads)
color_palette <- c("#999999","#E69F00","#56B4E9","#009E73",
                   "#F0E442","#0072B2","#D55E00","#CC79A7")
## read splicing matrix
#junction_matrix_all <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/junction_matrix.tsv")
junction_df <- tibble()
for (i in 1:nrow(sample_list_represent)) {
  runID <- sample_list_represent$runID[i]
  sampleID <- sample_list_represent$sequenceID[i]
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


mrsd_75_10 <- read_tsv("/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/code/deep_paper/MRSD/official/msrd_results_75_10.tsv") %>%
  mutate(count=10,percent=0.25)
mrsd_75_20 <- read_tsv("/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/code/deep_paper/MRSD/official/msrd_results_75_20.tsv") %>%
  mutate(count=20,percent=0.25)
mrsd_95_10 <- read_tsv("/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/code/deep_paper/MRSD/official/msrd_results_95_10.tsv") %>%
  mutate(count=10,percent=0.05)
mrsd_95_20 <- read_tsv("/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/code/deep_paper/MRSD/official/msrd_results_95_20.tsv") %>%
  mutate(count=20,percent=0.05)
mrsd_orig <- rbind(mrsd_75_10,mrsd_75_20,mrsd_95_10,mrsd_95_20)

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
combined_matrix <- tibble()
correlation_matrix <- tibble()
for (i in 1:length(c("LCL","Fibroblast","Blood"))){
  #geneID="ENSG00000187634"
  #i=1
  color=color_palette[i]
  sampletype=c("LCL","Fibroblast","Blood")[i]
  reference_list <- sample_list_represent %>% filter(gradient==1,deduped=='yes',sample_type==sampletype)
  junction_matrix <- junction_df[,c("chr","region_start","region_end","width","ensembl_transcript_id","ensembl_gene_id","hgnc_symbol",reference_list$sequenceID)]
  junction_matrix[,reference_list$sequenceID] <- junction_matrix[,reference_list$sequenceID]/reference_list$total_reads
  junction_matrix$median <-  apply(junction_matrix[, reference_list$sequenceID], 1, median)
  #for (count in c(10,20)){}
  count=10
  percentage=0.05 #loop manually
  matrix <- junction_matrix %>% group_by(ensembl_transcript_id,ensembl_gene_id) %>%
    dplyr::summarise(base_count=quantile(median,percentage),.groups = "drop")
  matrix$percentage <- percentage
  matrix$targetcount <- count
  matrix$mrsd <- sapply(matrix$base_count, calculate_MRSD, targetcount = count)
  matrix <- filter(matrix,mrsd!="Inf")
  mrsd_matrix <- mrsd_orig[mrsd_orig[["count"]]==count & mrsd_orig[["percent"]]==percentage,]
  mrsd_matrix <- mrsd_matrix[,c("ensembl_id",sampletype)]
  names(mrsd_matrix) <- c("ensembl_id","mrsd_orig")
  matrix <- left_join(matrix,mrsd_matrix,by=c("ensembl_gene_id"="ensembl_id")) 
  matrix <- filter(matrix,mrsd_orig!="-")
  matrix$mrsd_orig <- as.numeric(matrix$mrsd_orig)
  #matrix$mrsd_orig <- log10(matrix$mrsd_orig*1000000)
  matrix$sample_type <- sampletype
  correlation_matrix <- rbind(correlation_matrix,matrix)
  matrix <- pivot_longer(matrix,cols = c(mrsd, mrsd_orig), names_to = "mrsd_type", values_to = "mrsd")
  #count_vec <- c(count_vec,length(matrix$mrsd[matrix$mrsd!="NaN" & matrix$mrsd<10000000000]))
  combined_matrix <- rbind(combined_matrix,matrix)
}
p1 <- ggplot(combined_matrix, aes(x=sample_type, y=mrsd,fill=mrsd_type)) + # fill=name allow to automatically dedicate a color for each group
    geom_violin(position="dodge")+
    scale_y_log10() +
    scale_fill_simpsons()+
    xlab(NULL)+
    ylab("MRSD (million reads)")+
    theme_classic()+
    theme(text = element_text(size = 20))
  
#p1=do.call(grid.arrange, c(plot_list, ncol=1))
p1
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/MRSD_compare.pdf",plot =p1,
       width = 200,height = 100,units ="mm")

plot_list <- list()
for (i in 1:length(c("LCL","Fibroblast","Blood"))){
  color=color_palette[i]
  sampletype=c("LCL","Fibroblast","Blood")[i]
  data <- correlation_matrix[correlation_matrix$sample_type==sampletype,]
  print(cor(data$mrsd,data$mrsd_orig,method='spearman'))
  p <- ggplot(data,aes(x=mrsd,y=mrsd_orig))+
    geom_point() +
    scale_x_log10(labels = label_number(accuracy = 1))+
    scale_y_log10(labels = label_number(accuracy = 1))+
    #xlim(0,5000000000)+
    xlab("MRSD deep")+
    ylab("MRSD original")+
    #geom_text(aes_(label=current_cor, x = 1000, y = 1000000))+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
    theme_classic()
  #p
  #geom_smooth(method=lm,formula=y ~ poly(x,2), color="red", fill="#69b3a2", se=TRUE)
  plot_list[[i]]=p
}
p1=do.call(grid.arrange, c(plot_list, ncol=3))
p1
ggsave("/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/code/deep_paper/figures/MRSD_correlation_R1.pdf",plot =p1,
       width = 200,height = 100,units ="mm")
###############################################################
###### change percentage, feasible genes bar plot###############################################################
###############################################################

barplot_list=list()
#count_vec <- vector()
for (i in 1:length(c("LCL","Fibroblast","Blood"))){
  #i=1
  combined_matrix <- tibble( sample_type = character(),     # integer type column
                             mrsd_type = character(), # character type column
                             feasible = numeric(),
                             percentage = numeric())
  #geneID="ENSG00000187634"
  color=color_palette[i]
  sampletype=c("LCL","Fibroblast","Blood")[i]
  reference_list <- sample_list_represent %>% filter(gradient==1,deduped=='yes',sample_type==sampletype)
  junction_matrix <- junction_df[,c("chr","region_start","region_end","width","ensembl_transcript_id","ensembl_gene_id","hgnc_symbol",reference_list$sequenceID)]
  junction_matrix[,reference_list$sequenceID] <- junction_matrix[,reference_list$sequenceID]/reference_list$total_reads
  junction_matrix$median <-  apply(junction_matrix[, reference_list$sequenceID], 1, median)

  for (percentage in c(0.25,0.05)){
    #percentage=0.25
    count=10
    matrix <- junction_matrix %>% group_by(ensembl_transcript_id,ensembl_gene_id) %>%
      dplyr::summarise(base_count=quantile(median,percentage),.groups = "drop")
    matrix$percentage <- percentage
    matrix$targetcount <- count
    matrix$mrsd <- sapply(matrix$base_count, calculate_MRSD, targetcount = count)
    #matrix <- filter(matrix,mrsd!="Inf")
    mrsd_matrix <- mrsd_orig[mrsd_orig[["count"]]==count & mrsd_orig[["percent"]]==percentage,]
    mrsd_matrix <- mrsd_matrix[,c("ensembl_id",sampletype)]
    names(mrsd_matrix) <- c("ensembl_id","mrsd_orig")
    matrix <- left_join(matrix,mrsd_matrix,by=c("ensembl_gene_id"="ensembl_id")) %>% filter(!is.na(mrsd_orig))
    #matrix <- filter(matrix,mrsd_orig!="-")
    combined_matrix <- combined_matrix %>% add_row(sample_type = sampletype, mrsd_type = "mrsd", feasible = sum(matrix$mrsd!=Inf),percentage=percentage) #/nrow(matrix)
    combined_matrix <- combined_matrix %>% add_row(sample_type = sampletype, mrsd_type = "mrsd_orig", feasible = sum(matrix$mrsd_orig!="-"),percentage=percentage) #/nrow(matrix)
  }
  combined_matrix$mrsd_type <- factor(combined_matrix$mrsd_type)
  combined_matrix$mrsd_type <- factor(combined_matrix$mrsd_type, levels = rev(levels(combined_matrix$mrsd_type)))
  combined_matrix$coverage <- 1-combined_matrix$percentage
  p2 <- ggplot(combined_matrix, aes(x=factor(coverage), y=feasible,fill=mrsd_type)) +
    geom_bar(position="dodge", stat="identity")+
    #ylim(0,16000)+
    scale_fill_simpsons()+
    xlab(paste0("Percentage covered - ",sampletype))+
    ylab("Percent of feasible genes")+
    theme_classic()+
    theme(legend.position = "none")
  barplot_list[[i]]=p2
  ##save tables as source data
}
p2=do.call(grid.arrange, c(barplot_list, ncol=3))
p2
ggsave("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/figures/MRSD_compare_bar.pdf",plot =p2,
       width = 200,height = 100,units ="mm")
write_csv(combined_matrix,"/storage/liu/home/u247123/deep_RNAseq/paper_analysis/MRSD/percentage_compare.csv")
