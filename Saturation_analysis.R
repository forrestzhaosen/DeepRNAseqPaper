library(minpack.lm)

source("~/deep_RNAseq/code/util.R")

RNA_data_path <- "/storage/liu/home/u247123/deep_RNAseq/freeze2023Oct/qc_data/"
sample_list <- read_tsv("/storage/liu/home/u247123/deep_RNAseq/paper_analysis/gradiant/all_list.tsv")
sample_list_represent <- sample_list %>% filter(gradient==1,deduped=='yes')

### lapply and unlist to calculate raw count and deduped count for each sample
sample_list_represent$count_raw <- lapply(
  1:nrow(sample_list_represent), 
  function(i) {
    get_normal_gene_counts(sample_list_represent$pre_dedup_runID[i], sample_list_represent$sample_type[i])
  }
)
sample_list_represent$count_raw <- unlist(sample_list_represent$count_raw)

sample_list_represent$count_dedup <- lapply(
  1:nrow(sample_list_represent), 
  function(i) {
    get_normal_gene_counts(sample_list_represent$runID[i], sample_list_represent$sample_type[i])
  }
)
sample_list_represent$count_dedup <- unlist(sample_list_represent$count_dedup)

### calculate duplication rate
sample_list_represent$duplication_rate <- (sample_list_represent$count_raw - sample_list_represent$count_dedup)/sample_list_represent$count_raw
sample_list_represent$total_reads <- sapply(sample_list_represent$pre_dedup_runID,get_total_reads)
sample_list_represent$total_reads_dedup <- sapply(sample_list_represent$runID,get_total_reads)
sample_list_represent$total_reads_ratio <-sample_list_represent$total_reads_dedup/sample_list_represent$total_reads
write_csv(sample_list_represent,"~/deep_RNAseq/paper_analysis/gradiant/duplication_rate_new.csv")

p <- ggplot(sample_list_represent,aes(x=total_reads,y=duplication_rate,group=paper_ID, color=sample_type))+
  geom_point() +
  geom_line() +
  scale_color_jama(name = NULL)+
  #ylim(0,ylim[i]) + geom_hline(yintercept = yintercept[i],linetype = "dashed", color = "black", linewidth = 0.6)+
  #xlim(0,5000000000)+
  scale_x_log10(breaks=c(50000000,100000000,500000000,1000000000,3000000000)) +
  xlab("Sequencing depth")+
  ylab("Duplication rate")+
  #theme(text = element_text(size = 23),
  #      axis.text.x = element_text(angle = 45, hjust = 1))+
  #theme(plot.margin = unit(c(2, 2, 2,2), "lines")) +
  theme_classic()
p

### different amount of RNA in diffrent smaple/run
model_func <- function(Y,X,N) {
  B <- X/N
  expected_untouched <- B * ((B - 1) / B)^Y
  return(1 - (B - expected_untouched) / Y)
}

sample_list_represent$samples_per_run <- as.numeric(sample_list_represent$samples_per_run)
X_list=list()
for (sample in unique(sample_list_represent$paper_ID)) {
 depth_vector <- sample_list_represent$total_reads[sample_list_represent$paper_ID==sample]
 dup_vector <- sample_list_represent$duplication_rate[sample_list_represent$paper_ID==sample]
 N <- sample_list_represent$samples_per_run[sample_list_represent$paper_ID==sample]
 fit <- nlsLM(dup_vector ~ model_func(depth_vector, X,N), start=list(X=10000000000))
 params <- coef(fit)
 X_list[[sample]]=params
}
X_list

# Create a data frame with X values and a separate column for each 'a'
X_value <- seq(0, 20000000000, by = 200000)  # Adjust the range and step as needed
df <- data.frame(X = X_value)
#df <- data.frame(X = rep(X_value, each = length(X_list)))
for (i in 1:length(X_list)) {
  # Access the name of the first element
  col_name <- names(X_list)[i]
  # Access the value of the first element
  a <- X_list[[i]]
  df[, col_name] <- 1- ((a - 1) / a)^X_value
}
# Reshape the data frame from wide to long format
df_long <- tidyr::pivot_longer(df, cols = 2:10, names_to = "sampleID", values_to = "Y")
df_long$sample_type <- substr(df_long$sampleID, 1, nchar(df_long$sampleID) - 1)
# Create the plot
p <- ggplot(df_long, aes(x = X, y = Y, group=sampleID, color=sample_type)) +
  geom_line() +
  labs(x = "Sequencing depth", y = "Saturation rate", color = "Sample type") +
  theme_classic()+
  geom_hline(yintercept = 0.9,linetype = "dashed", color = "black", linewidth = 0.6)+
  scale_color_jama()
  #scale_color_discrete(name = "Sample type")

p  
  
