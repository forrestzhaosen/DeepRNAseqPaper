library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
library(readr)

RNA_data <- fread("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Phoenix/sample_run_data/sample_run_20220224.txt",na.strings="")
RNA_data_ultima <- filter(RNA_data,`Sequencing type`=="rnaseq_ultima")
RNA_data_illumina <- filter(RNA_data,`Sequencing type`=="rnaseq")

##some inspection
unique_counts_ultima <- sapply(RNA_data_ultima[,c("Family group","patientID","sample_ID",)], function(x) length(unique(x))) %>%
  as.data.frame()
setdiff(RNA_data_ultima$sample_ID,RNA_data_illumina$sample_ID)
View(table(RNA_data_illumina$sample_type))

##remove duplicate(RU1)
ultima_duplicate <- RNA_data_ultima$`External ID1`[grepl("RU2",RNA_data_ultima$`External ID1`)]
for (i in 1:length(ultima_duplicate)) {
  ultima_duplicate[i] <- paste0(substring(ultima_duplicate[i], 1, nchar(ultima_duplicate[i])-1), "1")
}
RNA_data_ultima <- RNA_data_ultima %>% 
  filter(!`External ID1` %in% ultima_duplicate)
RNA_data_illumina <- RNA_data_illumina %>%
  filter(`External ID1`!="UDN399289-NFR1-R1")
RNA_data_filtered <- RNA_data %>% 
  filter(!`External ID1` %in% ultima_duplicate) %>%
  filter(`External ID1`!="UDN399289-NFR1-R1")
  
##Filter samples sequenced by both techniques
sample_both <- intersect(RNA_data_ultima$sample_ID,RNA_data_illumina$sample_ID) %>% as.data.frame()
colnames(sample_both) <- c("sample_ID")

#sample_both <- left_join(sample_both,RNA_data_ultima[,c("sample_ID","External ID1")]) %>%
#  rename("ultima_id"="External ID1") %>%
#  left_join(RNA_data_illumina[,c("sample_ID","External ID1")]) %>%
#  rename("illumina_id"="External ID1")

sample_both <- left_join(sample_both,RNA_data_filtered) 

write_csv(RNA_data_filtered,"/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/freeze_2023feb/RNA_sample_list.csv")
write_csv(sample_both,"/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/freeze_2023feb/overlapping_sample_list.csv")
