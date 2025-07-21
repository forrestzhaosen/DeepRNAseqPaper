library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(biomaRt)
library(readxl)
library(ggplot2)

ensembl <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                      dataset="hsapiens_gene_ensembl")
hgnc_esemble <- getBM(attributes=c("ensembl_gene_id", "hgnc_id", "hgnc_symbol","band"), values="*", mart=ensembl) %>% filter(hgnc_symbol!="") 
hgnc_enst <- getBM(attributes=c("hgnc_symbol", "hgnc_id", "ensembl_transcript_id","band"), values="*", mart=ensembl)  %>% filter(hgnc_id!="") 
clingen <- read_excel("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/gene_list/Clingen-Gene-Disease-Summary-2023-02-17.xlsx") %>%
  unlist()
ID_genelist <- read_excel("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/gene_list/Intellectual disability.xlsx",sheet = "unique_genes")  %>%
  unlist()

ultima_data <- read_tsv("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-RU1_singleton.06d1960a-f593-4557-9651-ae4dfebb2ffc.qc_data/rsem/06d1960a-f593-4557-9651-ae4dfebb2ffc.rsem.genes.results") %>%
  mutate(gene_id = gsub("\\..*", "", gene_id)) %>%
  left_join(hgnc_esemble,by= c("gene_id"="ensembl_gene_id")) %>%
  mutate(sequence_type = "deep") %>%
  mutate(clingen = hgnc_symbol %in% clingen,ID = hgnc_symbol %in% ID_genelist)

illumina_data <- read_tsv("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-R1_singleton.d9296509-a8b6-4288-b26a-a9958a6e5899.qc_data/rsem/d9296509-a8b6-4288-b26a-a9958a6e5899.rsem.genes.results")  %>%
  mutate(gene_id = gsub("\\..*", "", gene_id)) %>%
  left_join(hgnc_esemble,by= c("gene_id"="ensembl_gene_id")) %>%
  mutate(sequence_type = "normal") %>%
  mutate(clingen = hgnc_symbol %in% clingen,ID = hgnc_symbol %in% ID_genelist)

ultima_data_filter <- ultima_data %>% filter(expected_count>100)
illumina_data_filter <- illumina_data %>% filter(expected_count>100)
count(filter(ultima_data,clingen == TRUE,expected_count>100))
intersect(filter(ultima_data,ID == TRUE,expected_count>100)$gene_id,filter(illumina_data,ID == TRUE,expected_count>100)$gene_id) %>% length

plot_data <- bind_rows(ultima_data,illumina_data) %>%
  filter(ID == TRUE)

ggplot(plot_data,aes(x=sequence_type, y=expected_count,fill=sequence_type)) + 
  geom_violin(position="dodge", alpha=0.5) +
  ylim(0,10000)+
  geom_hline(yintercept = 200)+
  theme_classic() + theme(legend.position="none")


