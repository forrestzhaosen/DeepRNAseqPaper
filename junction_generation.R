library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(readxl)
#library(biomaRt) 


###read star junction files in to a jointed table
#sample_list <- read_excel("/storage/liu/home/u247123/deep_RNAseq/freeze2023feb/RNA_sample_list_with_path.xlsx")
column_name <- c("chromosome","start","end","strand","motif","annotated","unique_map","multi_map","junction_alignment")
star_junction <- read_tsv("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/ultima/Pai, Abhav_UDN761469-FBR3-RU1_singleton.b556d2d3-b375-49f0-ad57-03ff7e27a020.qc_data/star/b556d2d3-b375-49f0-ad57-03ff7e27a020.SJ.out.tab",col_names = column_name) %>%
  filter(strand != 0) %>%
  mutate(strand=gsub(1,"+",strand)) %>%
  mutate(strand=gsub(2,"-",strand))
write_tsv(star_junction,"/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/ultima/Pai, Abhav_UDN761469-FBR3-RU1_singleton.b556d2d3-b375-49f0-ad57-03ff7e27a020.qc_data/star/junction.tsv")
