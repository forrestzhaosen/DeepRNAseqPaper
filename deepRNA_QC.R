library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(Hmisc) ##calculate weighted variance

## read in processed length distribution file
read_distribution <- read_tsv("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/52f64d6b-0c8d-4fe9-b0df-d4e278a2da40.raw.R1_fastqc.readlength.txt") %>%
  separate(col=`#Length`,into=c("start","end") , sep="-") %>%
  mutate(start=as.numeric(start),end=as.numeric(end)) %>%
  mutate(length=(start+end)/2-0.5)
sd <- sqrt(wtd.var(read_distribution$length,read_distribution$Count))
mean <- weighted.mean(read_distribution$length,read_distribution$Count)
read_tsv

##Star junction number count

column_name <- c("chromosome","start","end","strand","motif","annotated","unique_map","multi_map","junction_alignment")
star_junction_illumina <- fread("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-R1_singleton.d9296509-a8b6-4288-b26a-a9958a6e5899.qc_data/star/d9296509-a8b6-4288-b26a-a9958a6e5899.SJ.out.tab",col.names = column_name)
sum(star_junction_illumina$unique_map)
sum(star_junction_illumina$multi_map)

star_junction_ultima <- fread("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-RU1_singleton.06d1960a-f593-4557-9651-ae4dfebb2ffc.qc_data/star/06d1960a-f593-4557-9651-ae4dfebb2ffc.SJ.out.tab",col.names = column_name)
sum(star_junction_ultima$unique_map)
sum(star_junction_ultima$multi_map)
