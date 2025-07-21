library(tidyverse)

nonstranded <- read_tsv("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/ultima/P-UDN027588_UDN027588-FBR1-RU1_singleton.69bf6358-1de4-4d82-836e-721febe21060.qc_data/rsem/69bf6358-1de4-4d82-836e-721febe21060.rsem.genes.results")
stranded <- read_tsv("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/ultima/P-UDN027588_UDN027588-FBR1-RU1_singleton.347bb671-3df1-4998-b880-1f2ba6878db4.qc_new/rsem/347bb671-3df1-4998-b880-1f2ba6878db4.rsem.genes.results")

compare=inner_join(nonstranded,stranded)

ggplot(compare,aes(x=`69bf6358-1de4-4d82-836e-721febe21060`,y=`347bb671-3df1-4998-b880-1f2ba6878db4`)) +
  geom_point()

sum(nonstranded$expected_count)
sum(stranded$expected_count)
