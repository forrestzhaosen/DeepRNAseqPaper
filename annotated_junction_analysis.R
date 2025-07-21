library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(biomaRt)
library(readxl)
library(ggplot2)

OMIM_list <- read_tsv("/Users/zhaosen/PUMCH Dropbox/PUMCH团队文件夹/Public-Databases/OMIM/2022-02-24/genemap2.txt",skip = 3)  %>%
  filter(!is.na(Phenotypes),!is.na(`Approved Gene Symbol`))
OMIM_gene_list <- OMIM_list$`Approved Gene Symbol` %>% unique() %>% unlist()
clingen <- read_excel("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/gene_list/Clingen-Gene-Disease-Summary-2023-02-17.xlsx") %>%
  unlist()
ID_genelist <- read_excel("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/gene_list/Intellectual disability.xlsx",sheet = "unique_genes")  %>%
  unlist()


column_name <- c("chromosome","start","end","strand","motif","annotated","unique_map","multi_map","junction_alignment")
ultima_data <- fread("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-RU1_singleton.06d1960a-f593-4557-9651-ae4dfebb2ffc.qc_data/star/06d1960a-f593-4557-9651-ae4dfebb2ffc.SJ.out.tab",col.names = column_name) %>%
  mutate(total_reads = unique_map+multi_map,chrom = gsub("chr", "", chromosome), strand=gsub("1","+",gsub("2","-",strand))) %>%
  filter(junction_alignment > 2, unique_map != 0)
illumina_data <- fread("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-R1_singleton.d9296509-a8b6-4288-b26a-a9958a6e5899.qc_data/star/d9296509-a8b6-4288-b26a-a9958a6e5899.SJ.out.tab",col.names = column_name) %>%
  mutate(total_reads = unique_map+multi_map,chrom = gsub("chr", "", chromosome), strand=gsub("1","+",gsub("2","-",strand))) %>%
  filter(junction_alignment > 2, unique_map != 0)

ensembl_transcripts <- fread("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/annotation_file/ensembl_introns_exons_v104.tsv")
ensembl_transcripts_canonical_junction <- ensembl_transcripts[region_type == 'intron' & canonical==1 & nexon != 1] ##nexon != 1 not necessary because I filtered for intron
omim_junction <- ensembl_transcripts_canonical_junction[gene_name %in% OMIM_gene_list]
clingen_junction <- ensembl_transcripts_canonical_junction[gene_name %in% clingen]

ultima_omim <- omim_junction %>% left_join(ultima_data,by=c("region_start"="start","region_end"="end")) %>%
  mutate(total_reads = unique_map+multi_map) %>%
  mutate(unique_map=replace_na(unique_map,0),multi_map=replace_na(multi_map,0),total_reads=replace_na(total_reads,0)) %>%
  mutate(unique5 = unique_map > 10,total5 = total_reads > 10)

ultima_omim_collapse <- ultima_omim %>% group_by(tx_id,gene_name) %>%
  summarise(proportion_unique=sum(unique5)/n(),proportion_total=sum(total5)/n()) %>%
  mutate(type="deep")

illumina_omim <- omim_junction %>% left_join(illumina_data,by=c("region_start"="start","region_end"="end")) %>%
  mutate(total_reads = unique_map+multi_map) %>%
  mutate(unique_map=replace_na(unique_map,0),multi_map=replace_na(multi_map,0),total_reads=replace_na(total_reads,0)) %>%
  mutate(unique5 = unique_map > 10,total5 = total_reads > 10)

illumina_omim_collapse <- illumina_omim %>% group_by(tx_id,gene_name) %>%
  summarise(proportion_unique=sum(unique5)/n(),proportion_total=sum(total5)/n()) %>%
  mutate(type="normal")

omim_comparison <-  bind_rows(ultima_omim_collapse,illumina_omim_collapse)

ggplot(omim_comparison,aes(x=type, y=proportion_unique,fill=type)) + 
  geom_violin(position="dodge", alpha=0.5) +
  stat_summary(fun = "mean",
               geom = "point",size=5,
               color = "green") +
  theme(text = element_text(family = "Arial",size = 20,color="#000000"),axis.text = element_text(family = "Arial",size = 20,color="#000000"))+
  #geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  #ylim(0,10000)+
  theme(legend.position="none")
##See how many genes have more than two canonical transcripts
junction_reads_dup <- junction_reads[, .(tx_id, gene_name)]
junction_reads_unique <- junction_reads_dup[!duplicated(junction_reads_dup)]
gene_number <- groupby(junction_reads_unique,gene_name) %>% summarise(count=)
##ensembl starts and ends

ensembl_transcripts <- ensembl_transcripts[region_type %in% c('exon', 'intron')]
ensembl_ann_starts <- unique(ensembl_transcripts[, .(chrom, start = region_start)])
ensembl_ann_starts[, annotated_start := 1]
ensembl_ann_ends <- unique(ensembl_transcripts[, .(chrom, end = region_end)])
ensembl_ann_ends[, annotated_end := 1]
ensembl_ann_sjs <- unique(ensembl_transcripts[, .(chrom, start = region_start, end = region_end)])
ensembl_ann_sjs[, annotated_sj := 1]

##merge our junction with ensembl
junction_data <- ensembl_ann_starts[junction_data,on = .(chrom, start)]
junction_data <- ensembl_ann_ends[junction_data,on = .(chrom, end)]
junction_data <- ensembl_ann_sjs[junction_data,on = .(chrom, start,end)]
junction_data[is.na(junction_data)] <- 0
junction_data <- junction_data[annotated_start == 1 | annotated_end == 1]

#annotate with basic mis-splicing event categories
matrix <- data.table(strand = c(), annotated_intron = c(), annotated_start = c(), annotated_end = c(), splicing_event_class = c())
matrix[, strand := rep(c('+', '-'), 4)]
matrix[, annotated_sj := c(1,1,rep(0, 6))]
matrix[, annotated_start := c(1,1,1,1,1,0,0,1)]
matrix[, annotated_end := c(1,1,1,1,0,1,1,0)]
matrix[, splicing_event_class := c(rep('normal splicing',2), rep('exon skipping', 2), rep('cryptic acceptor', 2), rep('cryptic donor', 2))]

junction_data <- matrix[junction_data, on = .(strand, annotated_sj, annotated_start, annotated_end)]

junction_data[strand == '+' & annotated_start == 1, annotated_donor := 1]
junction_data[strand == '-' & annotated_end == 1, annotated_donor := 1]
junction_data[strand == '+' & annotated_end == 1, annotated_acceptor := 1]
junction_data[strand == '-' & annotated_start == 1, annotated_acceptor := 1]
junction_data[is.na(annotated_donor), annotated_donor := 0]
junction_data[is.na(annotated_acceptor), annotated_acceptor := 0]

fwrite(junction_data,"/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-RU1_singleton.06d1960a-f593-4557-9651-ae4dfebb2ffc.qc_data/ultima_junction_processed.tsv")

# Check for duplicates
if (nrow(junction_data[,.(N = .N), by = .(chrom, start, end)][N > 1]) != 0) { stop("Duplicate rows in merged data") }

#rename start/end to donor/acceptor
junction_data[strand == '+', `:=` (donor_pos = start, acceptor_pos = end)]
junction_data[strand == '-', `:=` (donor_pos = end, acceptor_pos = start)]

# Load Exons reference - engineering columns to help with variant processing
exons <- ensembl_transcripts[region_type == 'exon']
exons[, donor_pos := ifelse(strand == '+', region_end + 1, region_start - 1)]
exons[, acceptor_pos := ifelse(strand == '+', region_start - 1, region_end + 1)]
exons[, nextExon_acceptorPos :=  lead(acceptor_pos), by = .(tx_id)]
exons[, prevExon_donorPos :=  lag(donor_pos), by = .(tx_id)]
exons[, strand_mult := ifelse(strand == '+', 1, -1)]

# remove single exon transcripts
exons <- exons[!is.na(nextExon_acceptorPos) | !is.na(prevExon_donorPos)]

# for first and last exons annotate neighbouring splice sites slighlty differently
exons[is.na(nextExon_acceptorPos), nextExon_acceptorPos := donor_pos]
exons[is.na(prevExon_donorPos), prevExon_donorPos := acceptor_pos]

# get bounds of transcript
exons[, tx_start := acceptor_pos[region_no == 1], by = tx_id]
exons[, tx_end := donor_pos[region_no == nexon], by = tx_id]

# make helper columns
exons[, ES_d_lb := min(donor_pos,tx_end), by = 1:nrow(exons)]
exons[, ES_d_ub := max(donor_pos,tx_end), by = 1:nrow(exons)]
exons[, ES_a_lb := min(acceptor_pos,tx_start), by = 1:nrow(exons)]
exons[, ES_a_ub := max(acceptor_pos,tx_start), by = 1:nrow(exons)]
exons[, css_lb_d := min(acceptor_pos, nextExon_acceptorPos), by = 1:nrow(exons)]
exons[, css_ub_d := max(acceptor_pos, nextExon_acceptorPos), by = 1:nrow(exons)]
exons[, css_lb_a := min(donor_pos, prevExon_donorPos), by = 1:nrow(exons)]
exons[, css_ub_a := max(donor_pos, prevExon_donorPos), by = 1:nrow(exons)]

# filter to unique intron-exon pairs
exons_dup <- exons[, .(chrom, region_start, region_end, strand, nextExon_acceptorPos, prevExon_donorPos)]
exons_unique <- exons[!duplicated(exons_dup),]

exons_df_donor <- exons[last_region== 0, 
                        .(gene_name, tx_id,canonical,
                          chrom, region_start, region_end, region_width, region_no, 
                          nexon, strand, last_region,
                          exon_donor_pos = donor_pos, 
                          exon_acceptor_pos = acceptor_pos, 
                          nextExon_acceptorPos, 
                          prevExon_donorPos, 
                          strand_mult, tx_start, tx_end, 
                          ES_d_lb, ES_d_ub, ES_a_lb, ES_a_ub, 
                          css_lb_d, css_ub_d, css_lb_a, css_ub_a)]

ms_donors <- get_missplicing_table_donors(exons_df_donor,junction_data,exons)

exons_df_acceptor <- exons[region_no != 1, 
                           .(gene_name, tx_id,canonical,
                             chrom, region_start, region_end, region_width, region_no, 
                             nexon, strand, last_region,
                             exon_donor_pos = donor_pos, 
                             exon_acceptor_pos = acceptor_pos, 
                             nextExon_acceptorPos, 
                             prevExon_donorPos, 
                             strand_mult, tx_start, tx_end, 
                             ES_d_lb, ES_d_ub, ES_a_lb, ES_a_ub, 
                             css_lb_d, css_ub_d, css_lb_a, css_ub_a)]
ms_acceptors <- get_missplicing_table_acceptors(exons_df_acceptor,junction_data,exons)

ms_donors[, ss_type := 'donor']
ms_acceptors[, ss_type := 'acceptor']
ms_all <- distinct(rbind(ms_donors, ms_acceptors))

ms_all[, event_count := .N, by = list(splice_site_pos, tx_id)]
ms_all[, event_rank := rowid(splice_site_pos), by = list(tx_id)]
ms_all[splicing_event_class != 'normal splicing', missplicing_event_rank := rowid(splice_site_pos), by = list(tx_id)]
ms_all[missplicing_inframe == 1, missplicing_inframe := TRUE]
ms_all[missplicing_inframe == 0, missplicing_inframe := FALSE]
ms_all[splicing_event_class == 'normal splicing', missplicing_inframe := TRUE]
ms_all[, missplicing_inframe := as.logical(missplicing_inframe)]

ms_all2 <- ms_all[,.(splice_site_pos, 
                     gene_name,
                     tx_id, 
                     canonical,
                     ss_type,
                     exon_no = region_no,
                     strand, 
                     splicing_event_class, 
                     event_rank, 
                     missplicing_inframe,
                     total_reads,
                     skipped_exons_count,
                     skipped_exons_id,
                     cryptic_distance, 
                     chr = chrom, 
                     donor_pos, 
                     acceptor_pos, 
                     assembly = 'hg38',
                     transcript_type = 'ensembl')]

fwrite(ms_all2, '/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/example_data/P-UDN523337_UDN523337-FBR1-RU1_singleton.06d1960a-f593-4557-9651-ae4dfebb2ffc.qc_data/mis_splice_ensembl.tsv.gz', sep = '\t', nThread = 4, compress = 'gzip')
