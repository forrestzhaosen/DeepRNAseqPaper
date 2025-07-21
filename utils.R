get_missplicing_table_donors <- function(exons_df, sj_df, exons_source) {
  # step 1. get normal splicing
  NS <- sj_df[exons_df, 
              on = .(chrom, donor_pos = exon_donor_pos, acceptor_pos = nextExon_acceptorPos, strand)]
  NS[is.na(splicing_event_class), splicing_event_class := 'normal splicing']
  NS[,exon_donor_pos := donor_pos ]
  # step 2. get exon skipping events
  ## anchored at previous donor 
  ES_d <- sj_df[annotated_acceptor == 1][exons_df,allow.cartesian=TRUE,on = .(chrom, strand, 
                                                                              donor_pos = prevExon_donorPos)
  ][data.table::between(acceptor_pos,ES_d_lb,ES_d_ub,incbounds = F)]
  
  ES_d <- exons_source[, .(tx_id, donor_pos, donor_exon_no = region_no)][ES_d, on = .(tx_id, donor_pos)]
  ES_d <- exons_source[, .(tx_id, acceptor_pos, acceptor_exon_no = region_no)][ES_d, on = .(tx_id, acceptor_pos)]
  
  ## anchored at next acceptor
  ES_a <- sj_df[annotated_donor == 1][exons_df,allow.cartesian=TRUE,on = .(chrom, strand, 
                                                                           acceptor_pos = nextExon_acceptorPos)
  ][data.table::between(donor_pos,ES_a_lb,
                        ES_a_ub,incbounds = F)]
  
  ES_a <- exons_source[, .(tx_id, donor_pos, donor_exon_no = region_no)][ES_a, on = .(tx_id, donor_pos)]
  ES_a <- exons_source[, .(tx_id, acceptor_pos, acceptor_exon_no = region_no)][ES_a, on = .(tx_id, acceptor_pos)]
  # get exon skipping info - exons skipped, #nts omitted from transcript ( to calculate frame)
  ES <- rbind(ES_a, ES_d, fill = T)
  ES <- ES[!is.na(acceptor_exon_no) & !is.na(donor_exon_no)]
  ES[, skipped_exons_count := acceptor_exon_no - donor_exon_no - 1]
  
  ES[, paste_first := donor_exon_no + 1]
  ES[, paste_second := acceptor_exon_no - 1]
  
  nts_skipped_join <- exons_source[, .(tx_id, region_no, region_width)][ES[, .(tx_id, exon_donor_pos,exon_acceptor_pos, 
                                                                               paste_first, paste_second)], 
                                                                        on = .(tx_id, region_no >= paste_first, 
                                                                               region_no <= paste_second)]
  
  nts_skipped_join <- nts_skipped_join[, .(skipped_exons_nt = sum(region_width)), 
                                       by = .(tx_id, exon_donor_pos, exon_acceptor_pos, 
                                              paste_first = region_no, paste_second = region_no.1)]
  
  ES <- nts_skipped_join[ES, on = .(tx_id, exon_donor_pos,exon_acceptor_pos, 
                                    paste_first, paste_second)]
  
  ES[paste_second == paste_first, paste_second := NA]
  ES[, skipped_exons_id := do.call(paste, c(.SD, sep="-")), .SDcols= paste_first:paste_second]
  ES[, skipped_exons_id := gsub('-NA', '',skipped_exons_id )]
  
  ES[splicing_event_class == 'normal splicing', splicing_event_class := 'exon skipping (annotated)']
  ES[, missplicing_inframe := ifelse(skipped_exons_nt %% 3 == 0, 1, 0)]
  # step 3. get cryptic splicing events
  ## cryptic donor
  CSS <-sj_df[exons_df,allow.cartesian = TRUE, on = .(chrom, strand, 
                                                      acceptor_pos = nextExon_acceptorPos)
  ][donor_pos != exon_donor_pos & 
      data.table::between(donor_pos, css_lb_d, css_ub_d)]
  
  # calculate distance of cryptic to authentic splice site
  CSS[, cryptic_pos := donor_pos]
  CSS[, cryptic_distance := strand_mult * (donor_pos - exon_donor_pos) ]
  CSS[, missplicing_inframe := ifelse(cryptic_distance %% 3 == 0, 1, 0)]
  CSS[, cryptic_distance := ifelse(cryptic_distance > 0, cryptic_distance + 1, cryptic_distance)]
  CSS[splicing_event_class == 'normal splicing' | annotated_donor == 1, splicing_event_class := 'alternative donor (annotated)']
  
  # step 4. combine into one mis-splicing event table
  ms_table <- rbind(NS, ES, CSS, fill = T)
  setDT(ms_table)
  ms_table[, splice_site_pos := exon_donor_pos]
  cols <- c('gene_name', 'tx_id', 'canonical', 'chrom', 'splice_site_pos',
            'exon_donor_pos', 'exon_acceptor_pos', 
            'region_width', 'region_no', 'strand',
            'donor_pos', 'acceptor_pos', 
            'total_reads', 'splicing_event_class',
            'skipped_exons_count', 'skipped_exons_id', 'skipped_exons_nt',
            'cryptic_pos', 'cryptic_distance', 'missplicing_inframe')
  ms_table <- distinct(ms_table[, ..cols])
  
  setorder(ms_table, tx_id, region_no, -total_reads)
  
  
  return(ms_table)
}



get_missplicing_table_acceptors <- function(exons_df, sj_df, exons_source) {
  # step 1. get normal splicing
  NS <- sj_df[exons_df, 
              on = .(chrom, donor_pos = prevExon_donorPos, acceptor_pos = exon_acceptor_pos, strand)]
  NS[is.na(splicing_event_class), splicing_event_class := 'normal splicing']
  NS[, exon_acceptor_pos := acceptor_pos]
  
  # step 2. get exon skipping events
  ## anchored at previous donor 
  ES_d <- sj_df[annotated_acceptor == 1][exons_df,allow.cartesian=TRUE,on = .(chrom, strand, 
                                                                              donor_pos = prevExon_donorPos)
  ][data.table::between(acceptor_pos,ES_d_lb,ES_d_ub,incbounds = F)]
  
  ES_d <- exons_source[, .(tx_id, donor_pos, donor_exon_no = region_no)][ES_d, on = .(tx_id, donor_pos)]
  ES_d <- exons_source[, .(tx_id, acceptor_pos, acceptor_exon_no = region_no)][ES_d, on = .(tx_id, acceptor_pos)]
  
  ## anchored at next acceptor
  ES_a <- sj_df[annotated_donor == 1][exons_df,allow.cartesian=TRUE,on = .(chrom, strand, 
                                                                           acceptor_pos = nextExon_acceptorPos)
  ][data.table::between(donor_pos,ES_a_lb,
                        ES_a_ub,incbounds = F)]
  
  ES_a <- exons_source[, .(tx_id, donor_pos, donor_exon_no = region_no)][ES_a, on = .(tx_id, donor_pos)]
  ES_a <- exons_source[, .(tx_id, acceptor_pos, acceptor_exon_no = region_no)][ES_a, on = .(tx_id, acceptor_pos)]
  
  # get exon skipping info - exons skipped, #nts omitted from transcript ( to calculate frame)
  ES <- rbind(ES_a, ES_d, fill = T)
  ES <- ES[!is.na(acceptor_exon_no) & !is.na(donor_exon_no)]
  ES[, skipped_exons_count := acceptor_exon_no - donor_exon_no - 1]
  
  ES[, paste_first := donor_exon_no + 1]
  ES[, paste_second := acceptor_exon_no - 1]
  
  nts_skipped_join <- exons_source[, .(tx_id, region_no, region_width)][ES[, .(tx_id, exon_donor_pos,exon_acceptor_pos, 
                                                                               paste_first, paste_second)], 
                                                                        on = .(tx_id, region_no >= paste_first, region_no <= paste_second)]
  
  nts_skipped_join <- nts_skipped_join[, .(skipped_exons_nt = sum(region_width)), 
                                       by = .(tx_id, exon_donor_pos, exon_acceptor_pos, paste_first = region_no, paste_second = region_no.1)]
  
  ES <- nts_skipped_join[ES, on = .(tx_id, exon_donor_pos,exon_acceptor_pos, 
                                    paste_first, paste_second)]
  
  ES[paste_second == paste_first, paste_second := NA]
  ES[, skipped_exons_id := do.call(paste, c(.SD, sep="-")), .SDcols= paste_first:paste_second]
  ES[, skipped_exons_id := gsub('-NA', '',skipped_exons_id )]
  
  ES[splicing_event_class == 'normal splicing', splicing_event_class := 'exon skipping (annotated)']
  ES[, missplicing_inframe := ifelse(skipped_exons_nt %% 3 == 0, 1, 0)]
  
  # step 3. get cryptic splicing events
  ## cryptic acceptor
  CSS <-sj_df[exons_df,allow.cartesian=TRUE, on = .(chrom, strand, 
                                                    donor_pos = prevExon_donorPos)
  ][acceptor_pos != exon_acceptor_pos & 
      data.table::between(acceptor_pos, css_lb_a, css_ub_a)]
  
  # calculate distance of cryptic to authentic splice site
  CSS[, cryptic_pos := acceptor_pos]
  CSS[, cryptic_distance := strand_mult * (acceptor_pos - exon_acceptor_pos) ]
  CSS[, missplicing_inframe := ifelse(cryptic_distance %% 3 == 0, 1, 0)]
  CSS[, cryptic_distance := ifelse(cryptic_distance < 0, cryptic_distance - 1, cryptic_distance)]
  CSS[splicing_event_class == 'normal splicing'| annotated_acceptor == 1, splicing_event_class := 'alternative acceptor (annotated)']
  
  # step 4. combine into one mis-splicing event table
  ms_table <- rbind(NS, ES, CSS, fill = T)
  setDT(ms_table)
  ms_table[, splice_site_pos := exon_acceptor_pos]
  cols <- c('gene_name', 'tx_id','canonical', 'chrom', 'splice_site_pos',
            'exon_donor_pos', 'exon_acceptor_pos', 
            'region_width', 'region_no', 'strand',
            'donor_pos', 'acceptor_pos', 
            'total_reads', 'splicing_event_class',
            'skipped_exons_count', 'skipped_exons_id', 'skipped_exons_nt',
            'cryptic_pos', 'cryptic_distance', 'missplicing_inframe')
  ms_table <- distinct(ms_table[, ..cols])
  
  
  
  setorder(ms_table, tx_id, region_no, -total_reads)
  
  return(ms_table)
}
