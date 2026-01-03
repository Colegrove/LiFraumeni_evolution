
## table of mutations greater than VAF threshold
snp_table <- MAF_table_CHIP %>% 
  filter(NF <= inputs$maxNs & t_depth >= inputs$min_depth) %>%
  filter(VAF >= 0.3)

all_muts_somatic_germline <- MAF_table_CHIP %>%
  filter(NF <= inputs$maxNs & t_depth >= inputs$min_depth)
  

###
filt_maf_1_CHIP <- MAF_table_CHIP %>%
  # Apply filters
  ## Per sample depth filter
  # filter(! Tumor_Sample_Barcode %in% SampleDepthTest$Samp) %>%
  ## Per position depth filter
  filter(t_depth >= inputs$min_depth) %>% 
  
  ## Per position Ns filter
  filter(NF <= inputs$maxNs) %>%
  ## FILTERs from VCF
  filter(! grepl(filter_string, FILTER)) %>% 
  ## Masking filter
  filter(!inMask)


# count_SNPs_old_MAF_CHIP <- filt_maf_1_CHIP %>% 
#   group_by(mutPosition) %>% 
#   summarise(SNP_in = sum(VAF>= 0.9 | 
#                            (VAF >= 0.4 & VAF <= 0.6) 
#                          # | (coding == "non-coding" & VAF >= 0.2)
#                          ), 
#             mut_in = sum(VAF > 0), 
#             min_dp_pass = sum(t_depth >= inputs$min_depth)) %>% 
#   filter(mut_in > SNP_in) %>% 
#   filter(SNP_in > 0) 
# 
# count_SNPs_MAF_2_CHIP <- filt_maf_1_CHIP %>% 
#   #group_by(mutPosition, PatientID) %>% 
#   group_by(mutPosition) %>%
#   summarise(
#     high_VAF_in = sum(VAF>= 0.9 | (VAF >= 0.4 & VAF <= 0.6)),
#     min_dp_pass = sum(t_depth >= inputs$min_depth)) %>% 
#   filter(high_VAF_in > 0) %>%
#   filter(min_dp_pass == high_VAF_in) %>% 
#   #group_by(PatientID) %>% 
#   summarise(matching_high_VAF = sum(high_VAF_in == 2),
#             not_matching_high_VAF = sum(high_VAF_in != 2)) %>% 
#   filter(not_matching_high_VAF > 0)

filt_maf_2_CHIP <- filt_maf_1_CHIP %>% 
  # Apply filters
  ## Contamination Filter
  # filter(
  #   !mutPosition %in% count_SNPs_old_MAF$mutPosition
  # ) %>%
  # Filter to just variants in NEXUS
  # filter(in_NEXUS) %>% 
  ## Filter SNPs not matching in paired samples
  # filter(! Sample %in% count_SNPs_old_MAF_2$Sample) %>%
  ## Per-position VAF filter
  filter(!(VAF>= 0.3))

filt_maf_CHIP <- filt_maf_2_CHIP %>% 
  mutate(path_am_order = case_when(
    am_class == "likely_pathogenic" ~ 1,
    am_class == "ambiguous" ~ 2,
    am_class == "likely_benign" ~ 3, 
    T ~ NA
  )) %>% 
  ## need to add sample coding order in skyscraper plot script
  #group_by(Tumor_Sample_Barcode, coding) %>% 
  #arrange(desc(t_alt_count), path_am_order) %>% 
  #mutate(SampCodingOrder = rank(desc(t_alt_count), ties.method = "first")) %>% 
  #group_by(NULL) %>% 
  # left_join(CHIP_muts %>% mutate(isCHIP = "CHIP") %>% select(-c(mut_in))) %>% 
  write_delim(paste0(out_prefix(), ".mutations.filt.tsv"), delim = "\t", quote = "none")
