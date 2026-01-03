
###
filt_maf_1 <- MAF_table %>%
  # Apply filters
  filter(t_depth >= inputs$min_depth) %>%
  ## Per position Ns filter
  filter(NF <= inputs$maxNs) %>%
  ## FILTERs from VCF
  filter(! grepl(filter_string, FILTER)) %>% 
  ## Masking filter
  filter(!inMask)

# count_SNPs_old_MAF <- filt_maf_1 %>%
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
# count_SNPs_old_MAF %>% print(n=Inf)
# 
# count_SNPs_MAF_2 <- filt_maf_1 %>%
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

filt_maf_2 <- filt_maf_1 %>% 
  # Apply filters
  ## Contamination Filter
  # filter(
  #   !mutPosition %in% count_SNPs_old_MAF$mutPosition
  # ) %>%
  ## Filter SNPs not matching in paired samples
  # filter(! Sample %in% count_SNPs_old_MAF_2$Sample) %>%
  
  ## Per-position VAF filter
  filter(!(VAF>= 0.3))

filt_maf <- filt_maf_2 %>% 
  mutate(path_am_order = case_when(
    am_class == "likely_pathogenic" ~ 1,
    am_class == "ambiguous" ~ 2,
    am_class == "likely_benign" ~ 3, 
    T ~ NA
  )) %>% 
  group_by(Tumor_Sample_Barcode, coding, Hugo_Symbol) %>% 
  arrange(desc(t_alt_count), path_am_order) %>% 
  mutate(SampCodingOrder = rank(dplyr::desc(t_alt_count), ties.method = "first")) %>% 
  group_by(NULL) %>% 
  write_delim(paste0(out_prefix(), ".mutations.filt.tsv"), delim = "\t", quote = "none")

