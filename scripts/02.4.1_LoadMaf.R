
MAF_table <-
  loadMaf(inputs$concat_maf_file,
          alphamissense,
          intOGen
  )  %>% 
  mutate(VAF = t_alt_count / (t_ref_count + t_alt_count),
         t_NC = t_depth - t_ref_count - t_alt_count,
         NF = t_NC / (t_depth) ) %>%

  #left_join(OCHC, by=c("Hugo_Symbol"="Gene", "prot.pos")) %>% print()
  #replace_na(list(isHotspot = "Not Hotspot")) %>%
  # mutate(
  #   isDriver = if_else((
  #     ( Mutation_type == "Missense_Mutation" &
  #         isHotspot == "Hotspot") |
  #       (Mutation_type %in% c("Indel","Splice","Nonsense_Mutation") )
  #   ),"Driver","Non-driver","Non-driver")
  # ) %>%
  
  ## keep columns commented out below, but populate sampleID 18OCT24HLC
  mutate(SampID1 = Tumor_Sample_Barcode,
         SampID = Tumor_Sample_Barcode,    
         PatientID = Tumor_Sample_Barcode  
  ) %>%
  left_join(sample_data, by=c("PatientID")) 
#%>%
  #write_delim(paste0(out_prefix(), ".MafWithActivity.tsv"), delim = "\t")


MAF_table %>% print(width = Inf)
