
MAF_table_CHIP <-
  loadMaf(inputs$concat_maf_file,
          inputs$seshat_long_file,
          testBed_CHIP,
          seshat_col_trans,
          seshat_path_trans,
          variant_clasification_table,
          alphamissense,
          intOGen
  )  %>%
  # extract(Tumor_Sample_Barcode, into = TSB_parts,
  #         regex = TSB_pattern,
  #         remove=F) %>%
  
  ## t_NC depth will not be present in maf file from twinstrand vcf conversion
  ## t_depth is total depth, t_ref_count is reference reads
  ## t_alt_count is mutant reads
  ## can compute t_NC by subtracting reference/mutant reads from total depth
  ## Will need to compute no call frequency from genome mut 28OCT24HLC
  mutate(#VAF=t_alt_count/t_depth) %>%
         #NF = t_NC / (t_depth + t_NC)) 
          VAF = t_alt_count / (t_ref_count + t_alt_count),
          t_NC = t_depth - t_ref_count - t_alt_count,
          NF = t_NC / (t_depth)) %>%
  
  
  left_join(OCHC, by=c("Hugo_Symbol"="Gene", "prot.pos")) %>%
  replace_na(list(isHotspot = "Not Hotspot")) %>%
  mutate(
    isDriver = if_else((
      ( Mutation_type == "Missense_Mutation" &
          isHotspot == "Hotspot") |
        (Mutation_type %in% c("Indel","Splice","Nonsense_Mutation") )
    ),"Driver","Non-driver","Non-driver")
  ) %>%
  ## keep columns commented out below, but populate sampleID 18OCT24HLC
  mutate(SampID1 = Tumor_Sample_Barcode,
         SampID = Tumor_Sample_Barcode,    
         PatientID = Tumor_Sample_Barcode  
  ) %>%
  
  # separate_wider_delim(Tumor_Sample_Barcode,delim = ".",names = c("SampID1","LibID","SeqDate","ProcDate","clipping"),cols_remove = F) %>%
  # extract(SampID1, into = c("SampID","FragmentationMethod"),regex = "^(.*)_([EZFSON]{3})$",remove = F) %>%
  # extract(SampID, into = c("PatientID"), regex = "^(.*)_[A-Z]{3}_[0-9]{2}$", remove = F) %>%
  left_join(sample_data, by=c("PatientID")) %>%
  write_delim(paste0(out_prefix(), ".MafWithActivity.tsv"), delim = "\t")

