main_MFs <- 
  calculateMFs(
    filt_maf, 
    DP_filt
  ) %>% 
  write_delim(paste0(out_prefix(), ".MF_table_main.tsv"), delim = "\t", quote = "none") %>% 
  extract(SampID, into = c("PatientID"), regex = "^(.*)_[A-Z]{3}_[0-9]{2}$", remove = F) %>% 
  
  ## no sample_data? 18OCT24HLC
  #left_join(sample_data, by=c("PatientID")) %>% 
  
  write_delim(paste0(out_prefix(), ".MF_table_main_with_patient_data.tsv"), delim = "\t", quote = "none")
  

