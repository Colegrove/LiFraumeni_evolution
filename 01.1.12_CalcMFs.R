calculateMFs <- function(inMaf, inDepths) {
  MF_nums <- inMaf %>% 
    group_by(SampID, coding) %>%
    summarise(num_muts = n(), 
              am_total = sum(am_pathogenicity),
              am_burden = sum(am_pathogenicity * t_alt_count),
              mb = sum(t_alt_count),
              #num_path = sum(TP53_path2 == 1),
              num_lc = sum(t_alt_count > 1),
              num_hotspot = sum(isHotspot == "Hotspot"),
              mb_hotspot = sum(t_alt_count[isHotspot == "Hotspot"]),
              #num_common = sum(TP53_common2 == 1),
              num_exon58 = sum(
                Exon_Number %in% c("5/11",
                                   "6/11",
                                   "7/11",
                                   "8/11")),
              num_driver = sum(isDriver == "Driver"),
              num_boost_driver = sum(boostDM_class == T),
              mb_driver = sum(t_alt_count[isDriver == "Driver"]),
              mb_boost_driver = sum(t_alt_count[boostDM_class == T]),
              num_pathogenic_am =  sum(isPathogenic_AM == "Pathogenic"),
              num_lc_pathogenic_am = sum(isPathogenic_AM == "Pathogenic" & t_alt_count > 1),
              mb_pathogenic_am = sum(t_alt_count[isPathogenic_AM == "Pathogenic"]),
              #num_pathogenic_seshat = sum(isPathogenic_seshat == "Pathogenic"),
              #num_lc_pathogenic_seshat = sum(isPathogenic_seshat == "Pathogenic" & t_alt_count > 1),
              #mb_pathogenic_seshat = sum(t_alt_count[isPathogenic_seshat == "Pathogenic"]),
              #num_common_seshat = sum(isCommon_seshat == "Common"),
              #mb_common_seshat = sum(t_alt_count[isCommon_seshat == "Common"])
    )
  
  
  MF_nums <- 
    MF_nums %>% 
    bind_rows(
      MF_nums %>% 
        group_by(SampID) %>%
        summarise(num_muts = sum(num_muts), 
                  am_total = sum(am_total),
                  am_burden = sum(am_burden),
                  mb = sum(mb),
                  #num_path = sum(num_path),
                  num_lc = sum(num_lc),
                  num_hotspot = sum(num_hotspot),
                  mb_hotspot = sum(mb_hotspot),
                  #num_common = sum(num_common),
                  num_exon58 = sum(num_exon58),
                  num_driver = sum(num_driver),
                  num_boost_driver = sum(num_boost_driver),
                  mb_driver = sum(mb_driver),
                  mb_boost_driver = sum(mb_boost_driver),
                  num_pathogenic_am =  sum(num_pathogenic_am),
                  num_lc_pathogenic_am = sum(num_lc_pathogenic_am),
                  mb_pathogenic_am = sum(mb_pathogenic_am),
                  #num_pathogenic_seshat = sum(num_pathogenic_seshat),
                  #num_lc_pathogenic_seshat = sum(num_lc_pathogenic_seshat),
                  #mb_pathogenic_seshat = sum(mb_pathogenic_seshat),
                  #num_common_seshat = sum(num_common_seshat),
                  #mb_common_seshat = sum(mb_common_seshat)
        ) %>% 
        mutate(coding = "total")
    )
  
  MF_denoms <- inDepths  %>% 
    
    mutate(SampID1 = Samp,
           SampID = Samp,    
           PatientID = Samp  
    ) %>%
    
    # separate_wider_delim(Samp,delim = ".",names = c("SampID1","LibID","SeqDate","ProcDate","clipping"),cols_remove = F) %>% 
    # extract(SampID1, into = c("SampID","FragmentationMethod"),regex = "^(.*)_([EZFSON]{3})$",remove = F) %>%
    # extract(SampID, into = c("PatientID"), regex = "^(.*)_[A-Z]{3}_[0-9]{2}$", remove = F) %>% 
    
    ## no sample_data? 24OCT24HLC
    #left_join(sample_data, by=c("PatientID")) %>% 
    
    group_by(SampID, coding) %>%
    summarise(denom = sum(DP),
              meanDP = mean(DP),
              numPositions = n())
  MF_denoms <- 
    MF_denoms %>% 
    bind_rows(
      MF_denoms %>% 
        group_by(SampID) %>% 
        summarise(denom = sum(denom),
                  numPositions = sum(numPositions)) %>% 
        mutate(
          coding = "total",
          meanDP = denom / numPositions
        )
    )
  outData <- MF_nums %>%
    full_join(MF_denoms, 
              by = c("SampID", "coding")) %>% 
    replace_na(list("num_muts" = 0, 
                    "am_total" = 0,
                    "am_total" = 0,
                    "am_burden" = 0,
                    "mb" = 0,
                    "num_path" = 0,
                    "num_lc" = 0,
                    "num_hotspot" = 0,
                    "mb_hotspot" = 0,
                    "num_common" = 0,
                    "num_exon58" = 0,
                    "num_driver" = 0,
                    "num_boost_driver" = 0,
                    "mb_driver" = 0,
                    "mb_boost_driver" = 0,
                    "num_pathogenic_am" = 0,
                    "num_lc_pathogenic_am" = 0,
                    "mb_pathogenic_am" = 0,
                    "num_pathogenic_seshat" = 0,
                    "num_lc_pathogenic_seshat" = 0,
                    "mb_pathogenic_seshat" = 0,
                    "num_common_seshat" = 0,
                    "mb_common_seshat" = 0
    )) %>% 
    mutate(MF = num_muts / denom, 
           MF_AM = am_total / denom,
           MBF_AM = am_burden / denom,
           MBF = mb / denom, 
           #MF_path2 = num_path / denom,
           MF_LC = num_lc / denom,
           MF_hotspot = num_hotspot / denom,
           MBF_hotspot = mb_hotspot / denom,
           #MF_common = num_common / denom,
           MF_exon58 = num_exon58 / denom,
           MF_driver = num_driver / denom,
           MBF_driver = mb_driver / denom,
           MF_boost_driver = num_boost_driver / denom,
           MBF_boost_driver = mb_boost_driver / denom,
           MF_pathogenic_am = num_pathogenic_am / denom,
           MBF_pathogenic_am = mb_pathogenic_am / denom,
           #MF_pathogenic_seshat = num_pathogenic_seshat / denom,
           #MBF_pathogenic_seshat= mb_pathogenic_seshat / denom,
           #MF_common_seshat = num_common_seshat / denom,
           #MBF_common_seshat = mb_common_seshat / denom
    )
  return(outData)
}

