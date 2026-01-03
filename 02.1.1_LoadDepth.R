DP_table_full <- 
  read_delim(inputs$depth_file, delim="\t", col_types = cols(.default = "c" )) %>% 
  type_convert() %>% 
  ## twinstrand data Depth includes No_calls 22OCT24HLC
  mutate(NF = Ns/(DP))

if (! file.exists("inputs/bed_file_pan00367.csv") ) {
  testBed <- DP_table_full %>%
    select(Chr, Pos) %>%
    unique()
  ## Check for coding / non-coding position
  testBed$InBed <- check_bed_overlap(targ_bed6,
                                     testBed$Chr,
                                     testBed$Pos)
  ## Check for masking
  testBed$inMask <- check_bed_overlap(maskBed,
                                      testBed$Chr,
                                      testBed$Pos)
  
  ## Check for which gene
  testBed$Gene <- check_bed_overlap(targ_bed12,
                                    testBed$Chr,
                                    testBed$Pos,
                                    returnGenes = TRUE)
  
  # Save testBed so we don't need to rerun it
  testBed %>% write_csv("inputs/bed_file_pan00367.csv")
} else {
  testBed <- read_csv("inputs/bed_file_pan00367.csv")
}

# CHIP testBed
if (! file.exists("inputs/bed_file_CHIP.csv") ) {
  testBed_CHIP <- DP_table_full %>%
    dplyr::select(Chr, Pos) %>%
    unique()
  ## Check for coding / non-coding position
  testBed_CHIP$InBed <- check_bed_overlap(targ_bed6_CHIP,
                                     testBed_CHIP$Chr,
                                     testBed_CHIP$Pos)
  ## Check for masking
  testBed_CHIP$inMask <- check_bed_overlap(maskBed,
                                      testBed_CHIP$Chr,
                                      testBed_CHIP$Pos)
  
  ## Check for which gene
  testBed_CHIP$Gene <- check_bed_overlap(targ_bed12_CHIP,
                                    testBed_CHIP$Chr,
                                    testBed_CHIP$Pos,
                                    returnGenes = TRUE)
  
  # Save testBed so we don't need to rerun it
  testBed_CHIP %>% write_csv("inputs/bed_file_CHIP.csv")
} else {
  testBed_CHIP <- read_csv("inputs/bed_file_CHIP.csv")
}

# Mutagenesis testBed
if (! file.exists("inputs/bed_file_MUT.csv") ) {
  testBed_MUT <- DP_table_full %>%
    dplyr::select(Chr, Pos) %>%
    unique()
  ## Check for coding / non-coding position
  testBed_MUT$InBed <- check_bed_overlap(targ_bed6_MUT,
                                          testBed_MUT$Chr,
                                          testBed_MUT$Pos)
  ## Check for masking
  testBed_MUT$inMask <- check_bed_overlap(maskBed,
                                           testBed_MUT$Chr,
                                           testBed_MUT$Pos)
  
  ## Check for which gene
  testBed_MUT$Gene <- check_bed_overlap(targ_bed12_MUT,
                                         testBed_MUT$Chr,
                                         testBed_MUT$Pos,
                                         returnGenes = TRUE)
  
  # Save testBed so we don't need to rerun it
  testBed_MUT %>% write_csv("inputs/bed_file_MUT.csv")
} else {
  testBed_MUT <- read_csv("inputs/bed_file_MUT.csv")
}


# CHIP + Mutagenesis testBed
if (! file.exists("inputs/bed_file_CHIP_MUT.csv") ) {
  testBed_CHIP_MUT <- DP_table_full %>%
    dplyr::select(Chr, Pos) %>%
    unique()
  ## Check for coding / non-coding position
  testBed_CHIP_MUT$InBed <- check_bed_overlap(targ_bed6_CHIP_MUT,
                                         testBed_CHIP_MUT$Chr,
                                         testBed_CHIP_MUT$Pos)
  ## Check for masking
  testBed_CHIP_MUT$inMask <- check_bed_overlap(maskBed,
                                          testBed_CHIP_MUT$Chr,
                                          testBed_CHIP_MUT$Pos)
  
  ## Check for which gene
  testBed_CHIP_MUT$Gene <- check_bed_overlap(targ_bed12_CHIP_MUT,
                                        testBed_CHIP_MUT$Chr,
                                        testBed_CHIP_MUT$Pos,
                                        returnGenes = TRUE)
  
  # Save testBed so we don't need to rerun it
  testBed_CHIP_MUT %>% write_csv("inputs/bed_file_CHIP_MUT.csv")
} else {
  testBed_CHIP_MUT <- read_csv("inputs/bed_file_CHIP_MUT.csv")
}

DP_table <- DP_table_full %>%
  left_join(testBed) %>%
  mutate(coding = if_else(InBed, "coding", "non-coding"))
DP_table$NF


DP_summary <- DP_table %>%
  mutate(coding = if_else(InBed, "coding", "non-coding")) %>%
  filter(inMask == FALSE) %>%
  filter(NF <= inputs$maxNs) %>%
  group_by(Samp, coding) %>%
  summarise(meanDP = mean(DP),
            NTs = sum(DP)) %>% 
  group_by(coding) %>% 
  # extract(Tumor_Sample_Barcode, into = TSB_parts,
  #         regex = TSB_pattern,
  #         remove=F) %>% 
  
  mutate(averageDepth = mean(meanDP), stDevDepth = sd(meanDP),
         estNTs = meanDP * inputs$targetSize)%>% 
  mutate(highDepth = meanDP > averageDepth + 1.96 * stDevDepth,
         lowDepth = meanDP < averageDepth - 1.96 * stDevDepth) %>% 
  mutate(downsamplingRatio = if_else(highDepth, averageDepth / meanDP,if_else(lowDepth,0,1)),
         diffNTs = (NTs - estNTs)/estNTs)

show_excluded_depth_samples <- DP_summary %>%
  left_join(sample_data, by=c( "Samp" = "PatientID")) %>%
  filter(coding == 'coding') %>%
  filter(meanDP < 8000 | meanDP > 19000 | highDepth == TRUE | lowDepth == TRUE)
show_excluded_depth_samples

SampleDepthTest <- DP_summary %>%
  filter(coding == "coding") %>%
  # left_join(sampleDepthLimits) %>% 
  # mutate(DP_lim = case_when(
  #   is.na(DP_lim) == TRUE ~ (sampleDepthLimits %>% filter(Tissue == "TRUE"))$DP_lim,
  #   T ~ DP_lim
  # )) %>%
  # mutate(DP_lim = 1000) %>% 
  
  filter(meanDP < 1 | meanDP > 30000)
  #filter(meanDP < 8000 | meanDP > 19000)
SampleDepthTest

################## CHIP


DP_table_CHIP <- DP_table_full %>%
  left_join(testBed_CHIP) %>%
  mutate(coding = if_else(InBed, "coding", "non-coding"))


DP_summary_CHIP <- DP_table_CHIP %>%
  mutate(coding = if_else(InBed, "coding", "non-coding")) %>%
  filter(inMask == FALSE) %>%
  filter(NF <= inputs$maxNs) %>%
  group_by(Samp, coding) %>%
  summarise(meanDP = mean(DP),
            NTs = sum(DP)) %>% 
  group_by(coding) %>% 
  # extract(Tumor_Sample_Barcode, into = TSB_parts,
  #         regex = TSB_pattern,
  #         remove=F) %>% 
  
  mutate(averageDepth = mean(meanDP), stDevDepth = sd(meanDP),
         estNTs = meanDP * inputs$targetSize)%>% 
  mutate(highDepth = meanDP > averageDepth + 1.96 * stDevDepth,
         lowDepth = meanDP < averageDepth - 1.96 * stDevDepth) %>% 
  mutate(downsamplingRatio = if_else(highDepth, averageDepth / meanDP,if_else(lowDepth,0,1)),
         diffNTs = (NTs - estNTs)/estNTs)

DP_summary_CHIP


################## CHIP + MUT


DP_table_CHIP_MUT <- DP_table_full %>%
  left_join(testBed_CHIP_MUT) %>%
  mutate(coding = if_else(InBed, "coding", "non-coding"))


DP_summary_CHIP_MUT <- DP_table_CHIP_MUT %>%
  mutate(coding = if_else(InBed, "coding", "non-coding")) %>%
  filter(inMask == FALSE) %>%
  filter(NF <= inputs$maxNs) %>%
  group_by(Samp, coding) %>%
  summarise(meanDP = mean(DP),
            NTs = sum(DP)) %>% 
  group_by(coding) %>% 
  # extract(Tumor_Sample_Barcode, into = TSB_parts,
  #         regex = TSB_pattern,
  #         remove=F) %>% 
  
  mutate(averageDepth = mean(meanDP), stDevDepth = sd(meanDP),
         estNTs = meanDP * inputs$targetSize)%>% 
  mutate(highDepth = meanDP > averageDepth + 1.96 * stDevDepth,
         lowDepth = meanDP < averageDepth - 1.96 * stDevDepth) %>% 
  mutate(downsamplingRatio = if_else(highDepth, averageDepth / meanDP,if_else(lowDepth,0,1)),
         diffNTs = (NTs - estNTs)/estNTs)

DP_summary_CHIP_MUT