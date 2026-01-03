DP_filt <- 
  DP_table %>% 
  filter(! Samp %in% SampleDepthTest$Samp) %>%
  filter(inMask == FALSE) %>%
  filter(NF <= inputs$maxNs) %>% 
  filter(DP >= inputs$min_depth)

## CHIP
DP_filt_CHIP <- 
  DP_table_CHIP %>% 
  filter(! Samp %in% SampleDepthTest$Samp) %>%
  filter(inMask == FALSE) %>%
  filter(NF <= inputs$maxNs) %>% 
  filter(DP >= inputs$min_depth)

## CHIP + MUT
DP_filt_CHIP_MUT <- 
  DP_table_CHIP_MUT %>% 
  filter(! Samp %in% SampleDepthTest$Samp) %>%
  filter(inMask == FALSE) %>%
  filter(NF <= inputs$maxNs) %>% 
  filter(DP >= inputs$min_depth)

