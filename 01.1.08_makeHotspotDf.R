makeHotspotDf <- function(inData, inHotspotsData, 
                          data_label_col = F, 
                          data_label_text = F) {
  if (!xor(data_label_col != F, data_label_text != F)) {
    stop("Specify one of data_label_col or data_label_text")
  }
  outData <- inData %>%
    left_join(inHotspotsData %>%
                select(prot.pos,isHotspot)) %>%
    replace_na(list(isHotspot = "Not Hotspot"))
  if(data_label_col != F) {
    outData <- outData %>%
      mutate(
        SpectrumGroup2=!!as.name(data_label_col)
      )
  }
  if(data_label_text != F) {
    outData <- outData %>%
      mutate(
        SpectrumGroup2=data_label_text
      )
  }
  outData <- outData %>%
    mutate(isHotspot = case_when(
      prot.ref == "X" ~ "Not Hotspot",
      T ~ isHotspot
    )) %>% 
    group_by(NULL) %>%
    filter(!is.na(prot.pos)) %>%
    group_by(SpectrumGroup2, isHotspot) %>% 
    mutate(Total_n=n()) %>%
    mutate(n=n()) %>% 
    group_by(SpectrumGroup2) %>% 
    mutate(Total_n=n()) %>% 
    select(SpectrumGroup2, isHotspot, n, Total_n) %>% 
    group_by(NULL) %>%
    unique()
  return(outData)
}

