makeLargeCloneDf <- function(inData,  
                          data_label_col = F, 
                          data_label_text = F) {
  if (!xor(data_label_col != F, data_label_text != F)) {
    stop("Specify one of data_label_col or data_label_text")
  }
    
  if(data_label_col != F) {
    outData <- inData %>%
      mutate(
        SpectrumGroup1=!!as.name(data_label_col)
      )
  }
  if(data_label_text != F) {
    outData <- inData %>%
      mutate(
        SpectrumGroup2=data_label_text
      )
  }
  
  group_totals <- outData %>%
    mutate(isLargeClone = case_when(
      t_alt_count > 1 ~ "LC",
      T ~ "Not LC"
    )) %>% 
    group_by(SpectrumGroup1) %>% 
    summarise(Total_n=n())
  outData <- outData %>%
    mutate(isLargeClone = case_when(
      t_alt_count > 1 ~ "LC",
      T ~ "Not LC"
    )) %>% 
    group_by(NULL) %>%
    group_by(SpectrumGroup1, isLargeClone) %>% 
    mutate(n=n()) %>% 
    select(SpectrumGroup1, isLargeClone, n) %>% 
    group_by(NULL) %>%
    complete(SpectrumGroup1, isLargeClone) %>% 
    left_join(group_totals) %>%
    replace_na(list("n"=0)) %>% 
    unique()
  return(outData)
}
