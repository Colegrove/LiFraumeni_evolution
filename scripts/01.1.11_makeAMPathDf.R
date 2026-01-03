makeAMPathDf <- function(inData, 
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
        SpectrumGroup1=data_label_text
      )
  }
  outData <- outData %>%
    group_by(SpectrumGroup1) %>% 
    mutate(Total_n=n()) %>%
    group_by(isPathogenic_AM,SpectrumGroup1) %>% 
    mutate(n=n()) %>% 
    select(SpectrumGroup1, isPathogenic_AM,Total_n, n) %>% 
    unique() %>% 
    mutate(isPathogenic_AM=factor(isPathogenic_AM)) %>% 
    group_by(NULL)
}

makeAMPathDf2 <- function(inData, 
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
        SpectrumGroup1=data_label_text
      )
  }
  outData <- outData %>%
    group_by(SpectrumGroup1) %>% 
    mutate(Total_n=n()) %>%
    group_by(am_class,SpectrumGroup1) %>% 
    mutate(n=n()) %>% 
    select(SpectrumGroup1, am_class,Total_n, n) %>% 
    unique() %>% 
    mutate(am_class=factor(am_class)) %>% 
    group_by(NULL)
}
