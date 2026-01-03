makeMutTypeDf <- function(inData, 
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
    group_by(SpectrumGroup1,Mutation_type) %>% 
    mutate(n=n()) %>% 
    select(SpectrumGroup1, Mutation_type,Total_n, n) %>% 
    unique()%>% 
    mutate(Mutation_type = case_when(
      Mutation_type == "Missense" ~ "Missense_Mutation", 
      Mutation_type == "Nonsense" ~ "Nonsense_Mutation",
      T ~ Mutation_type
    ))
  return(outData)
}
