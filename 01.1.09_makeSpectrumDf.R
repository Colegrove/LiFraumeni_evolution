subs_simplification <- tibble(
  BaseVal = c("C>A","C>T","C>G","G>A","G>T","G>C","T>A","T>C","T>G","A>T","A>C","A>G"),
  ModVal  = c("C>A","C>T","C>G","C>T","C>A","C>G","T>A","T>C","T>G","T>A","T>G","T>C"))

makeSpectrumDf <- function(inData, 
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
    # filter for substitutions only
    filter(Variant_Type == "SNP") %>%
    mutate(Mutational_Event = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) %>% 
    left_join(subs_simplification, by=c("Mutational_Event"="BaseVal")) %>% 
    rename("Change" = "ModVal") %>%
    group_by(SpectrumGroup1) %>% 
    mutate(Total_n=n()) %>%
    group_by(SpectrumGroup1,Change,Total_n) %>% 
    count() %>% 
    select(SpectrumGroup1, Change, Total_n,n) %>% 
    group_by(NULL) %>% 
    unique()
  return(outData)
}
