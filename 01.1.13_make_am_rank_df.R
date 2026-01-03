make_am_rank_df <- function(inData, 
                            data_label_col = F, 
                            data_label_text = F) {
  if (!xor(data_label_col != F, data_label_text != F)) {
    stop("Specify one of data_label_col or data_label_text")
  }
  if(data_label_col != F) {
    outData <- inData %>%
      mutate(
        group_name=!!as.name(data_label_col)
      )
  }
  if(data_label_text != F) {
    outData <- inData %>%
      mutate(
        group_name=data_label_text
      )
  }
  outData <- outData %>% 
    filter(!is.na(am_pathogenicity)) %>% 
    group_by(group_name) %>% 
    mutate(order = rank(am_pathogenicity,ties.method = "first")) %>% 
    mutate(maxorder = max(order)) %>% 
    mutate(order2 = order / maxorder) %>% 
    select(group_name, am_pathogenicity, order, maxorder, order2)
  return(outData)
}
