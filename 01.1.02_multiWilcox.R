multi_wilcox <- function(data, groupVars, x_var, y_var) {
  outerGroups <- data %>% select(any_of(groupVars)) %>% unique() %>% mutate(groupNumber = row_number())
  print(outerGroups)
  data2 <- data %>% left_join(outerGroups)
  innerGroups <- unique(data %>% pull(x_var))
  outData = list()
  groups_compared = c()
  for (group1 in outerGroups$groupNumber) {
    print(group1)
    for (group2 in innerGroups) {
      for (group3 in innerGroups) {
        if (group2 > group3 & ! paste(group1, group2, group3) %in% groups_compared) {
          groups_compared = c(groups_compared, paste(group2, group3))
          comparison = paste(group1, group2, group3)
        }
        else if (group3 > group2 & ! paste(group1, group3, group2) %in% groups_compared) {
          groups_compared = c(groups_compared, paste(group3, group2))
          comparison = paste(group1, group3, group2)
        }
        if (group2 != group3) {
          test_data <- data2 %>% 
            filter(groupNumber == group1)
          for (yVar in y_var) {
            test_data_1 <- pull(test_data, yVar)[pull(test_data, x_var) == group2]
            test_data_2 <- pull(test_data, yVar)[pull(test_data, x_var) == group3]
            outData[[paste(comparison, yVar)]] = tryCatch(wilcox.exact(
              test_data_1, test_data_2, exact = T, conf.int = T) %>% 
                tidy() %>% 
                mutate(OuterGroup = paste(outerGroups %>% filter(groupNumber == group1) %>% select(-groupNumber) %>% toString()),
                       Group1 = group2,
                       Group2 = group3, 
                       yVar = yVar) %>% 
                select(OuterGroup, Group1, Group2, yVar, everything()),
              error = function(e) {
                message("An error occurred:\n", e)
                tibble(
                  "OuterGroup" = paste(outerGroups %>% 
                                         filter(groupNumber == group1) %>% 
                                         select(-groupNumber) %>% 
                                         toString()),
                  "Group1" = group2,
                  "Group2" = group3,
                  "yVar"=NA,
                  "estimate"=NA,
                  "statistic"=NA,
                  "p.value"=NA,
                  "conf.low"=NA,
                  "conf.high"=NA,
                  "method"=NA,
                  "alternative"=NA
                )
              },
              warning = function(w){
                message("A warning occured:\n", w)
              }
            )
            
          } 
        }
      }
    }
  }
  outData <- outData %>% 
    map_df(bind_rows)
  return(outData)
}