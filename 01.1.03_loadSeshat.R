loadSeshat <- function(inSeshatFile, col_trans, path_trans) {
  outData <- 
    read_delim(inSeshatFile, delim="\t") %>% 
    # remove unnecessary columns
    select(cDNA_Variant, Comment_1_Frequency, Pathogenicity) %>% 
    # Pull out first non-space letter of Comment_1_Frequency
    mutate(C1F = str_sub(str_replace_all(Comment_1_Frequency, " ",""), 1,1)) %>%
    # Join on the classification by C1F
    left_join(col_trans, by="C1F") %>% 
    mutate(Pathogenicity = str_to_upper(Pathogenicity)) %>%
    # Join on the classification by Pathogenicity
    left_join(path_trans, by="Pathogenicity") %>% 
    mutate(C1F = factor(C1F,
                        levels = c("V","F","N","R","U","T")), 
           TP53_common1 = factor(TP53_common1), 
           TP53_common2 = factor(TP53_common2),
           TP53_path5 = factor(TP53_path5),
           TP53_path2 = factor(TP53_path2))
  return(outData)
}