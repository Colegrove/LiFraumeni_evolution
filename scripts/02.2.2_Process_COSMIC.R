

COSMIC_table <- COSMIC_file  %>% 
  extract(MUTATION_AA, c("prot.ref","prot.pos","prot.alt"),
          "p.([A-Z])([0-9]+)([A-Z=*])", remove=FALSE, convert = TRUE) %>% 
  filter(!is.na(prot.pos)) %>% 
  group_by(GENE_NAME, prot.pos) %>% 
  summarize(Count=n()) %>% 
  group_by(GENE_NAME) %>%
  mutate(propMuts = Count/sum(Count))
COSMIC_transTable <- tibble(GENE_NAME=unique(COSMIC_table$GENE_NAME)) %>% 
  mutate(Gene = str_extract(GENE_NAME, "([A-Z0-9]*)"))
COSMIC_table <- COSMIC_table %>% 
  left_join(COSMIC_transTable, by=c("GENE_NAME")) %>% 
  mutate(Gene = factor(
    Gene, 
    levels=c("TP53")))

# Extract hotspot codons
OCHC<-COSMIC_table %>% 
  filter(propMuts >= .01 & Count >= 2) %>% 
  mutate(isHotspot = "Hotspot") %>% 
  write_delim(paste(out_prefix(),".OCHC.tsv", sep = ""), 
              delim = "\t")

COSMIC_table_forLolipop <- COSMIC_file  %>% 
  filter(!is.na(prot.pos)) %>% 
  filter(coding == "coding") %>% 
  filter(inMask == FALSE) %>% 
  filter(Variant_Type == "SNP") %>%
  group_by(GENE_NAME, prot.pos) %>% 
  summarize(Count=n()) %>% 
  group_by(GENE_NAME) %>%
  mutate(propMuts = Count/sum(Count))
COSMIC_transTable_forLolipop <- tibble(GENE_NAME=unique(COSMIC_table_forLolipop$GENE_NAME)) %>% 
  mutate(Gene = str_extract(GENE_NAME, "([A-Z0-9]*)"))
COSMIC_table_forLolipop <- COSMIC_table_forLolipop %>% 
  left_join(COSMIC_transTable, by=c("GENE_NAME")) %>% 
  mutate(Gene = factor(
    Gene, 
    levels=c("TP53")))
