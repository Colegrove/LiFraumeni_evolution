tissue_filtering <- filt_maf %>% 
  filter(!is.na(prot.pos)) %>% 
  filter(Subject == "Patient") %>%
  filter(Tissue %in% all_samples) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(coding_from_maf == "coding") %>%
  filter(Variant_Classification != "Splice_Site") %>%
  filter(inMask == FALSE) %>% 
  filter(Variant_Type=="SNP")



## proportion of mutations in the DBD compared to non-DBD
tissue_filtering %>% 
  mutate(DBD = if_else(prot.pos >= 109 & prot.pos <= 288, "DBD", "non-DBD")) %>% 
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  group_by(DBD, LFS) %>%
  summarise(count = n()) %>% print(width = Inf)

exon2 = abs(7676521-	7676594)
exon3 = abs(7676382-	7676403)
exon4 = abs(7675994-	7676272)
exon5 = abs(7675053-7675236) 
exon6 = abs(7674859	- 7674971)
exon7 = abs(7674181- 7674290)
exon8 = abs(7673701 - 	7673837)
exon9 = abs(7673535	-7673608)
exon10 = abs(7670609-7670715)
exon11 = abs(7669609 -	7669690)

exon5 + exon6 + exon7 + exon8
exon2 + exon3 + exon4 + exon5 + exon6 + exon7 + exon8 + exon9 + exon10 + exon11
1172/3
114/(114+47)