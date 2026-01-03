install.packages("fuzzyjoin")
library("fuzzyjoin")

SNP_list <- MAF_table %>% 
  mutate(isHighVaf = (VAF>= 0.9 | (VAF >= 0.4 & VAF <= 0.6))) %>%
  select(Tumor_Sample_Barcode, 
         Chromosome, 
         Start_Position, 
         End_Position, 
         Reference_Allele, Tumor_Seq_Allele2, 
         VAF, 
         isHighVaf) %>% 
  filter(isHighVaf == TRUE) %>% 
  group_by(Tumor_Sample_Barcode) %>%
  select(-c(isHighVaf))

maf_for_exam_SNP <- MAF_table %>%
  group_by(Tumor_Sample_Barcode) %>%
  select(Tumor_Sample_Barcode, 
         Chromosome, 
         Start_Position, 
         End_Position, 
         VAF, 
         FILTER, 
         Reference_Allele, Tumor_Seq_Allele2, 
         t_depth) %>% 
  mutate(isHighVaf = (VAF>= 0.9 | (VAF >= 0.4 & VAF <= 0.6))) %>%
  genome_left_join(
    SNP_list, 
    by=c("Tumor_Sample_Barcode","Start_Position", 
         "End_Position")) %>% 
  filter(!is.na(Tumor_Sample_Barcode.y)) %>% 
  mutate(joined_length = End_Position.y - Start_Position.y + 1,
         base_length = End_Position.x - Start_Position.x + 1) %>% 
  filter(!isHighVaf) %>% 
  filter(!grepl("masked", FILTER)) 
maf_for_exam_clusterd <- MAF_table %>%
  group_by(Tumor_Sample_Barcode) %>%
  select(Tumor_Sample_Barcode.x = Tumor_Sample_Barcode, 
         Chromosome.x = Chromosome, 
         Start_Position.x = Start_Position, 
         End_Position.x = End_Position, 
         VAF.x = VAF, 
         FILTER, 
         Reference_Allele.x = Reference_Allele, 
         Tumor_Seq_Allele2.x = Tumor_Seq_Allele2, 
         t_depth) %>% 
  mutate(isHighVaf = (VAF.x>= 0.9 | (VAF.x >= 0.4 & VAF.x <= 0.6))) %>%
  filter(!isHighVaf) %>% 
  filter(!grepl("masked", FILTER)) %>% 
  filter(grepl("clustered", FILTER))

maf_for_exam_SNP %>% 
  rbind(maf_for_exam_clusterd) %>% 
  write_delim(paste0(out_prefix(), ".Vars_for_manual_examination.tsv"), delim = "\t")
