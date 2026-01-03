## Summarise sequencing depth information
## For results section in ms. 

# depth file by gene/region, subject, and panel
out_path <- "/Users/huntc10/Library/Mobile Documents/com~apple~CloudDocs/GenomePhD/Dissertation/proj/li_fraumeni/manuscript/supplemental_data/"
out_file <- "supplemental_depth_panels.csv"
supp_file_depth <- read_delim(paste0(out_path,out_file))

## average sequencing depth across all samples (all regions)
## sorted from maximum to minimum
supp_file_depth %>%
  group_by(Subject, Tissue) %>%
  summarise(DP = mean(MeanDepth)) %>% 
  arrange(desc(DP)) %>%
  ungroup() %>%
  mutate(average = mean(DP)) %>%
  print(n=Inf)


## average sequencing depth across blood samples (all regions)
blood_tissues <- c("PBMC", "Buffy", "Plasma", "WB", "BM")

supp_file_depth %>%
  filter(Tissue %in% blood_tissues) %>%
  group_by(Subject, Tissue) %>%
  summarise(DP = mean(MeanDepth)) %>% 
  arrange(desc(DP)) %>%
  ungroup() %>%
  mutate(average = mean(DP)) %>%
  print(n=Inf)


## average sequencing depth across tissue samples (all regions)
supp_file_depth %>%
  filter(Subject == "LFS01") %>%
  filter(!(Tissue %in% blood_tissues)) %>%
  group_by(Subject, Tissue) %>%
  summarise(DP = mean(MeanDepth)) %>% 
  arrange(desc(DP)) %>%
  ungroup() %>%
  mutate(average = mean(DP)) %>%
  print(n=Inf)


#########
#### Total mutations coding/non-coding
#########


count_all_coding <- maf_masked_coding %>% 
  filter(!(Samp %in% samples_all_exclude)) %>%
  filter(!is.na(gene_name) & (!is.na(exon_number) | Variant_Classification == "Splice_Site")) %>%
  summarise(count_all_coding = n())
count_all_coding


maf_masked_coding %>% 
  filter(!(Samp %in% samples_all_exclude)) %>%
  filter(is.na(gene_name)) %>% print(width = Inf)
filter(!is.na(gene_name) & (!is.na(exon_number) | Variant_Classification == "Splice_Site")) %>%
  summarise(count_all_coding = n())
count_all_coding

count_all_non_coding <- maf_masked_coding %>%
  filter(!(Samp %in% samples_all_exclude)) %>%
  filter(is.na(gene_name) | (!is.na(gene_name) & is.na(exon_number))) %>% 
  filter(Variant_Classification == "Splice_Site") %>%
  print(width = Inf)

unique(maf_masked_coding$Tissue)
maf_masked_coding %>% filter(Tissue %in% c("Colon", "Thyroid", "Mainstem bronchus", "Skin, non-sun-exposed")) %>% print(n=Inf, width = Inf)


## mutation count in CHIP (and TP53) coding
maf_masked_coding %>% print(width = Inf)

CHIP_genes
count_coding_all <- maf_masked_coding %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(is.na(MUT_region_StartPosition) | is.na(MUT_region_EndPosition)) %>% 
  filter(!(Tissue %in% c("Urine cells"))) %>% 
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>% print(n = Inf)


unique(coding_check$Variant_Classification)
maf_masked_coding %>% filter(Variant_Classification == "Intron")  %>% print(width = Inf)
