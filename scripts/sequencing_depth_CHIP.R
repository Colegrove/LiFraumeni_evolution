### Visualize sequencing depths across Persons and genes in CHIP panel


sample_id_mapping_path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sample_map <- read_delim(sample_id_mapping_path, delim = "\t", quote="\"") %>%
  mutate(tissue = str_trim(str_replace_all(tissue, '"', '')), 
         subject = str_trim(str_remove(subject, ":$")))
sample_map %>% print(n=Inf)
# all blood
# samples_with_chip <- c("DNA03977_S14.1", "DNA03975_S12.1", "DNA03976_S13.1", 
#                        "DNA03978_S15.1", "DNA03962_S2.1", "DNA03963_S3.1", 
#                        "DNA03964_S4.1", "DNA03968_S5.1", 
#                        "DNA03969_S6.1", "DNA03970_S7.1", "DNA03971_S8.1", 
#                        "DNA03972_S9.1", "DNA03973_S10.1", "DNA03974_S11.1")
# buffy coat only
samples_with_chip <- c("DNA03975_S12.1", "DNA03962_S2.1", "DNA03963_S3.1", 
                       "DNA03964_S4.1", "DNA03968_S5.1", "DNA03969_S6.1", 
                       "DNA03970_S7.1", "DNA03971_S8.1", "DNA03972_S9.1", 
                       "DNA03973_S10.1", "DNA03974_S11.1")

chip_depths <- DP_filt_CHIP %>% 
  filter(InBed == TRUE & Samp %in% samples_with_chip) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  #filter(Samp == "DNA03975_S12.1") %>%
  mutate(full_gene = sub("_.*","",Gene)) %>%
  filter(full_gene != "region") %>%
  group_by(Samp, subject, tissue, full_gene) %>%
  summarise(mean_gene_DP = mean(DP), .groups = "drop") 
  #ungroup()
chip_depths %>% print(width = Inf)


chip_depths_sum <- group_by(Samp, ) %>% summarise(meanDP = sum(mean_gene_DP))

depth_per_subject <- ggplot(chip_depths, aes(x = full_gene, y = mean_gene_DP)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ subject, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "Gene", y = "mean depth")

depth_per_subject
ggsave("./results/depth_by_subject_chip.png", depth_per_subject, width = 12, height = 12)


depth_by_gene <- ggplot(chip_depths, aes(x = subject, y = mean_gene_DP)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ full_gene, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "Subject", y = "Mean depth (DP)")

ggsave("./results/depth_by_gene_chip.png", depth_by_gene, width = 12, height = 12)




filt_maf_CHIP %>% filter(Hugo_Symbol == "NRAS") %>% arrange(Start_Position) %>% print(width = Inf)
