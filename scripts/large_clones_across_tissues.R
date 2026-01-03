
tissue_order <- c("Whole blood", 
                  "Buffy coat", 
                  "Plasma", 
                  "Bone marrow", 
                  "Buccal mucosa", 
                  "Thyroid", 
                  "Mainstem bronchus",
                  "Lung", 
                  "Esophagus 1", 
                  "Esophagus 2", 
                  "Gastric 1",
                  "Gastric 2",
                  "Cardiac muscle",
                  "Spleen",
                  "Liver",
                  "Colon",
                  "Omentum",
                  "Peritoneum",
                  "Renal",
                  "Testis",
                  "Skeletal muscle",
                  "Skin",
                  "Skin, non-sun-exposed",
                  "Mediastinal metastasis",
                  "Lung metastasis",
                  "Esophageal cancer 1",
                  "Esophageal cancer 2",
                  "Liver metastasis 1",
                  "Liver metastasis 2")

cancer_samples <- c(                  
                    "Mediastinal metastasis",
                    "Lung metastasis",
                    "Esophageal cancer 1",
                    "Esophageal cancer 2",
                    "Liver metastasis 1",
                    "Liver metastasis 2")

#### tumor vs not
small_clones <- maf_masked_coding %>%
  filter(Tissue %in% tissue_order) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(Subject == "Patient") %>% 
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>% 
  group_by(Tissue) %>% 
  summarise(count = n(), 
            burden = sum(t_alt_count)) %>% print()

chip_depths_tissue <- final_masked_depth %>% 
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  mutate(region_split = sub("_.*", "", Gene),
         region = ifelse(region_split == "region",
                         sub(".*_", "", Gene), NA)
  ) %>%
  filter((region_split == "region") | (!is.na(gene_name))) %>%
  mutate(plot_region = if_else(is.na(gene_name), region, gene_name)) %>%
  group_by(plot_region, Samp) %>%
  filter(gene_name == "TP53" & subject == "Patient") %>% 
  group_by(tissue) %>% 
  summarise(depth = sum(DP))

small_clones
chip_depths_tissue
small_clones_MF <- small_clones %>% left_join(chip_depths_tissue, by = c("Tissue" = "tissue")) %>%
  mutate(MF = count/depth,
         MB = burden/depth) %>%
  mutate(cancer = if_else(Tissue %in% cancer_samples, "cancer", "non-cancer")) %>%
  filter(Tissue != "Skin")
  
df_long <- small_clones_MF %>%
  pivot_longer(cols = c(MF, MB), names_to = "Metric", values_to = "Value")
ggplot(df_long, aes(x = cancer, y = Value)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )

wilcox_MF <- wilcox.test(MF ~ cancer, data = small_clones_MF)
wilcox_MF

wilcox_MB <- wilcox.test(MB ~ cancer, data = small_clones_MF)
wilcox_MB


### small clones
small_clones <- maf_masked_coding %>%
  filter(Tissue %in% tissue_order) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(Subject == "Patient") %>% 
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>% 
  filter(t_alt_count == 1) %>% 
  group_by(Tissue) %>% 
  summarise(count = n()) %>% print()
quantiles <- quantile(small_clones$count, probs = c(0.25, 0.75))
iqr_range <- c(lower = quantiles[1], upper = quantiles[2])
iqr_range
