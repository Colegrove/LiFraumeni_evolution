#### TP53 mutation number regression in LFS01 tissues

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

sample_id_mapping_path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sample_map <- read_delim(sample_id_mapping_path, delim = "\t", quote="\"") %>%
  mutate(tissue = str_trim(str_replace_all(tissue, '"', '')), 
         subject = str_trim(str_remove(subject, ":$")))

## number of coding mutations
tissue_counts <-
  maf_masked_coding %>%
  filter(Tissue %in% tissue_order) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(Subject == "Patient") %>% 
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>% 
  group_by(Tissue) %>% summarise(mut_count = n())

## tp53 coding depth by tissue
tissue_tp53_depths <- final_masked_depth %>% 
  filter(gene_name == "TP53") %>%
  filter(!is.na(exon_number)) %>%
  filter(!inRepeatMask) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  
  filter(tissue %in% tissue_order) %>%
  group_by(subject, tissue, gene_name) %>%
  summarise(denominator_coding = sum(DP), .groups = "drop") 

## combine mutations and depth
tissue_counts
combine_df_all <- tissue_counts %>% left_join(tissue_tp53_depths, by = c("Tissue" = "tissue"))
combine_df_no_skin <- combine_df %>% filter(Tissue != "Skin")

## linear regression
lm_all <- lm(mut_count ~ denominator_coding, data = combine_df_all)
lm_no_skin <- lm(mut_count ~ denominator_coding, data = combine_df_no_skin)

r2_all <- summary(lm_all)$r.squared
p_all <- coef(summary(lm_all))[2, 4]

r2_no_skin <- summary(lm_no_skin)$r.squared
p_no_skin <- coef(summary(lm_no_skin))[2, 4]

plot_all <- ggplot(combine_df_all, aes(x = denominator_coding, y = mut_count)) +
  geom_point(size = 2, color = "#3182bd", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R² = ", round(r2_all, 3), "\n",
                          "p = ", signif(p_all, 3)),
           size = 2.8) +
  labs(
    x = "Sequencing depth",
    y = "Mutation count",
    title = "All tissues"
  ) +
  theme_classic(base_size = 8)

plot_no_skin <- ggplot(combine_df_no_skin, aes(x = denominator_coding, y = mut_count)) +
  geom_point(size = 2, color = "#3182bd", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R² = ", round(r2_no_skin, 3), "\n",
                          "p = ", signif(p_no_skin, 3)),
           size = 2.8) +
  labs(
    x = "Sequencing depth",
    y = "Mutation count",
    title = "Excluding sun-exposed skin"
  ) +
  theme_classic(base_size = 8)

combined_plot <- plot_all + plot_no_skin
combined_plot
ggsave("./results/tissues_depth_regression.png", combined_plot, height = 2.5, width = 4, units = "in", dpi = 300)


################################################################################
#################### MF vs sequencing depth
################################################################################

## combine mutations and depth
tissue_counts
combine_df_all <- tissue_counts %>% left_join(tissue_tp53_depths, by = c("Tissue" = "tissue"))
combine_df_no_skin <- combine_df %>% filter(Tissue != "Skin")

## linear regression
lm_all <- lm(mut_count/denominator_coding ~ denominator_coding, data = combine_df_all)
lm_no_skin <- lm(mut_count/denominator_coding ~ denominator_coding, data = combine_df_no_skin)

r2_all <- summary(lm_all)$r.squared
p_all <- coef(summary(lm_all))[2, 4]

r2_no_skin <- summary(lm_no_skin)$r.squared
p_no_skin <- coef(summary(lm_no_skin))[2, 4]

plot_all <- ggplot(combine_df_all, aes(x = denominator_coding, y = mut_count/denominator_coding)) +
  geom_point(size = 2, color = "#3182bd", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R² = ", round(r2_all, 3), "\n",
                          "p = ", signif(p_all, 3)),
           size = 2.8) +
  labs(
    x = "Sequencing depth",
    y = "Mutation frequency",
    title = "All tissues"
  ) +
  theme_classic(base_size = 8)

plot_no_skin <- ggplot(combine_df_no_skin, aes(x = denominator_coding, y = mut_count/denominator_coding)) +
  geom_point(size = 2, color = "#3182bd", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R² = ", round(r2_no_skin, 3), "\n",
                          "p = ", signif(p_no_skin, 3)),
           size = 2.8) +
  labs(
    x = "Sequencing depth",
    y = "Mutation frequency",
    title = "Excluding sun-exposed skin"
  ) +
  theme_classic(base_size = 8)

combined_plot_MF <- plot_all + plot_no_skin
combined_plot_MF
ggsave("./results/tissues_depth_regression_MF.png", combined_plot_MF, height = 2.5, width = 4, units = "in", dpi = 300)

