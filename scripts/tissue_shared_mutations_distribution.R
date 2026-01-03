## Li-Fraumeni shared tissue mutation distribution

non_cancer_samples = c("Whole blood", 
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
                       "Skin, non-sun-exposed")

cancer_samples = c("All tissue-types",
                   "Mediastinal metastasis",
                   "Lung metastasis",
                   "Esophageal cancer 1",
                   "Esophageal cancer 2",
                   "Liver metastasis 1",
                   "Liver metastasis 2")


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
 

skyscraper_prep <-
  filt_maf %>%
  filter(Tissue %in% tissue_order) %>%
  #filter(Tissue %in% cancer_samples | Tissue %in% non_cancer_samples) %>%
  filter(coding == "coding") %>%
  filter(Subject == "Patient") %>%
  
  # remove alternate splice site
  filter(Variant_Classification != "Intron") %>%

  mutate(
    color_group = case_when(
      am_class == "likely_benign" & t_alt_count == 1 ~ "likely_benign",
      am_class == "likely_benign" & t_alt_count > 1 ~ "likely_benign_LC",
      am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
      am_class == "likely_pathogenic" & t_alt_count > 1 ~ "likely_pathogenic_LC",
      am_class == "ambiguous" & t_alt_count == 1 ~ "ambiguous",
      am_class == "ambiguous" & t_alt_count > 1 ~ "ambiguous_LC"
    )
  ) 
print(skyscraper_prep, n=Inf) %>% arrange(desc(VAF)) %>% dplyr::select(VAF)

skyscraper_prep <- skyscraper_prep %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))

###################################
###### Shared mutations distribution
###################################

skyscraper_prep_shared <- skyscraper_prep %>% 
  group_by(Start_Position, Tumor_Seq_Allele2) %>%
  mutate(shared_muts = if (n() > 1) as.character(cur_group_id()) else "999") %>%
  mutate(tissues_shared = n()) %>%
  ungroup() 
skyscraper_prep_shared %>% print(width = Inf)
mut_idxs <- unique(skyscraper_prep_shared$shared_muts)

skyscraper_prep_shared <- skyscraper_prep_shared %>%
  group_by(Tissue, Subject) %>%
  arrange(desc(tissues_shared), SampCodingOrder, .by_group = TRUE) %>%
  mutate(NewOrder = row_number()) %>%
  ungroup()

shared_dist_prep <- skyscraper_prep_shared %>% distinct(Start_Position, .keep_all = TRUE)

## barplot of shared mutations
tissues_shared_plot <- ggplot(shared_dist_prep, aes(x = tissues_shared)) + 
  geom_bar(fill = "steelblue") +
  labs(
    x = "Tissues containing mutation",
    y = "Count"
  ) +
  theme_minimal()
tissues_shared_plot
ggsave("results/tissues_shared_distribution.png", tissues_shared_plot, width = 3, height = 2, units = "in", dpi = 300)
