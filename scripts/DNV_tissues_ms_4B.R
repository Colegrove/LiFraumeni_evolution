## Li-Fraumeni skyscraper plot

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
  maf_masked_coding %>%
  filter(Tissue %in% tissue_order) %>%
  #filter(Tissue %in% cancer_samples | Tissue %in% non_cancer_samples) %>%
  filter(gene_name == "TP53" & !is.na(exon_number)) %>%
  filter(Subject == "Patient") %>%

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

skyscraper_prep <- skyscraper_prep %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))

custom_label <- function(x) {
  sapply(x, function(t) {
    if (t %in% c("Mediastinal metastasis",
                 "Lung metastasis",
                 "Esophageal cancer 1",
                 "Esophageal cancer 2",
                 "Liver metastasis 1",
                 "Liver metastasis 2")){
      paste0("<span style='color:red;'>", t, "</span>")
    } else {
      t
    }
  })
}

################################################################################
######### Mutation Types
################################################################################
skyscraper_prep
mutation_type_prep <- skyscraper_prep %>%
  mutate(Mutation_type = factor(Mutation_type, levels = c(
    "Silent",
    "Splice",
    "Indel",
    "Nonsense_Mutation",
    "Missense_Mutation"
  )))


dnv <- mutation_type_prep %>% filter(Variant_Type == "DNP") %>%
  mutate(DNV_type = paste0(Reference_Allele,">",Tumor_Seq_Allele2)) %>%
  mutate(
    DNV_group = case_when(
      DNV_type %in% c("CC>TT", "GG>AA") ~ "CC>TT",
      DNV_type %in% c("GT>AA", "AC>TT") ~ "GT>AA",
      DNV_type %in% c("CT>AC", "GT>TG") ~ "CT>AC",
      DNV_type %in% c("GC>TT", "AA>GC") ~ "GC>TT",
      TRUE ~ DNV_type
    )
  ) %>%
  mutate(DNV_group = factor(DNV_group, levels = c("CT>AC", "CC>TT", "GC>TT", "GT>AA")))

dnv_colors <- c(
  "CC>TT" = "#E69F00",
  "CT>AC" = "#5F0873", 
  "GC>TT" = "#009E73",  
  "GT>AA" = "#56B4E9" 
)

dnv_colors <- c(
  "CC>TT" = "#000000",
  "CT>AC" = "#D9D9D9", 
  "GC>TT" = "#5d5d5d",  
  "GT>AA" = "#999999" 
)

dnv$DNV_group <- factor(dnv$DNV_group, 
                        levels = c("CT>AC", "GC>TT", "GT>AA", "CC>TT"))

dnv_counts <- ggplot(dnv, aes(x = Tissue, fill = DNV_group)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = dnv_colors) +
  
  #guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  
  ylim(0,15) +
  labs(
    x = "Tissue",
    y = "Count",
    fill = "DNV type"
  ) +
  theme_minimal() +
  theme(axis.text = element_text(size=8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=8),
        legend.text = element_text(size=8, margin=margin(0,0,0,0)),
        legend.title = element_blank(),
        legend.key.size=unit(8,"pt"),
        legend.margin = margin(0,0,0,0),
        
        #legend.position = c(.3,.7)
        legend.position = c(.5, 1.05)
        )


#ggsave("results/dnv_counts_tissues.png", dnv_counts, width = 3.75, height = 0.65, units= "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_4/dnv_counts_tissues.png", dnv_counts, width = 3.75, height = 0.65, units= "in", dpi = 300)


