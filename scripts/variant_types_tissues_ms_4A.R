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
 

################################################################################
######### Mutation Types
################################################################################


tissue_abbreviations <- tribble(
  ~Tissue,                   ~Tissue_abbr,
  "Whole blood",             "WB",
  "Buffy coat",              "Buffy",
  "Plasma",                  "Plasma",
  "Bone marrow",             "BM",
  #"Buccal mucosa",           "Bucc",
  "Thyroid",                 "Thyroid",
  "Mainstem bronchus",       "Bronchus",
  "Lung",                    "Lung",
  "Esophagus 1",             "Esoph1",
  "Esophagus 2",             "Esoph2",
  "Gastric 1",               "Gast1",
  "Gastric 2",               "Gast2",
  "Cardiac muscle",          "Cardiac",
  "Spleen",                  "Spleen",
  "Liver",                   "Liver",
  "Colon",                   "Colon",
  "Omentum",                 "Omentum",
  "Peritoneum",              "Peritoneum",
  "Renal",                   "Renal",
  "Testis",                  "Testis",
  "Skeletal muscle",         "Skeletal",
  "Skin",                    "Skin",
  "Skin, non-sun-exposed",   "SkinNS",
  "Mediastinal metastasis",  "Med Met",
  "Lung metastasis",         "Lung Met",
  "Esophageal cancer 1",     "Esoph Ca1",
  "Esophageal cancer 2",     "Esoph Ca2",
  "Liver metastasis 1",      "Liver Met1",
  "Liver metastasis 2",      "Liver Met2",
  "PBMC",                    "PBMC"
)

custom_label <- function(x) {
  sapply(x, function(t) {
    abbr <- tissue_abbreviations$Tissue_abbr[match(t, tissue_abbreviations$Tissue)]
    if (is.na(abbr)) abbr <- t
    if (t %in% cancer_samples) {
      paste0("<span style='color:red;'>", abbr, "</span>")
    } else {
      abbr
    }
  })
}

mutation_type_prep <- skyscraper_prep %>%
  mutate(Mutation_type = factor(Mutation_type, levels = c(
    "Silent",
    "Splice",
    "Indel",
    "Nonsense_Mutation",
    "Missense_Mutation"
  )))

# variant_colors <- c(
#   "Indel"             = "#E69F00",
#   "Missense_Mutation" = "#F0E442", 
#   "Nonsense_Mutation" = "#009E73", 
#   "Silent"            = "#56B4E9", 
#   "Splice"            = "#CC79A7"  
# )
variant_colors <- c(
  "Indel"             = "#CC6677",
  "Missense_Mutation" = "#DDCC77", 
  "Nonsense_Mutation" = "#117733", 
  "Silent"            = "#88CCEE", 
  "Splice"            = "#999933"  
)

variant_type_proportion <- ggplot(mutation_type_prep, aes(x = Tissue_ordered, fill = Mutation_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = variant_colors,     
    labels = c(
    "Indel"             = "Indel",
    "Missense_Mutation" = "Missense",
    "Nonsense_Mutation" = "Nonsense",
    "Silent"            = "Silent",
    "Splice"            = "Splice"
  )
  ) +
  labs(
    x = "Tissue",
    #y = "Proportion",
    y = "Prop.",
    fill = "Variant Classification"
  ) +
  theme_minimal() +
  scale_x_discrete(labels = custom_label) +
  scale_y_continuous(breaks = c(0, 1), labels = c("0","1")) +
  theme(#axis.text.x = element_blank(),
        axis.text.x.bottom = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size=8),
        
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y = element_markdown(size=8, hjust=0.5),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        #legend.position = "none",
        legend.text = element_text(size=8, margin=margin(r=2)),
        legend.key.size = unit(8,"pt"),
        legend.key.spacing.x = unit(3,"pt"),
        legend.margin = margin(-15,0,0,0),
        axis.text.y.left = element_text(size=8, margin = margin(r=2)),
        plot.margin=margin(2,4,2,4))

variant_type_proportion
#ggsave("results/variant_type_proportion_tissues.png", variant_type_proportion, width = 4, height = .65, units= "in", dpi = 300)
#ggsave("results/variant_type_proportion_tissues.png", variant_type_proportion, width = 4, height = .5, units= "in", dpi = 300)
#ggsave("results/variant_type_proportion_tissues.png", variant_type_proportion, width = 4, height = .5, units= "in", dpi = 300)

