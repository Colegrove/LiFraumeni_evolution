## Lollipop plot

# default plot size
source("scripts/01.1.05_makeLolipopComparison.R")
# plot for slide
source("scripts/makeLolipopComparison_slide.R")

# TP53 domains
filteredDomains <- read_delim(inputs$domains_file)

tissue = "All tissue-types"
#tissue = "Skin"
non_cancer_samples = c("All tissue-types",
                       "Whole blood", 
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

all_samples = c("All tissue-types",
                "Whole blood", 
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

tissue_filtering <- filt_maf %>% 
  filter(!is.na(prot.pos)) %>% 
  filter(Subject == "Patient") %>%
  filter(Tissue %in% all_samples) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(coding_from_maf == "coding") %>%
  filter(Variant_Classification != "Splice_Site") %>%
  filter(inMask == FALSE) %>% 
  filter(Variant_Type=="SNP")

tissue_filtering %>%  print(width = Inf)
print(unique(tissue_filtering$Tissue))

#### Slide size

tissue <- "All tissues"
lollipop_src <- lolipopComparison_slide(
  tissue_filtering, 
  tissue,
  COSMIC_table_forLolipop, 
  filteredDomains, 
  data_name = "Sample", 
  comp_data_name = "COSMIC")



# if tissue specific
lollipop_plot <- lollipop_src[["TP53"]]
# if not tissue specific
lollipop_plot_peices <- lollipop_src[["TP53"]]

slide_theme <- theme(
  axis.text.x = element_text(size=8),
  axis.text.y = element_text(size = 8),
  axis.title.y = element_text(size=8)
)


test_plot <- lollipop_plot_peices[["test_plot"]] + slide_theme + theme(
  plot.margin = margin(2,2,2,2)
)
test_domains_plot <- lollipop_plot_peices[["test_domains_plot"]]
test_plot_2 <- lollipop_plot_peices[["test_plot_2"]] + slide_theme

domain_legend <- get_legend(test_domains_plot +
    theme(legend.key.size = unit(8, "pt"),
          legend.text = element_text(size = 8, margin= margin(0,0,0,1)), 
          legend.box.margin = margin(0, 0, 0, 0),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.key.spacing.x = unit(3, "pt"))
    )

#domain_legend <- get_legend(test_domains_plot)

lollipop_slide <- cowplot::plot_grid(
  test_plot,
  test_domains_plot + theme(legend.position="none"),
  test_plot_2,
  domain_legend,
  align='v', axis="lr", rel_heights = c(3.25,1,2,1.25), ncol = 1) +
  theme(plot.margin = margin(b=-5))

show(lollipop_slide)
ggsave("results/lollipop_tissues_ms_4B.png", lollipop_slide, width = 3, height = 1.5, units = "in", dpi = 300)


