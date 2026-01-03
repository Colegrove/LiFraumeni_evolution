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

filt_maf %>% filter(Start_Position == 7673535) %>% print(width=Inf)

tissue_filtering <- filt_maf %>% 
  filter(!is.na(prot.pos)) %>% 
  filter(Subject == "Patient") %>%
  #filter(Tissue %in% cancer_samples) %>%
  filter(Tissue %in% all_samples) %>%
  filter(coding == "coding") %>% 
  filter(inMask == FALSE) %>% 
  filter(Variant_Type=="SNP")

print(unique(tissue_filtering$Tissue))
#tissue_filtering
#tissue_filtering %>% filter(tissue == "Skin") %>% select("HGVSc", "HGVSp", "prot.ref", "prot.pos", "prot.alt") %>% arrange(prot.pos) %>% print(n=Inf)
# if(tissue != "All tissue-types"){
#   tissue_filtering <- tissue_filtering %>%
#     filter(Tissue == tissue) 
# }

# tissue_filtering
# tissue_filtering %>% 
#   filter(prot.pos == 248) %>%
#   dplyr::select(Hugo_Symbol, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, HGVSp, Tissue, t_alt_count)
# tissue_filtering %>% 
#   filter(prot.pos == 220)
# 
# tissue_filtering %>%
#   group_by(prot.pos) %>%
#   summarise(freq = n()) %>%
#   arrange(desc(freq))

## uncomment to plot each tissue individually

# all_tissues_plot <- function(tissue_df){
#   lollipop_src <- lolipopComparison(
#     #tissue_filtering, 
#     tissue_df,
#     tissue,
#     COSMIC_table_forLolipop, 
#     filteredDomains, 
#     data_name = "Sample", 
#     comp_data_name = "COSMIC")
#   
#   # if tissue specific
#   lollipop_plot <- lollipop_src[["TP53"]]
#   #show(lollipop_plot[["aligned"]])
#   # if not tissue specific
#   lollipop_plot_peices <- lollipop_src[["TP53"]]
#   
#   slide_theme <- theme(
#     axis.text.x = element_text(size=10),
#     axis.text.y = element_text(size = 10),
#     axis.title.y = element_text(size=10)
#   )
#   
#   Gene <- "TP53"
#   plot_title <- ggdraw() +
#     draw_label(
#       glue("{Gene} - {tissue}"),
#       #GeneIter,
#       fontface = 'bold.italic',
#       x = 0,
#       hjust = 0,
#       size = 12
#     ) +
#     theme(
#       # add margin on the left of the drawing canvas,
#       # so title is aligned with left edge of first plot
#       plot.margin = margin(0, 0, 0, 7),
#     )
#   
#   test_plot <- lollipop_plot_peices[["test_plot"]] + slide_theme
#   test_domains_plot <- lollipop_plot_peices[["test_domains_plot"]]
#   test_plot_2 <- lollipop_plot_peices[["test_plot_2"]] + slide_theme
#   domain_legend <- lollipop_plot_peices[["domain_legend"]] + slide_theme
#   
#   
#   domain_legend <- get_legend(
#     test_domains_plot +
#       theme(legend.box.margin = margin(7, 7, 7, 7),
#             legend.direction = "horizontal",
#             legend.background = element_blank(), 
#             legend.text = element_text(size = 10)) +
#       guides(
#         fill=guide_legend(
#           nrow=min(3, length(filteredDomains$Label))
#         )
#       )
#   )
#   
#   lollipop_slide <- cowplot::plot_grid(
#     plot_title,
#     test_plot,
#     test_domains_plot + theme(legend.position="none"),
#     test_plot_2,
#     domain_legend,
#     align='v', axis="lr", rel_heights = c(1,4,1,3,3), ncol = 1)
#   
#   show(lollipop_slide)
#   ggsave(glue("results/lollipop_plot_slide_{tissue}.png"), lollipop_slide, width = 6, height = 3, units = "in", dpi = 300)
# }
# 
# 
# for(tissue in cancer_samples){
#   if(tissue == "All tissue-types"){
#     all_tissues_plot(tissue_filtering)
#   }
#   else{
#     tissue_df <- tissue_filtering %>% 
#       filter(Tissue == tissue)
#     all_tissues_plot(tissue_df)
#   }
# }



#### Slide size

tissue_filtering %>% group_by(prot.pos) %>% summarise(count = n()) %>% arrange(desc(count))

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
#show(lollipop_plot[["aligned"]])
# if not tissue specific
lollipop_plot_peices <- lollipop_src[["TP53"]]

slide_theme <- theme(
  axis.text.x = element_text(size=10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size=10)
)

Gene <- "TP53"
plot_title <- ggdraw() +
  draw_label(
    glue("{Gene} - {tissue}"),
    #GeneIter,
    fontface = 'bold.italic',
    x = 0,
    hjust = 0,
    size = 12
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7),
  )

test_plot <- lollipop_plot_peices[["test_plot"]] + slide_theme
test_domains_plot <- lollipop_plot_peices[["test_domains_plot"]]
test_plot_2 <- lollipop_plot_peices[["test_plot_2"]] + slide_theme
#domain_legend <- lollipop_plot_peices[["domain_legend"]] + slide_theme
domain_legend <- get_legend(test_domains_plot +
    theme(legend.key.size = unit(1.5, "lines"),
          legend.text = element_text(size = 10), 
          #legend.title = element_text(size = 10),
          #legend.box.margin = margin(7, 7, 7, 7),
          legend.direction = "horizontal",
          legend.background = element_blank()
          )
    )

#domain_legend <- get_legend(test_domains_plot)

lollipop_slide <- cowplot::plot_grid(
  plot_title,
  test_plot,
  test_domains_plot + theme(legend.position="none"),
  test_plot_2,
  domain_legend,
  align='v', axis="lr", rel_heights = c(1,4,1,3,2), ncol = 1)

show(lollipop_plot_peices$aligned)
show(lollipop_slide)
ggsave("results/lollipop_plot_patient_cancer.png", lollipop_slide, width = 6, height = 3, units = "in", dpi = 300)

MAF_table_CHIP %>% filter(Hugo_Symbol == "TET2") %>% print(width = Inf)
  distinct(transcript_id) %>% print(width = Inf)
filt_maf_CHIP
