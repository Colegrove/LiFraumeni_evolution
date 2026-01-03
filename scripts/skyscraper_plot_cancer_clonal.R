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
MAF_table
high_freq_tp53 <- MAF_table %>% filter(VAF >= 0.01) %>% filter(Subject == "Patient") %>% filter(Hugo_Symbol == "TP53")
high_freq_tp53
duplicates_summary <- high_freq_tp53 %>%
  group_by(Start_Position, Reference_Allele, Tumor_Seq_Allele2, HGVSp_Short) %>%
  summarize(
    duplicate_count = n(),
    tissues         = toString(unique(Tissue)),  # or paste(unique(Tissue), collapse = ", ")
    .groups         = "drop"
  )

duplicates_summary
duplicates_summary_1 <- duplicates_summary %>% filter(Start_Position == 7676301)

cat(duplicates_summary_1$tissues)


# 
# 
# skyscraper_prep
# skyscraper_prep <-
#   filt_maf %>%
#   filter(Tissue %in% cancer_samples) %>%
#   filter(coding == "coding") %>%
#   filter(Subject == "Patient") %>%
# 
#   mutate(
#     color_group = case_when(
#       am_class == "likely_benign" & t_alt_count == 1 ~ "likely_benign",
#       am_class == "likely_benign" & t_alt_count > 1 ~ "likely_benign_LC",
#       am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
#       am_class == "likely_pathogenic" & t_alt_count > 1 ~ "likely_pathogenic_LC",
#       am_class == "ambiguous" & t_alt_count == 1 ~ "ambiguous",
#       am_class == "ambiguous" & t_alt_count > 1 ~ "ambiguous_LC"
#     )
#   ) 
# print(skyscraper_prep, n=Inf)
# 
# #skyscraper_prep )
# skyscraper_prep %>% 
#   filter(Tissue == "Mediastinal metastasis") %>% 
#   #select(Start_Position, Tumor_Sample_Barcode, am_class, t_alt_count, 
#          #Mutation_Class, Mutation_type, HGVSp) %>% 
#   #filter(Mutation_Class == "SNP") %>% 
#   #filter(is.na(color_group)) %>%
#   print(n=Inf)
# 
# 
# skyscraper <- skyscraper_prep %>%
#   ggplot(
#     aes(
#       x = Tissue,
#       y = SampCodingOrder,
#       fill = color_group,
#       label = t_alt_count
#     )
#   ) +
#   geom_tile() +
#   
#   # uncomment
#   geom_text( size=font.subscript.size*25.4/72.27) +
#   
#   geom_tile(data = skyscraper_prep,
#             aes(
#               x = Tissue,
#               y = SampCodingOrder
#             ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
#   scale_y_continuous(expand = c(0,0),
#                      breaks=c(0.5,5.5,10.5,15.5,20.5,25.5,30.5,35.5,40.5,45.5,50.5,55.5,60.5),
#                      labels=c(0,5,10,15, 20, 25, 30,35,40,45,50,55,60),
#                      limits = c(0,20),
#                      name = "Number of <i>TP53</i> mutations") +
#   Col.amClass.fill +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size = 8),
#         panel.spacing = unit(0,"pt"),
#         strip.background = element_rect(fill = "white", color = "black"),
#         strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
#         axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
#         axis.title.x = element_blank(),
#         legend.position = c(0.25,0.7),
#   )
# 
# show(skyscraper)
# 
# skyscraper_slide <- skyscraper +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5,size = 20),
#     axis.text.y = element_text(size = 20),
#     panel.spacing = unit(0,"pt"),
#     strip.background = element_rect(fill = "white", color = "black"),
#     strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
#     axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
#     axis.title.x = element_blank(),
#     legend.position = c(0.25,0.7),
#     legend.text = (element_text(size = 20)),
#     legend.title = (element_text(size=20))
#   ) + geom_text(size=16*25.4/72.27)
# 
# show(skyscraper_slide)
# ggsave("results/skyscraper_plot_slide_patient_cancer.png", skyscraper_slide, width = 17, height = 13, units = "in", dpi = 300)
