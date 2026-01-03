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
print(skyscraper_prep, n=Inf) %>% arrange(desc(VAF)) %>% select(VAF)

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


slide = TRUE
if(slide == FALSE){

  skyscraper <- skyscraper_prep %>%
    ggplot(
      aes(
        x = Tissue,
        y = SampCodingOrder,
        fill = color_group,
        label = t_alt_count
      )
    ) +
    geom_tile() +
    
    # uncomment
    #geom_text( size=font.subscript.size*25.4/72.27) +
    geom_text( size=9*25.4/72.27) +
    
    geom_tile(data = skyscraper_prep,
              aes(
                x = Tissue,
                y = SampCodingOrder
              ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
    scale_y_continuous(expand = c(0,0),
                       #breaks=c(0.5,5.5,10.5,15.5,20.5,25.5,30.5,35.5,40.5,45.5,50.5,55.5,60.5),
                       #labels=c(0,5,10,15, 20, 25, 30,35,40,45,50,55,60),
                       limits = c(0,55),
                       name = "Number of <i>TP53</i> mutations") +
    Col.amClass.fill +
    theme_bw() + 
    scale_x_discrete(labels = custom_label) +
    theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 8),
          panel.spacing = unit(0,"pt"),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
          axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
          axis.title.x = element_blank(),
          legend.position = c(0.25,0.7),
    )
  show(skyscraper)
  ggsave("results/skyscraper_plot_patient_cancer.png", skyscraper, width = 7, height = 6, units = "in", dpi = 300)
  
} else if(slide == TRUE){
  skyscraper <- skyscraper_prep %>%
    ggplot(aes(x = Tissue,y = SampCodingOrder,fill = color_group,label = t_alt_count)
    ) +
    geom_tile(width = 1, height = 1) +
    
    # uncomment
    #geom_text( size=font.subscript.size*25.4/72.27) +
    #geom_text( size=9*25.4/72.27) +
    
    geom_tile(data = skyscraper_prep,
              aes(
                x = Tissue,
                y = SampCodingOrder, 
                height = 1, 
                width = 1
              ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
    scale_y_continuous(expand = c(0,0),
                       #breaks=c(0.5,5.5,10.5,15.5,20.5,25.5,30.5,35.5,40.5,45.5,50.5,55.5,60.5),
                       #labels=c(0,5,10,15, 20, 25, 30,35,40,45,50,55,60),
                       limits = c(0,55),
                       name = "Number of <i>TP53</i> mutations") +
    Col.amClass.fill +
    theme_bw() + 
    scale_x_discrete(labels = custom_label, expand = c(0,0)) +
    theme(#axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 8),
          panel.spacing = unit(0,"pt"),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
          axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
          axis.title.x = element_blank(),
          legend.position = c(0.25,0.7),
    )
  
  
  skyscraper_slide <- skyscraper +

    theme(
      axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 20),
      axis.text.y = element_text(size = 20),
      panel.spacing = unit(0,"pt"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
      axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
      axis.title.x = element_blank(),
      legend.position = c(0.25,0.7),
      legend.text = (element_text(size=32)),
      legend.title = (element_text(size=32))
    )

  show(skyscraper_slide)
  ggsave("results/skyscraper_plot_slide_patient_cancer_noCount.png", skyscraper_slide, width = 22, height = 12, units = "in", dpi = 300)

}


skyscraper_prep %>% filter(is.na(am_pathogenicity)) %>% print(width = Inf)
skyscraper_prep
ggplot(skyscraper_prep, aes(x = VAF)) +
  geom_histogram() + 
  scale_x_log10()

#####################################################################
########### Research Reports subsection
#####################################################################


##############################
############### number of mutations
##############################

skyscraper_nMuts <- skyscraper_prep %>%
  ggplot(aes(x = Tissue, y=SampCodingOrder)) +
  geom_tile(width = 1, height = 1) +
  geom_tile(data = skyscraper_prep,
            aes(x = Tissue,y = SampCodingOrder, height = 1, width = 1), 
            inherit.aes = F, color = "black", fill = "#CCCCCC", linewidth = .1) + 
  #geom_bar(fill = 'black') + 
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,55),
                     name = "Number of <i>TP53</i> mutations") +
  theme_bw() + 
  scale_x_discrete(labels = custom_label, expand = c(0,0)) +
  theme(
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7),
  )


skyscraper_nMuts_slide <- skyscraper_nMuts +
  
  theme(
    axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 20),
    axis.text.y = element_text(size = 20),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7),
    legend.text = (element_text(size = 20)),
    legend.title = (element_text(size=20))
  )

show(skyscraper_nMuts_slide)
ggsave("results/skyscraper_nMuts_slide.png", skyscraper_nMuts_slide, width = 22, height = 12, units = "in", dpi = 300)


##############################
############### size of mutations
##############################

Col.cloneSize.breaks <- c("small", "large")
Col.cloneSize.values <- c("#CCCCCC","#555555")
Col.cloneSize.labels <- c("1",">1")
Col.cloneSize.fill <- scale_fill_manual(
  breaks = Col.cloneSize.breaks,
  values = Col.cloneSize.values,
  labels = Col.cloneSize.labels,
  na.value = NA,
  name = "Duplex reads"
)

skyscraper_prep_sizeMuts <- skyscraper_prep %>%
  mutate(cloneSize = case_when(t_alt_count == 1 ~ "small", 
                               t_alt_count > 1 ~ "large"))
#"#B0B0B0"
skyscraper_sizeMuts <- skyscraper_prep_sizeMuts %>%
  ggplot(aes(x = Tissue, y=SampCodingOrder, fill=cloneSize)) +
  geom_tile(width = 1, height = 1) +
  geom_tile(data = skyscraper_prep,
            aes(x = Tissue,y = SampCodingOrder, height = 1, width = 1), 
            inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
  #geom_bar(fill = 'black') + 
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,55),
                     name = "Number of <i>TP53</i> mutations") +
  Col.cloneSize.fill +
  theme_bw() + 
  scale_x_discrete(labels = custom_label, expand = c(0,0)) +
  theme(
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7),
  )


skyscraper_sizeMuts_slide <- skyscraper_sizeMuts +
  
  theme(
    axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 20),
    axis.text.y = element_text(size = 20),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7),
    legend.text = (element_text(size = 32)),
    legend.title = (element_text(size=32))
  )

show(skyscraper_sizeMuts_slide)
ggsave("results/skyscraper_sizeMuts_slide.png", skyscraper_sizeMuts_slide, width = 22, height = 12, units = "in", dpi = 300)

