## Li-Fraumeni skyscraper plot


### family member samples
family <- c("Family member A", "Family member B", "Family member C")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")
skyscraper_prep <-
  filt_maf %>%
  filter(Tissue != "Urine cells") %>%
  filter(coding == "coding") %>%
  #filter(Subject %in% family) %>%
  filter(Subject %in% family_patient_blood) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  
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

skyscraper_prep2 <- skyscraper_prep %>%
  group_by(Subject) %>%
  mutate(n_tissues_for_subject = n_distinct(Tissue)) %>%
  ungroup() %>%
  mutate(
    x_axis_var = if_else(n_tissues_for_subject > 1, Tissue, Subject)
  )
skyscraper_prep2 <- skyscraper_prep2 %>%
  mutate(x_axis_var = factor(
    x_axis_var,
    levels = c("Family member A", "Family member B", "Family member C", "Buffy coat", "Plasma", "Whole blood")
  ))


skyscraper <- skyscraper_prep2 %>%
  ggplot(
    aes(
      x = x_axis_var,
      y = SampCodingOrder,
      fill = color_group,
      label = t_alt_count
    )
  ) +
  geom_tile() +
  
  # uncomment
  geom_text( size=font.subscript.size*25.4/72.27) +
  
  geom_tile(data = skyscraper_prep2,
            aes(
              x = x_axis_var,
              y = SampCodingOrder
            ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
  scale_y_continuous(expand = c(0,0),
                     breaks=c(0.5,5.5,10.5,15.5,20.5,25.5,30.5,35.5,40.5,45.5,50.5,55.5,60.5),
                     labels=c(0,5,10,15, 20, 25, 30,35,40,45,50,55,60),
                     limits = c(0,25),
                     name = "Number of <i>TP53</i> mutations") +
  Col.amClass.fill +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size = 8),
        panel.spacing = unit(0,"pt"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        legend.position = c(0.25,0.7),
  )

show(skyscraper)

skyscraper_slide <- skyscraper +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5,size = 20),
    axis.text.y = element_text(size = 20),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7),
    legend.text = (element_text(size = 20)),
    legend.title = (element_text(size=20))
  ) + geom_text(size=16*25.4/72.27)

show(skyscraper_slide)
ggsave("results/skyscraper_plot_family.png", skyscraper_slide, width = 17, height = 13, units = "in", dpi = 300)
