## Li-Fraumeni skyscraper plot


### uw mstp volunteers
mstp <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")

skyscraper_prep <-
  filt_maf %>%
  filter(coding == "coding") %>%
  filter(Subject %in% mstp) %>%
  
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


skyscraper <- skyscraper_prep %>%
  ggplot(
    aes(
      x = Subject,
      y = SampCodingOrder,
      fill = color_group,
      label = t_alt_count
    )
  ) +
  geom_tile() +
  
  # uncomment
  #geom_text( size=font.subscript.size*25.4/72.27) +
  
  geom_tile(data = skyscraper_prep,
            aes(
              x = Subject,
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
ggsave("results/skyscraper_plot_mstp.png", skyscraper_slide, width = 15, height = 9, units = "in", dpi = 300)
