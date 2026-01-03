

vaf_comparison_prep <-
  maf_masked_coding %>%
  filter(Tissue %in% tissue_order) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(Subject == "Patient") %>%
  filter(!is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  mutate(
    color_group = case_when(
      am_class == "likely_benign" & t_alt_count == 1 ~ "likely_benign",
      am_class == "likely_benign" & t_alt_count > 1 ~ "likely_benign_LC",
      am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
      am_class == "likely_pathogenic" & t_alt_count > 1 ~ "likely_pathogenic_LC",
      am_class == "ambiguous" & t_alt_count == 1 ~ "ambiguous",
      am_class == "ambiguous" & t_alt_count > 1 ~ "ambiguous_LC"
    )
  ) %>%
  mutate(hotspot_color = if_else((prot.pos == 248) , "red", "black")) %>% 
  mutate(hotspot_color = if_else(is.na(hotspot_color), "black", hotspot_color))

vaf_comparison_prep <- vaf_comparison_prep %>%
  mutate(
    hotspot_status = case_when(
      protein_variant == "R248Q" ~ "R248Q",
      protein_variant == "R248W" ~ "R248W",
      TRUE ~ "Non-248 mutation"
    )
  ) %>%
  filter(Tissue != "Skin")

vaf_duplex_reads <- ggplot(data=vaf_comparison_prep, aes(x= t_alt_count, y= VAF, color=hotspot_status)) +
  geom_point(size = 2) +
  scale_y_log10() +
  scale_x_log10() +
  #scale_color_identity() +
  scale_color_manual(
    values = c(
      "R248Q" = "#E69F00",         # red
      "R248W" = "#56B4E9",         # blue
      "Non-248 mutation" = "black" # default
    ),
    name = "Mutation type"
  ) +
  labs(
    x = "Duplex reads",
    y = "Variant allele frequency (VAF)"
  ) +
  theme_minimal(base_size = 8) +        # overall font size
  theme(
    legend.position = c(0.3, 0.9),
    legend.title = element_blank(),
    legend.text  = element_text(size = 8),
    axis.title   = element_text(size = 8),
    axis.text    = element_text(size = 8)
  )
  
vaf_duplex_reads

ggsave("results/vaf_duplex_reads_supp.png", vaf_duplex_reads, width = 3, height = 4, units = "in", dpi = 300)

