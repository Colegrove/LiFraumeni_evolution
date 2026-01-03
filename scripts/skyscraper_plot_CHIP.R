## Li-Fraumeni chip panel skyscraper plot

family <- c("Family member A", "Family member B", "Family member C")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")

skyscraper_prep %>% filter(Subject == "UW volunteer 7") %>% filter(coding_from_bed == "coding" & coding_from_maf == "coding") %>% print(n=20, width=Inf)
filt_maf_CHIP %>% filter(Subject == "UW volunteer 7") %>% mutate(SampCodingOrder = rank(desc(t_alt_count), ties.method = "first")) %>% print(n=500, width=Inf)
filt_maf_CHIP %>% filter(Subject == "UW volunteer 7") %>% print(n=20, width=Inf)

filt_maf_CHIP %>% filter(coding_from_bed == "non-coding") %>% filter(Subject == "Patient" & Tissue == "Buffy coat")  %>% filter(VAF <0.3) %>% filter(Hugo_Symbol != "Unknown") %>% print(width=Inf)

filt_maf_CHIP %>% filter(VAF <0.3) %>% print(width = Inf)


skyscraper_prep <-
  filt_maf_CHIP %>%
  filter(coding_from_bed == "coding") %>%
  #filter(coding_from_bed == "coding" & coding_from_maf == "coding") %>%
  #filter(coding_from_bed == "coding" & coding_from_maf == "non-coding") %>%
  filter(VAF < 0.3) %>% 
  group_by(Subject, Tissue) %>%
  mutate(SampCodingOrder = rank(desc(t_alt_count), ties.method = "first")) %>%
  ungroup() %>%
  filter(Tissue == "PBMC" | Tissue == "Buffy coat") #%>%
  # 
  # mutate(
  #   color_group = case_when(
  #     am_class == "likely_benign" & t_alt_count == 1 ~ "likely_benign",
  #     am_class == "likely_benign" & t_alt_count > 1 ~ "likely_benign_LC",
  #     am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
  #     am_class == "likely_pathogenic" & t_alt_count > 1 ~ "likely_pathogenic_LC",
  #     am_class == "ambiguous" & t_alt_count == 1 ~ "ambiguous",
  #     am_class == "ambiguous" & t_alt_count > 1 ~ "ambiguous_LC"
  #   )
  # )

skyscraper_prep %>% filter(Mutation_type == "Splice") %>% print(width)

skyscraper <- skyscraper_prep %>%
  ggplot(
    aes(
      x = Subject,
      y = SampCodingOrder,
      label = t_alt_count
    )
  ) +
  geom_tile(fill = "#DDDDDD") +
  
  # uncomment
  geom_text( size=font.subscript.size*25.4/72.27) +
  
  geom_tile(data = skyscraper_prep,
            aes(
              x = Subject,
              y = SampCodingOrder
            ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
  scale_y_continuous(expand = c(0,0),
                     breaks=c(25,75,125,175,225),
                     labels=c(0,50,100,150, 200),
                     #limits = c(0,250),
                     limits = c(0,20),
                     name = "Number of CHIP coding mutations") +
  
  #Col.amClass.fill +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size = 8),
        panel.spacing = unit(0,"pt"),
        #strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.x = element_blank()
        #legend.position = c(0.25,0.7),
        #legend.position = element_blank(),
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
  ) + geom_text(size=10*25.4/72.27)

show(skyscraper_slide)
ggsave("results/skyscraper_plot_slide_CHIP_introns.png", skyscraper_slide, width = 16, height = 9, units = "in", dpi = 300)


#########

library(ggforce)

skyscraper_violin <- skyscraper_prep %>%
  ggplot(aes(x = Subject, y = t_alt_count)) +
  
  # Violin plot to show distribution of mutation counts
  geom_violin(fill = "#DDDDDD", color = "black", alpha = 0.7) +
  
  # Jittered points to show individual mutations
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
  facet_zoom(ylim = c(0,150)) +
  
  # scale_y_continuous(expand = c(0, 0),
  #                    name = "TP53 Mutation Count") +
  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
        panel.spacing = unit(0, "pt"),
        strip.text.x.top = element_markdown(family = "sans", size = 8, face = "bold", color = "black"),
        axis.title.y = element_markdown(family = "sans", size = 8, face = "bold", color = "black"),
        axis.title.x = element_blank())

# Show the plot
show(skyscraper_violin)


#####################################################################
# heatmaps
#####################################################################


## heatmap by gene and unique mutation
heatmap_data <- skyscraper_prep %>%
  group_by(Subject, Hugo_Symbol) %>%
  filter(Hugo_Symbol != "Unknown") %>%
  summarise(Unique_Mutations = n_distinct(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2), 
            .groups = "drop") %>%
  group_by(Hugo_Symbol) %>%
  mutate(total_mut = sum(Unique_Mutations)) %>%
  ungroup() %>%
  mutate(Hugo_Symbol = fct_reorder(Hugo_Symbol, total_mut)) %>%
  mutate(Subject = factor(
    Subject, levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
  ))


heatmap_plot <- ggplot(heatmap_data, aes(x = Subject, y = Hugo_Symbol, fill = Unique_Mutations)) +
  
  geom_tile(color = "white") +  # Add white borders between tiles
  geom_text(aes(label = Unique_Mutations), color = "black", size = 3, fontface = "bold") +
  
  scale_fill_gradient(low = "white", high = "red", na.value = "white", oob = scales::squish) +  # Heatmap color scale
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(0, "pt"),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  
  labs(fill = "Unique Mutations")

# Show the heatmap
show(heatmap_plot)

subjects <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
ages <- c("20-30", "20-30", "20-30", "20-30", "30s", "39", "61", "60s", "69", "70s", "70s")
LFS <- c("-", "-", "-", "-", "+", "+", "-","-", "+","-", "-")
CTx <- c("-", "-", "-", "-", "+", "-", "-","-", "-","-", "+")

age_lf <- data.frame(
  x_axis_var = subjects, 
  CTx = CTx,
  LFS = LFS,
  age = ages
)

age_lf_long <- age_lf %>%
  pivot_longer(cols = c(age, LFS, CTx), names_to = "variable", values_to = "value") %>%
  mutate(y = case_when(
    variable == "age" ~ 2,
    variable == "LFS" ~ 1.5,
    variable == "CTx" ~ 1, 
  )) %>% 
  mutate(x_axis_var = factor(x_axis_var, levels = unique(x_axis_var)),
         x_axis_var = factor(x_axis_var, levels = unique(x_axis_var))  # preserve original order
  )


annotation_plot <- ggplot(age_lf_long, aes(x = x_axis_var, y = y)) +
  geom_text(aes(label = value), size = 4, vjust = 0) +
  scale_y_continuous(
    breaks = c(1.1, 1.6, 2.1),
    labels = c("CTx", "LFS", "age"),
    #expand = c(0, 0),
    limits = c(1, 2.3)
  ) +
  scale_x_discrete(expand = c(0.07, 0)) +
  theme_minimal(base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
show(annotation_plot)


combined_plot <- cowplot::plot_grid(
  heatmap_plot + theme(plot.margin = margin(0, 0, 2, 0, "pt"), legend.position = "top"),
  annotation_plot,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE
)
show(combined_plot)
ggsave("results/CHIP_gene_heatmap_unique_muts.png", combined_plot, width = 6, height = 6, units = "in", dpi = 300)



## heatmap by gene and total reads
heatmap_data <- skyscraper_prep %>%
  group_by(Subject, Hugo_Symbol) %>%
  filter(Hugo_Symbol != "Unknown") %>%
  summarise(Total_t_alt_count = sum(t_alt_count, na.rm = TRUE), .groups = "drop")   %>%
  group_by(Hugo_Symbol) %>%
  mutate(sum_count = sum(Total_t_alt_count)) %>%
  ungroup() %>%
  mutate(Hugo_Symbol = fct_reorder(Hugo_Symbol, sum_count)) %>%
  mutate(Subject = factor(
    Subject, levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
  ))

heatmap_plot <- ggplot(heatmap_data, aes(x = Subject, y = Hugo_Symbol, fill = Total_t_alt_count)) +
  
  geom_tile(color = "white") +  # Add white borders between tiles
  geom_text(aes(label = Total_t_alt_count), color = "black", size = 3, fontface = "bold") +
  
  scale_fill_gradient(low = "white", high = "red", na.value = "white", limits=c(0,400), oob = scales::squish) +  # Heatmap color scale
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(0, "pt"),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  
  labs(fill = "Total alt count")

# Show the heatmap
show(heatmap_plot)


combined_plot <- cowplot::plot_grid(
  heatmap_plot + theme(plot.margin = margin(0, 0, 2, 0, "pt"), legend.position = "top"),
  annotation_plot,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE
)
show(combined_plot)
ggsave("results/CHIP_gene_heatmap_total_count.png", combined_plot, width = 6, height = 6, units = "in", dpi = 300)


MAF_table %>% filter(prot.pos == 248) %>% filter(Hugo_Symbol == "TP53") %>% print(n=Inf, width = Inf)
