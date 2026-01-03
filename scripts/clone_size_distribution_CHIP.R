## CHIP panel clone size distribution

family <- c("Family member A", "Family member B", "Family member C")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")

#######################################################
##### Number of mutations per subject
#######################################################


mutationCounts <- 
  filt_maf_CHIP %>%
  filter(coding_from_bed == "coding") %>%
  filter(coding_from_bed == "coding" & coding_from_maf == "coding") %>%
  #filter(coding_from_bed == "coding" & coding_from_maf == "non-coding") %>%
  filter(VAF < 0.3) %>% 
  group_by(Subject, Tissue) %>%
  mutate(SampCodingOrder = rank(desc(t_alt_count), ties.method = "first")) %>%
  ungroup() %>%
  filter(Tissue == "PBMC" | Tissue == "Buffy coat") %>%
  group_by(Subject, Hugo_Symbol) %>%
  summarise(total_mut_count = n(), .groups = "drop") %>%
  mutate(Subject = factor(
  Subject, levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
  ))


gene_colors = c("#800000", # Maroon
                            "#9A6324", # Brown
                            "#808000", # Olive
                            "#469990", # Teal
                            "#000075", # Navy
                            "#e6194B", # Red
                            "#f58231", # Orange
                            "#ffe119", # Yellow
                            "#bfef45", # Lime
                            "#3cb44b", # Green
                            "#42d4f4", # Cyan
                            "#4363d8", # Blue
                            "#911eb4", # Purple
                            "#f032e6", # Magenta
                            "#a9a9a9", # Grey
                            "#fabed4", # Pink
                            "#ffd8b1", # Apricot
                            "#fffac8", # Beige
                            "#aaffc3", # Mint
                            "#dcbeff", # Lavender
                            "#000000", 
                            "#DDDDDD", 
                            "#AA4499", 
                            "#44aa99", 
                            "#b0e0e6", # Powder Blue
                            "#ffb6c1", # Light Pink
                            "#dda0dd", # Plum
                            "#add8e6"  # Light Blue
                            )
mutationCounts
mutationCounts %>% filter(Subject == "Patient")
mutationCount_plot <- mutationCounts %>% 
  ggplot(aes(x = Subject,y = total_mut_count, fill = Hugo_Symbol)) + 
  geom_col() +
  theme_bw()  +
  scale_fill_manual(values = gene_colors) +
  scale_y_continuous(limits = c(0,250)) +
  labs(y = "Number of CHIP mutations", fill = "Gene") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size = 8), 
        axis.title.x = element_blank())
show(mutationCount_plot)

##### Add age and LFS information
subjects <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Buffy coat", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
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

legend <- cowplot::get_legend(
  mutationCount_plot + theme(legend.position = "right")
)
legend_aligned <- cowplot::plot_grid(
  legend,
  NULL,  # spacer to push legend upward
  ncol = 1,
  rel_heights = c(1, .5)  # stretch legend space
)

mutationCount_plot_nolegend <- mutationCount_plot +
  theme(legend.position = "none")

combined_plot <- cowplot::plot_grid(
  mutationCount_plot_nolegend + theme(plot.margin = margin(0, 0, 2, 0, "pt")),
  annotation_plot,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE, 
  align = "v", 
  axis = "lr"
)

combined_plot_legend <- cowplot::plot_grid(
  combined_plot, 
  legend_aligned,
  ncol = 2, 
  rel_widths = c(1, 0.5), 
  align = "h",
  axis = "tb"
  
)

show(combined_plot_legend)

#ggsave("results/CHIP_numberMutations_annotated.png", combined_plot, width = 7, height = 7, units = "in", dpi = 300)
ggsave("results/CHIP_numberMutations_annotated.png", combined_plot_legend, width = 9, height = 5, units = "in", dpi = 300)



#######################################################
##### Clone size distributions color by age
#######################################################

cloneSize_prep <-
  filt_maf_CHIP %>%
  filter(coding_from_bed == "coding") %>%
  filter(coding_from_bed == "coding" & coding_from_maf == "coding") %>%
  #filter(coding_from_bed == "coding" & coding_from_maf == "non-coding") %>%
  filter(VAF < 0.3) %>% 
  group_by(Subject, Tissue) %>%
  mutate(SampCodingOrder = rank(desc(t_alt_count), ties.method = "first")) %>%
  ungroup() %>%
  filter(Tissue == "PBMC" | Tissue == "Buffy coat")

cloneSizes = cloneSize_prep %>%
  group_by(Subject) %>%
  summarize(
    count = n(),
    t_alt_counts = list(t_alt_count),
    clone_size_counts = list({
      counts <- t_alt_count
      ks <- min(counts):max(counts)
      names(ks) <- ks
      sapply(ks, function(k) sum(counts >= k))
    })
  )
cloneSizes$clone_size_counts

ge_long <- cloneSizes %>%
  dplyr::select(Subject, clone_size_counts) %>%
  unnest_longer(clone_size_counts, indices_to = "k") %>%
  rename(count_clone_size = clone_size_counts) %>%
  mutate(k = as.integer(k))

age_df <- age_lf %>%
  mutate(age = if_else(age == "20-30", "25", age)) %>%
  mutate(age = if_else(age == "30s", "35", age)) %>%
  mutate(age = if_else(age == "60s", "65", age)) %>%
  mutate(age = if_else(age == "70s", "75", age)) %>%
  mutate(Subject = if_else(x_axis_var == "Buffy coat", "Patient", x_axis_var)) %>%
  dplyr::select(Subject,age)

ge_long <- ge_long %>%
  left_join(age_df)

custom_colors <- c("#FFA07A","#E6194B","#C70039", 
                   #"#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#deebf7", "#000000", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594")

custom_colors_age <- c(
  "#f0f921","#fcce25","#fca636","#f2844b", "#e16462",
  "#cc4778","#b12a90","#8f0da4","#6a00a8","#41049d","#0d0887")
custom_colors_age <- c(
  "#fcce25","#fcce25","#fcce25","#fcce25","#fca636",
  "#f2844b","#b12a90","#8f0da4","#6a00a8","#0d0887","#0d0887")
custom_colors_age <- c(
  "#fcce25","#fcce25","#fcce25","#fcce25","cyan", #"#fca636",
  "#f2844b","#b12a90","#8f0da4","#6a00a8","#0d0887","#0d0887")





#manual_subject_order <- c("Patient", "Family member A", "Family member B", "Family member C", "UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")
#manual_subject_order <- c("Family member A", "Family member B", "Family member C", "UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7", "Patient")
manual_subject_order <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7", "Patient")
legend_subject_order <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")


ge_long <- ge_long %>%
  mutate(Subject = factor(Subject, levels = manual_subject_order)) %>%
  mutate(subject_age_label = paste0(Subject, " (", age, ")")) 

legend_levels <- ge_long %>%
  distinct(Subject, subject_age_label) %>%
  mutate(Subject = factor(Subject, levels = legend_subject_order)) %>%
  arrange(Subject) %>%
  pull(subject_age_label)

legend_values <- ge_long %>%
  distinct(Subject, subject_age_label) %>%
  mutate(Subject = factor(Subject, levels = manual_subject_order)) %>%
  arrange(Subject) %>%
  pull(subject_age_label)

ge_long <- ge_long %>%
  mutate(subject_age_label = factor(subject_age_label, levels = legend_levels)) %>%
  arrange(if_else(as.character(Subject) == "Patient", 1, 0))  # Patient drawn last

ge_long <- ge_long %>%
  mutate(age_numeric = as.numeric(age))

clone_size_plot <- ggplot() +
  geom_line(
    data = ge_long %>% filter(Subject != "Patient"),
    aes(x = k, y = count_clone_size, color = subject_age_label),
    size = 1.2, 
    linetype = 12
  ) +
  geom_line(
    data = ge_long %>% filter(Subject == "Patient"),
    aes(x = k, y = count_clone_size, color = subject_age_label),
    size = 1.2
  ) +
  labs(
    x = "clone size (x)",
    y = "Number of clones ≥ x"
  ) +  
  #scale_color_viridis_c(option = "plasma", direction = -1) +
  scale_color_manual(values = custom_colors_age, breaks = legend_levels) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
                      labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(#limits = c(1,1000),
                #breaks = 10^seq(0,3),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal() + 
  theme(legend.position = c(0.8, 0.6), 
        legend.title = element_blank())
show(clone_size_plot)
ggsave("results/clone_size_distribution_blood_age.png", clone_size_plot, width = 8, height = 4, units = "in", dpi = 300)


##############################################################################
##### Clone size distributions color by LFS status
##############################################################################

cloneSize_prep <-
  filt_maf_CHIP %>%
  filter(coding_from_bed == "coding") %>%
  filter(coding_from_bed == "coding" & coding_from_maf == "coding") %>%
  #filter(coding_from_bed == "coding" & coding_from_maf == "non-coding") %>%
  filter(VAF < 0.3) %>% 
  group_by(Subject, Tissue) %>%
  mutate(SampCodingOrder = rank(desc(t_alt_count), ties.method = "first")) %>%
  ungroup() %>%
  filter(Tissue == "PBMC" | Tissue == "Buffy coat")

cloneSizes = cloneSize_prep %>%
  group_by(Subject) %>%
  summarize(
    count = n(),
    t_alt_counts = list(t_alt_count),
    clone_size_counts = list({
      counts <- t_alt_count
      ks <- min(counts):max(counts)
      names(ks) <- ks
      sapply(ks, function(k) sum(counts >= k))
    })
  )
cloneSizes$clone_size_counts

ge_long <- cloneSizes %>%
  dplyr::select(Subject, clone_size_counts) %>%
  unnest_longer(clone_size_counts, indices_to = "k") %>%
  rename(count_clone_size = clone_size_counts) %>%
  mutate(k = as.integer(k))

age_df <- age_lf %>%
  mutate(age = if_else(age == "20-30", "25", age)) %>%
  mutate(age = if_else(age == "30s", "35", age)) %>%
  mutate(age = if_else(age == "60s", "65", age)) %>%
  mutate(age = if_else(age == "70s", "75", age)) %>%
  mutate(Subject = if_else(x_axis_var == "Buffy coat", "Patient", x_axis_var)) %>%
  dplyr::select(Subject,age)

ge_long <- ge_long %>%
  left_join(age_df)

custom_colors <- c("#FFA07A","#E6194B","#C70039", 
                   #"#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#deebf7", "#000000", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594")

custom_colors_age <- c(
  "#f0f921","#fcce25","#fca636","#f2844b", "#e16462",
  "#cc4778","#b12a90","#8f0da4","#6a00a8","#41049d","#0d0887")
custom_colors_age <- c(
  "#fcce25","#fcce25","#fcce25","#fcce25","#fca636",
  "#f2844b","#b12a90","#8f0da4","#6a00a8","#0d0887","#0d0887")
custom_colors_age <- c(
  "#fcce25","#fcce25","#fcce25","#fcce25","cyan", #"#fca636",
  "#f2844b","#b12a90","#8f0da4","#6a00a8","#0d0887","#0d0887")

subject_colors_LFS <- c("Patient" = "red", 
                    "Family member A" = "red", 
                    "Family member B" = "blue", 
                    "Family member C" = "red", 
                    "UW volunteer 1" = "blue", 
                    "UW volunteer 2" = "blue", 
                    "UW volunteer 3" = "blue", 
                    "UW volunteer 4" = "blue", 
                    "UW volunteer 5" = "blue", 
                    "UW volunteer 6" = "blue", 
                    "UW volunteer 7" = "blue")




#manual_subject_order <- c("Patient", "Family member A", "Family member B", "Family member C", "UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")
#manual_subject_order <- c("Family member A", "Family member B", "Family member C", "UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7", "Patient")
manual_subject_order <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7", "Patient")
legend_subject_order <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")


ge_long <- ge_long %>%
  mutate(Subject = factor(Subject, levels = manual_subject_order)) %>%
  mutate(subject_age_label = paste0(Subject, " (", age, ")")) 

legend_levels <- ge_long %>%
  distinct(Subject, subject_age_label) %>%
  mutate(Subject = factor(Subject, levels = legend_subject_order)) %>%
  arrange(Subject) %>%
  pull(subject_age_label)

legend_values <- ge_long %>%
  distinct(Subject, subject_age_label) %>%
  mutate(Subject = factor(Subject, levels = manual_subject_order)) %>%
  arrange(Subject) %>%
  pull(subject_age_label)

ge_long <- ge_long %>%
  mutate(subject_age_label = factor(subject_age_label, levels = legend_levels)) %>%
  arrange(if_else(as.character(Subject) == "Patient", 1, 0))  # Patient drawn last

ge_long <- ge_long %>%
  mutate(age_numeric = as.numeric(age))

ge_long$LFS_status_color <- ifelse(ge_long$Subject %in% c("Patient", "Family member A", "Family member C"), 
                                              "red", "blue")

clone_size_plot <- ggplot() + 
  geom_line(
    data = ge_long %>% filter(Subject != "Patient"),
    aes(x = k, y = count_clone_size, color = Subject),  # Use custom_color
    size = 1.2,
    linetype = 1
  ) + 
  geom_line(
    data = ge_long %>% filter(Subject == "Patient"),
    aes(x = k, y = count_clone_size, color = Subject),  # Use custom_color
    size = 1.2
  ) + 
  labs(
    x = "Clone Size (x)",
    y = "Number of Clones ≥ x"
  ) + 
  scale_color_manual(values = subject_colors_LFS) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x), 
    labels = trans_format("log10", math_format(10^.x))
  ) + 
  theme_minimal() + 
  theme(legend.position = c(0.8, 0.6), 
        legend.title = element_blank())

show(clone_size_plot)
ggsave("results/clone_size_distribution_blood_age.png", clone_size_plot, width = 8, height = 4, units = "in", dpi = 300)




#######################################################################
########### Synonymous vs non-synonymous mutations
#######################################################################




syn_nonsyn <-
  filt_maf_CHIP %>%
  filter(coding_from_bed == "coding") %>%
  #filter(coding_from_bed == "coding" & coding_from_maf == "coding") %>%
  #filter(coding_from_bed == "coding" & coding_from_maf == "non-coding") %>%
  filter(VAF < 0.4) %>% 
  group_by(Subject, Tissue) %>%
  mutate(SampCodingOrder = rank(desc(t_alt_count), ties.method = "first")) %>%
  ungroup() %>%
  filter(Tissue == "PBMC" | Tissue == "Buffy coat") %>%
  mutate(Subject = factor(Subject, levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
  ))

custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
#custom_colors <- c("#332288", "#44AA99", "#DDCC77", "#88CCEE", "#882255", "#AA4499", "#DDDDDD")

mutation_types <- ggplot(syn_nonsyn, aes(x = Subject, fill = Mutation_type)) +
  geom_bar(position = "fill") +
  labs(
    x = "Subject",
    y = "Count",
    fill = "Mutation_type"
  ) +
  scale_fill_manual(values = c(custom_colors)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

syn_nonsyn %>% filter(Mutation_type == "Splice") %>% print(width = Inf)

subjects <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Buffy coat", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
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
  mutate(x_axis_var = factor(x_axis_var, levels = unique(subjects)),
         x_axis_var = factor(x_axis_var, levels = unique(subjects))  # preserve original order
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
  mutation_types + theme(plot.margin = margin(0, 0, 2, 0, "pt"), legend.position = "top"),
  annotation_plot,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE
)
show(combined_plot)
ggsave("results/mutation_type_synonymous_CHIP.png", combined_plot, width = 7, height = 4, units = "in", dpi = 300)
