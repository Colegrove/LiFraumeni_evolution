
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_CHIP/output/SBS/li_fraumeni.SBS6.region"
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"


#####################################################################
######### SBS 6 CHIP 
#####################################################################

df <- read_delim(sigPath, delim = "\t")
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)

SBS_6_counts_CHIP <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")
  #filter(count > 0)

SBS_6_counts_CHIP <- left_join(SBS_6_counts_CHIP, sampleID, by = "sample")

custom_colors <- c("#000000", "#C70039", "#E6194B", "#FFA07A",
                   #"#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
                   "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")

manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")

SBS_6_counts_CHIP$subject <- factor(SBS_6_counts_CHIP$subject, levels = manual_subject_order)

SBS_6_counts_CHIP <- SBS_6_counts_CHIP %>%
  filter(tissue == "PBMC" | tissue == "Buffy coat") 

# mutation_type_plot <- ggplot(SBS_6_counts_CHIP, aes(x = MutationType, y = Count, fill = subject)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_fill_manual(values = custom_colors) +
#   labs(x = "Mutation Type", y = "Count") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# show(mutation_type_plot)

# path1 <- "results/SBS6_counts_CHIP.png"
# ggsave(path1, mutation_type_plot, width=10, height=4)


SBS_6_counts_proportion <- SBS_6_counts_CHIP %>%
  group_by(subject) %>%
  mutate(proportion = Count/(sum(Count))) %>%
  mutate(subject = factor(
    subject, levels = c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Patient:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")
  ))

SBS_6_counts_proportion

#custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
#custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#DDCC77", "#CC6677", "#882255", "#DDDDDD")
custom_colors <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")

mutation_type_proportion_plot <- ggplot(SBS_6_counts_proportion, aes(x = subject, y = proportion, fill = MutationType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "CHIP mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
show(mutation_type_proportion_plot)

path1 <- "results/SBS6_stacked_proportion_CHIP.png"
ggsave(path1, mutation_type_proportion_plot, width=8, height=4)

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

combined_plot <- cowplot::plot_grid(
  mutation_type_proportion_plot + theme(plot.margin = margin(0, 0, 2, 0, "pt"), legend.position = "top") + guides(fill = guide_legend(nrow = 1)),
  annotation_plot,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE
)
show(combined_plot)

path1 <- "results/SBS6_stacked_proportion_CHIP.png"
ggsave(path1, combined_plot, width=8, height=4)




#####################################################################
######### SBS 6 Control CHIP vs Mutagenesis percent changes
#####################################################################

# load family and control mutagenesis

sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_family_controls/output/SBS/li_fraumeni.SBS6.region"
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"

df <- read_delim(sigPath, delim = "\t")
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)

SBS_6_counts_all_control_mut <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")

SBS_6_counts_all_control_mut <- left_join(SBS_6_counts_all_control_mut, sampleID, by = "sample")
SBS_6_counts_control_mut <- SBS_6_counts_all_control_mut %>%
  filter(tissue == "PBMC" | tissue == "Buffy coat")
SBS_6_counts_control_mut %>% filter(subject == "Patient:")
SBS_6_counts_mutagenesis <- SBS_6_counts_control_mut %>%
  mutate(panel = "MUT")

SBS_6_counts_mutagenesis
## load control and family CHIP from above

SBS_6_counts_CHIP_panel <- SBS_6_counts_CHIP %>%
  mutate(panel = "CHIP")

df_combined <- rbind(SBS_6_counts_mutagenesis, SBS_6_counts_CHIP_panel)
df_combined %>% print(n=Inf)

df_combined <- df_combined %>%
  group_by(panel, subject) %>%
  mutate(TotalCount = sum(Count)) %>%
  ungroup()
df_combined %>% print(n=Inf)
# Compute the ratio for each MutationType
df_combined <- df_combined %>%
  mutate(Ratio = Count / TotalCount)
df_combined %>% print(n=Inf)

df_wide <- df_combined %>%
  select(MutationType, subject, panel, Ratio) %>%
  pivot_wider(names_from = panel, values_from = Ratio) %>%
  mutate(Percent_change = ((CHIP - MUT) / MUT) * 100)
df_wide %>% print(n=Inf)


percent_change_controls <- ggplot(df_wide, aes(x = MutationType, y = Percent_change, fill = MutationType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ subject) +  # Separate plots for each tissue
  labs(title = "Percentage Change in Mutation Ratios",
       x = "Mutation Type", y = "Percentage Change (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
show(percent_change_controls)

path2 <- "results/SBS6_percent_change_CHIP.png"
ggsave(path2, percent_change_controls, width=12, height=6)


custom_colors <- c("#000000", "#C70039", "#E6194B", "#FFA07A",
                   #"#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
                   "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")

manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")

df_wide$subject <- factor(df_wide$subject, levels = manual_subject_order)


percent_change_controls <- ggplot(df_wide, aes(x = subject, y = Percent_change, fill = subject)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ MutationType) +  # Separate plots for each tissue
  labs(title = "Percentage Change in Mutation Ratios",
       x = "Mutation Type", y = "Percentage Change (%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
show(percent_change_controls)

path2 <- "results/SBS6_percent_change_CHIP.png"
ggsave(path2, percent_change_controls, width=12, height=6)
