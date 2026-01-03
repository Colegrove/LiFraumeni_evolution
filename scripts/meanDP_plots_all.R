## Li-Fraumeni meanDP vs number of unique mutations and total mutant reads

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


family <- c("Family member A", "Family member B", "Family member C")
family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
DP_summary %>% print(n=Inf)
DP_summary_coding <- DP_summary %>% 
  filter(coding == "coding")


all_samples <-
  filt_maf %>%
  #filter(Subject %in% family) %>%
  filter(coding == "coding") %>%
  filter(Tissue != "Urine cells")
  #filter(Subject == "Patient")

unique_p53_all <- all_samples %>% 
  group_by(Subject, Samp) %>% 
  summarize(unique_p53_count = n(), .groups = "drop") %>% 
  complete(
    Subject = family,
    #Samp = unique(family_pbmc$Samp),
    fill = list(unique_p53_count = 0)
  ) %>%
  print(n=Inf)

p53_reads_all <- all_samples %>%
  group_by(Subject, Tissue, Samp) %>%
  summarize(total_p53_reads = sum(t_alt_count), .groups = 'drop') %>%
  complete(
    Subject = family,
    fill = list(total_p53_reads = 0)
  )
p53_reads_all

unique_p53_all
mutations_depth <- right_join(DP_summary_coding, unique_p53_all, by = "Samp")
mutations_depth
mutations_depth <- right_join(mutations_depth, p53_reads_all, by = c("Samp", "Subject"))
mutations_depth <- mutations_depth %>%
  mutate(background = ifelse(Subject == "Patient", TRUE, FALSE))


###############################################################################
############### depth summary of samples
###############################################################################

mutations_depth <- mutations_depth %>% 
  mutate(x_label = if_else(Tissue == "PBMC", Subject, Tissue)) %>%
  print(n=Inf)

ggplot(mutations_depth, aes(x = x_label, y = meanDP)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Tissue / Subject",
    y = "Average Sequencing Depth"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10)
  )

ggplot(mutations_depth, aes(x = reorder(x_label, meanDP), y = meanDP)) +
  geom_segment(aes(xend = x_label, y = 0, yend = meanDP), color = "gray70") +
  geom_point(size = 4, color = "steelblue") +
  coord_flip() +
  labs(x = "Sample", y = "Mean Depth") +
  theme_minimal()

ggplot(mutations_depth, aes(x = reorder(x_label, meanDP), y = meanDP)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(meanDP, 1)), hjust = -0.1, size = 3) +
  coord_flip() +
  theme_minimal() +
  ylim(0,32000) +
  labs(x = "Sample", y = "Mean Sequencing Depth")

###############################################################################
############### depth and tp53 mutations
###############################################################################

#linear_unique <- lm(unique_p53_count ~ meanDP, data=mutations_depth)
#summary(linear_unique)
unique_plot <- ggplot() +
  geom_point(
    data = filter(mutations_depth, background == TRUE),
    #aes(x = meanDP, y = unique_p53_count, color = Subject),
    aes(x = meanDP, y = unique_p53_count),
    size = 2,
  ) +
  geom_point(
    data = filter(mutations_depth, background == FALSE),
    #aes(x = meanDP, y = unique_p53_count, color = Subject),
    aes(x = meanDP, y = unique_p53_count),
    size = 2, 
  ) +
  labs(
    x = "Mean sequencing depth",
    y = "Unique TP53 mutation count"#,
    #color = "Subject"
  ) +
  #ylim(0,150) +
  ylim(0,150) +
  coord_cartesian(xlim = c(0, 30000)) + 
  theme_minimal() + 
  theme(axis.line.x=element_line(color="black"), 
        axis.line.y=element_line(color="black"), 
        axis.text = element_text(size=12)) 

unique_plot
#linear_unique <- lm(total_p53_reads ~ meanDP, data=mutations_depth)
#summary(linear_unique)
total_reads_plot <- ggplot() +
  geom_point(
    data = filter(mutations_depth, background == TRUE),
    #aes(x = meanDP, y = total_p53_reads, color = Subject),
    aes(x = meanDP, y = total_p53_reads),
    size = 2,
  ) +
  geom_point(
    data = filter(mutations_depth, background == FALSE),
    #aes(x = meanDP, y = total_p53_reads, color = Subject),
    aes(x = meanDP, y = total_p53_reads),
    size = 2, 
  ) +
  labs(
    x = "Mean sequencing depth",
    y = "Unique TP53 mutation count"#,
    #color = "Subject"
  ) +
  
  #geom_point(shape=1, color='black', size = 3) +
  labs(
    x = "Mean sequencing depth",
    y = "Total TP53 mutation read count",
    color = "Subject"
  ) +
  #ylim(0,150) +
  ylim(0,150) + 
  coord_cartesian(xlim = c(0, 30000)) +
  theme_minimal() + 
  theme(axis.line.x=element_line(color="black"), 
        axis.line.y=element_line(color="black"),
        axis.text = element_text(size=12)) 

total_reads_plot

shared_legend <- get_legend(
  unique_plot + 
    theme(legend.position = "bottom")
)

plot1_no_legend <- unique_plot + theme(legend.position = "none")
plot2_no_legend <- total_reads_plot + theme(legend.position = "none")

combined_plots <- plot_grid(plot1_no_legend, plot2_no_legend, ncol = 2)
final_plot <- plot_grid(
  combined_plots, 
  #shared_legend, 
  ncol = 1,
  rel_heights = c(1, 0.5)  # Adjust the legend space
)

final_plot

ggsave("results/DP_comparisons_all_participants.png", final_plot, width = 10, height = 5, units = "in")


#ggsave("results/skyscraper_plot_slide.png", skyscraper_slide, width = 17, height = 13, units = "in", dpi = 300)
