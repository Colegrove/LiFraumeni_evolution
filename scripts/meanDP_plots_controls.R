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

DP_summary_coding <- DP_summary %>% 
  filter(coding == "coding")

control_pbmc <-
  filt_maf %>%
  filter(!(Subject %in% family)) %>%
  filter(coding == "coding") %>%
  filter(Tissue == "PBMC")

unique_p53_control_pbmc <- control_pbmc %>% 
  group_by(Subject, Samp) %>% 
  summarize(unique_p53_count = n(), .groups = "drop") %>% 
  # complete(
  #   Subject = family,
  #   #Samp = unique(family_pbmc$Samp),
  #   fill = list(unique_p53_count = 0)
  # ) %>%
  print(n=Inf)

p53_reads_control_pbmc <- control_pbmc %>%
  group_by(Subject) %>%
  summarize(total_p53_reads = sum(t_alt_count), .groups = 'drop') #%>%
  # complete(
  #   Subject = family,
  #   fill = list(total_p53_reads = 0)
  # )


mutations_depth <- right_join(DP_summary_coding, unique_p53_control_pbmc, by = "Samp")
mutations_depth <- right_join(mutations_depth, p53_reads_control_pbmc, by = "Subject")
mutations_depth %>% print(n=Inf)


#linear_unique <- lm(unique_p53_count ~ meanDP, data=mutations_depth)
#summary(linear_unique)
unique_plot <- ggplot(mutations_depth, aes(x = meanDP, y = unique_p53_count, color = Subject)) +
  geom_point(size = 3) +
  #geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    x = "Mean DP",
    y = "Unique TP53 mutation count",
    color = "Subject"
  ) +
  #ylim(0,150) +
  ylim(0,150) +
  coord_cartesian(xlim = c(0, 30000)) + 
  theme_minimal() + 
  theme(axis.line.x=element_line(color="black"), 
        axis.line.y=element_line(color="black")) 

unique_plot
#linear_unique <- lm(total_p53_reads ~ meanDP, data=mutations_depth)
#summary(linear_unique)
total_reads_plot <- ggplot(mutations_depth, aes(x = meanDP, y = total_p53_reads, color = Subject)) +
  geom_point(size = 3) +
  labs(
    x = "Mean DP",
    y = "Total TP53 mutation read count",
    color = "Subject"
  ) +
  #ylim(0,150) +
  ylim(0,150) + 
  coord_cartesian(xlim = c(0, 30000)) +
  theme_minimal() + 
  theme(axis.line.x=element_line(color="black"), 
        axis.line.y=element_line(color="black")) 

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
  shared_legend, 
  ncol = 1,
  rel_heights = c(1, 0.5)  # Adjust the legend space
)

final_plot

ggsave("results/DP_comparisons_controls.png", final_plot, width = 10, height = 5, units = "in")


#ggsave("results/skyscraper_plot_slide.png", skyscraper_slide, width = 17, height = 13, units = "in", dpi = 300)
