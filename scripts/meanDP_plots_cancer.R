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

cancer_samples = c("All tissue-types",
                   "Mediastinal metastasis",
                   "Lung metastasis",
                   "Esophageal cancer 1",
                   "Esophageal cancer 2",
                   "Liver metastasis 1",
                   "Liver metastasis 2")

DP_summary_coding <- DP_summary %>% 
  filter(coding == "coding")

patient_cancer <-
  filt_maf %>%
  filter(Tissue %in% cancer_samples) %>%
  filter(coding == "coding") %>%
  filter(Subject == "Patient")

unique_p53_cancer <- patient_cancer %>% 
  group_by(Tissue, Samp) %>% 
  summarize(unique_p53_count = n()) %>% 
  print(n=Inf)

p53_reads_cancer <- patient_cancer %>%
  group_by(Tissue) %>%
  summarize(total_p53_reads = sum(t_alt_count))

mutations_depth <- right_join(DP_summary_coding, unique_p53_cancer, by = "Samp")
mutations_depth <- right_join(mutations_depth, p53_reads_cancer, by = "Tissue")
mutations_depth %>% print(n=Inf)


#linear_unique <- lm(unique_p53_count ~ meanDP, data=mutations_depth)
#summary(linear_unique)
unique_plot <- ggplot(mutations_depth, aes(x = meanDP, y = unique_p53_count, color = Tissue)) +
  geom_point(size = 3) +
  #geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    x = "Mean DP",
    y = "Unique TP53 mutation count",
    color = "Tissue"
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
total_reads_plot <- ggplot(mutations_depth, aes(x = meanDP, y = total_p53_reads, color = Tissue)) +
  geom_point(size = 3) +
  labs(
    x = "Mean DP",
    y = "Total TP53 mutation read count",
    color = "Tissue"
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

ggsave("results/DP_comparisons_patient_cancer.png", final_plot, width = 10, height = 5, units = "in")


#ggsave("results/skyscraper_plot_slide.png", skyscraper_slide, width = 17, height = 13, units = "in", dpi = 300)
