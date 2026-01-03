
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/output/SBS/li_fraumeni.SBS6.region"
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"

df <- read_delim(sigPath, delim = "\t")
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID, -subject)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)

SBS_6_counts <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")
  #filter(count > 0)
SBS_6_counts %>% print(n=Inf)

SBS_6_counts <- left_join(SBS_6_counts, sampleID, by = "sample")
SBS_6_counts
custom_colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
                   "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")

custom_colors <- c("#E6194B", "#A40000", "#FF5733", "#C70039",  # 4 blood related
                   "#FFE119", "#911EB4", "#F032E6", "#3CB44B", "#9A6324",
                   "#004C99", "#0086C3", "#009688", "#33A1C9", "#006D77", "#00A6A6")  # 6 cancer

manual_tissue_order <- c("Whole blood", "Buffy coat", "Plasma", "Bone marrow", "Colon", "Skin", "Skin, non-sun-exposed", "Liver", "Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis")

SBS_6_counts$tissue <- factor(SBS_6_counts$tissue, levels = manual_tissue_order)

mutation_type_plot <- ggplot(SBS_6_counts, aes(x = MutationType, y = Count, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_plot)
mutation_type_plot_scaled <- ggplot(SBS_6_counts, aes(x = MutationType, y = Count, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Count") +
  ylim(0,200) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_plot_scaled)

path1 <- "results/SBS6_proportion_patient.png"
path2 <- "results/SBS6_proportion_patient_scaled.png"
ggsave(path1, mutation_type_plot, width=10, height=4)
ggsave(path2, mutation_type_plot_scaled, width=10, height=4)


SBS_6_counts_proportion <- SBS_6_counts %>%
  group_by(tissue) %>%
  mutate(proportion = Count/(sum(Count)))

# custom_colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
#                    "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
#                    "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")

#custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
#custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#DDCC77", "#CC6677", "#882255", "#DDDDDD")
custom_colors <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")

mutation_type_proportion_plot <- ggplot(SBS_6_counts_proportion, aes(x = tissue, y = proportion, fill = MutationType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Proportion", title = "Mutagenesis mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_proportion_plot)

path1 <- "results/SBS6_stacked_proportion_patient.png"
ggsave(path1, mutation_type_proportion_plot, width=8, height=4)
