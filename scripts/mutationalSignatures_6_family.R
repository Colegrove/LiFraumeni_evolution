
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_family_controls/output/SBS/li_fraumeni.SBS6.region"
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"

df <- read_delim(sigPath, delim = "\t")
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID, -tissue)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)
sampleID
SBS_6_counts <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")
  #filter(count > 0)
SBS_6_counts %>% print(n=Inf)

SBS_6_counts <- left_join(SBS_6_counts, sampleID, by = "sample")
SBS_6_counts
custom_colors <- c("#000000", "#C70039", "#E6194B", "#FFA07A",
                  #"#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
                   "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")


SBS_6_counts$subject
manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")

SBS_6_counts$subject <- factor(SBS_6_counts$subject, levels = manual_subject_order)


mutation_type_plot <- ggplot(SBS_6_counts, aes(x = MutationType, y = Count, fill = subject)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_plot)

path1 <- "results/SBS6_proportion_family_controls.png"
ggsave(path1, mutation_type_plot, width=10, height=4)



SBS_6_counts_proportion <- SBS_6_counts %>%
  group_by(subject) %>%
  mutate(proportion = Count/(sum(Count)))

# custom_colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
#                    "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
#                    "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")

# custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
# custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#DDCC77", "#CC6677", "#882255", "#DDDDDD")
custom_colors <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
mutation_type_proportion_plot <- ggplot(SBS_6_counts_proportion, aes(x = subject, y = proportion, fill = MutationType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Proportion", title = "Mutagenesis mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_proportion_plot)

path1 <- "results/SBS6_stacked_proportion_family_controls.png"
ggsave(path1, mutation_type_proportion_plot, width=8, height=4)
