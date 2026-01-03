
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_TP53/output/SBS/li_fraumeni.SBS6.region"
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"

#####################################################################
######### SBS 6 Control Blood TP53
#####################################################################

df <- read_delim(sigPath, delim = "\t")
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)

SBS_6_counts_all <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")


SBS_6_counts_all <- left_join(SBS_6_counts_all, sampleID, by = "sample")
SBS_6_counts_control <- SBS_6_counts_control_all %>%
  filter(tissue == "PBMC" | tissue == "Buffy coat")
SBS_6_counts_control %>% filter(subject == "Patient:")

custom_colors <- c("#000000", "#C70039", "#E6194B", "#FFA07A",
                  #"#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
                   "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")


##### TP53 proband and Family member and control 

SBS_6_counts_control$subject
manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")

SBS_6_counts_control$subject <- factor(SBS_6_counts_control$subject, levels = manual_subject_order)

mutation_type_plot <- ggplot(SBS_6_counts_control, aes(x = MutationType, y = Count, fill = subject)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_plot)

path1 <- "results/SBS6_proportion_TP53_family_control.png"
ggsave(path1, mutation_type_plot, width=10, height=4)


SBS_6_counts_proportion <- SBS_6_counts_control %>%
  group_by(subject) %>%
  mutate(proportion = Count/(sum(Count)))

# custom_colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
#                    "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
#                    "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")

custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
custom_colors <- c("#332288", "#88CCEE", "#44AA99", "#DDCC77", "#CC6677", "#882255", "#DDDDDD")
custom_colors <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")


mutation_type_proportion_plot <- ggplot(SBS_6_counts_proportion, aes(x = subject, y = proportion, fill = MutationType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Proportion", title = "TP53 mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_proportion_plot)

path1 <- "results/SBS6_stacked_proportion_TP53_family_controls.png"
ggsave(path1, mutation_type_proportion_plot, width=8, height=4)



#####################################################################
######### SBS 6 Patient Tissue TP53
#####################################################################


SBS_6_counts_patient <- SBS_6_counts_all %>% filter(subject == "Patient:")

custom_colors <- c("#E6194B", "#A40000", "#FF5733", "#C70039",  # 4 blood related
                   "#FFE119", "#911EB4", "#F032E6", "#3CB44B", "#9A6324",
                   "#004C99", "#0086C3", "#009688", "#33A1C9", "#006D77", "#00A6A6")  # 6 cancer
manual_tissue_order <- c("Whole blood", "Buffy coat", "Plasma", "Bone marrow", "Colon", "Skin", "Skin, non-sun-exposed", "Liver", "Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis")
SBS_6_counts_patient$tissue <- factor(SBS_6_counts_patient$tissue, levels = manual_tissue_order)

mutation_type_plot_patient_tissue <- ggplot(SBS_6_counts_patient, aes(x = MutationType, y = Count, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Count") +
  theme_minimal() +
  #ylim(0,20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_plot_patient_tissue)

path2 <- "results/SBS6_proportion_TP53_patient_unscaled.png"
ggsave(path2, mutation_type_plot_patient_tissue, width=10, height=4)

## scale to remove Skin C>T
mutation_type_plot_patient_tissue <- ggplot(SBS_6_counts_patient, aes(x = MutationType, y = Count, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Mutation Type", y = "Count") +
  theme_minimal() +
  ylim(0,20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_plot_patient_tissue)

path2 <- "results/SBS6_proportion_TP53_patient_scaled.png"
ggsave(path2, mutation_type_plot_patient_tissue, width=10, height=4)


SBS_6_counts_proportion <- SBS_6_counts_patient %>%
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
  labs(x = "Mutation Type", y = "Proportion", title = "TP53 mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
show(mutation_type_proportion_plot)

path1 <- "results/SBS6_stacked_proportion_TP53_patient.png"
ggsave(path1, mutation_type_proportion_plot, width=8, height=4)











#####################################################################
######### SBS 6 Patient TP53 vs Mutagenesis Percent Change
#####################################################################

# load patient mutagenesis
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/output/SBS/li_fraumeni.SBS6.region"
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"

df <- read_delim(sigPath, delim = "\t")
df
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID, -subject)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)

SBS_6_counts_mutagenesis <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")
#filter(count > 0)
SBS_6_counts_mutagenesis %>% print(n=Inf)

SBS_6_counts_mutagenesis <- left_join(SBS_6_counts_mutagenesis, sampleID, by = "sample")
SBS_6_counts_mutagenesis %>% print(n=Inf)

SBS_6_counts_mutagenesis <- SBS_6_counts_mutagenesis %>%
  mutate(panel = "MUT")
SBS_6_counts_mutagenesis <- left_join(SBS_6_counts_mutagenesis, sampleID, by = c("sample", "tissue"))
SBS_6_counts_mutagenesis

# load patient TP53 (from above module)

SBS_6_counts_patient_TP53 <- SBS_6_counts_patient %>%
  mutate(panel = "TP53") %>%
  dplyr::select(-subject)
SBS_6_counts_patient_TP53

df_combined <- rbind(SBS_6_counts_mutagenesis, SBS_6_counts_patient_TP53)
df_combined

df_combined <- df_combined %>%
  group_by(panel, tissue) %>%
  mutate(TotalCount = sum(Count)) %>%
  ungroup()

# Compute the ratio for each MutationType
df_combined <- df_combined %>%
  mutate(Ratio = Count / TotalCount)

df_wide <- df_combined %>%
  select(MutationType, tissue, panel, Ratio) %>%
  pivot_wider(names_from = panel, values_from = Ratio) %>%
  mutate(Percent_change = ((TP53 - MUT) / MUT) * 100)
df_wide

custom_colors <- c("#E6194B", "#A40000", "#FF5733", "#C70039",  # 4 blood related
                   "#FFE119", "#911EB4", "#F032E6", "#3CB44B", "#9A6324",
                   "#004C99", "#0086C3", "#009688", "#33A1C9", "#006D77", "#00A6A6")  # 6 cancer
manual_tissue_order <- c("Whole blood", "Buffy coat", "Plasma", "Bone marrow", "Colon", "Skin", "Skin, non-sun-exposed", "Liver", "Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis")
df_wide$tissue <- factor(df_wide$tissue, levels = manual_tissue_order)


percent_change_patient <- ggplot(df_wide, aes(x = MutationType, y = Percent_change, fill = MutationType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ tissue) +  # Separate plots for each tissue
  labs(title = "Percentage Change in Mutation Ratios",
       x = "Mutation Type", y = "Percentage Change (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
percent_change_patient

path2 <- "results/SBS6_percent_change_patient.png"
ggsave(path2, percent_change_patient, width=12, height=6)

#### by mutation
percent_change_patient <- ggplot(df_wide, aes(x = tissue, y = Percent_change, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ MutationType) +  # Separate plots for each tissue
  labs(title = "Percentage Change in Mutation Ratios",
       x = "Mutation Type", y = "Percentage Change (%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
show(percent_change_patient)

path2 <- "results/SBS6_percent_change_patient.png"
ggsave(path2, percent_change_patient, width=12, height=6)


#####################################################################
######### SBS 6 Control TP53 vs Mutagenesis Percent change
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

## load control and family TP53

sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_TP53/output/SBS/li_fraumeni.SBS6.region"
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
df <- read_delim(sigPath, delim = "\t")
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)

SBS_6_counts_family_controls_TP53 <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")


SBS_6_counts_family_controls_TP53 <- SBS_6_counts_family_controls_TP53 %>%
  mutate(panel = "TP53") 
SBS_6_counts_family_controls_TP53

SBS_6_counts_family_controls_TP53 <- left_join(SBS_6_counts_family_controls_TP53, sampleID, by = "sample")
SBS_6_counts_family_controls_TP53 <- SBS_6_counts_family_controls_TP53 %>%
  filter(tissue == "PBMC" | tissue == "Buffy coat")
df_combined <- rbind(SBS_6_counts_mutagenesis, SBS_6_counts_family_controls_TP53)
df_combined

df_combined <- df_combined %>%
  group_by(panel, subject) %>%
  mutate(TotalCount = sum(Count)) %>%
  ungroup()
df_combined
# Compute the ratio for each MutationType
df_combined <- df_combined %>%
  mutate(Ratio = Count / TotalCount)
df_combined %>% print(n=Inf)

df_wide <- df_combined %>%
  select(MutationType, subject, panel, Ratio) %>%
  pivot_wider(names_from = panel, values_from = Ratio) %>%
  mutate(Percent_change = ((TP53 - MUT) / MUT) * 100)
df_wide %>% print(n=Inf)


percent_change_controls <- ggplot(df_wide, aes(x = MutationType, y = Percent_change, fill = MutationType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ subject) +  # Separate plots for each tissue
  labs(title = "Percentage Change in Mutation Ratios",
       x = "Mutation Type", y = "Percentage Change (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
show(percent_change_controls)

path2 <- "results/SBS6_percent_change_family_controls.png"
ggsave(path2, percent_change_controls, width=12, height=6)


manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")

df_wide$subject <- factor(df_wide$subject, levels = manual_subject_order)

custom_colors <- c("#000000", "#C70039", "#E6194B", "#FFA07A",
                   #"#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
                   "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
                   "#008080", "#E6BEFF", "#9A6324", "#800000", "#808000")

percent_change_controls <- ggplot(df_wide, aes(x = subject, y = Percent_change, fill = subject)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ MutationType) +  # Separate plots for each tissue
  labs(title = "Percentage Change in Mutation Ratios",
       x = "Mutation Type", y = "Percentage Change (%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
show(percent_change_controls)
df_wide
path2 <- "results/SBS6_percent_change_family_controls.png"
ggsave(path2, percent_change_controls, width=12, height=6)
